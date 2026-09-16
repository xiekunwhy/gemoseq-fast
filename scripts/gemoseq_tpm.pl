#!/usr/bin/env perl
#=============================================================================
# gemoseq_tpm.pl -- 依据 GeMoSeq 输出的 GFF3 计算并添加 TPM 属性
#
# 用法:
#   perl gemoseq_tpm.pl merged.gff3 > merged.tpm.gff3
#   cat chr01.gff3 chr02.gff3 ... | perl gemoseq_tpm.pl - > all.tpm.gff3
#   perl gemoseq_tpm.pl -o out.gff3 chr01.gff3 chr02.gff3
#   perl gemoseq_tpm.pl -f score merged.gff3   # 旧口径：用 score（EM 分配的 reads 数）
#
# 原理（默认 -f avgcov，与 StringTie 同构的覆盖度口径）:
#   TPM = avgCov / sum(avgCov) * 1e6
#   avgCov 是转录本的平均每碱基覆盖度（已含长度归一），来自剪接图节点计数，
#   不受 ra=true（丰度还原）影响——热点区饱和而非爆炸，分布稳健。
#
#   -f score 旧口径: RPK = score / (外显子总长 kb)，TPM = RPK/sum(RPK)*1e6。
#   注意: 若 gemoseq 用 ra=true 跑过，极端热点位点的 score 会被 1/sampleProb
#   放大到天文数字并压扁其他基因的 TPM——用 ra=true 跑的数据务必用默认口径。
#
# 注意:
#   1) avgCov（或 score）是必需属性；缺失时该转录本不加 TPM 并计数警告。
#   2) 分染色体单独跑 GeMoSeq 时，请使用相同的参数。
#   3) 基因/转录本 ID 以 (染色体, ID) 联合键控，分染色体结果 cat 到一起时
#      跨染色体的重复 ID（如 gnwc=false 时的 G1,G2,...）不会串。
#   4) TPM 追加在 mRNA 行的属性列末尾（TPM=<值>;，与 GeMoSeq 属性结尾带分号的风格一致），其余行原样输出。
#=============================================================================
use strict;
use warnings;

my $out = "";
my $field = "avgcov";
my @files;
while (@ARGV) {
	my $a = shift @ARGV;
	if ($a eq "-o") { $out = shift @ARGV; }
	elsif ($a eq "-f") { $field = lc(shift @ARGV); }
	elsif ($a eq "-h" || $a eq "--help") { exec "perldoc $0"; }
	else { push @files, $a; }
}
push @files, "-" if (!@files);
die "ERROR: -f 只接受 avgcov 或 score\n" if ($field ne "avgcov" && $field ne "score");

# ---- 第一遍：收集 mRNA 的 (seqid,ID) -> 度量值，exon 归属累长度 ----
my (%score, %exonlen, %order);
for my $f (@files) {
	my $fh;
	if ($f eq "-") { $fh = *STDIN; }
	else { open($fh, "<", $f) or die "ERROR: 无法打开 $f: $!\n"; }
	while (my $line = <$fh>) {
		next if ($line =~ /^#/ || $line !~ /\S/);
		chomp($line);
		my @F = split(/\t/, $line);
		next if (@F < 9);
		if ($F[2] eq "mRNA" || $F[2] eq "transcript") {
			my ($id) = $F[8] =~ /ID=([^;]+)/;
			next if (!defined $id);
			my ($sc);
			if ($field eq "avgcov") {
				($sc) = $F[8] =~ /avgCov=([\d.eE+-]+)/;
			} else {
				($sc) = $F[8] =~ /score=([\d.eE+-]+)/;
			}
			$score{$F[0]}{$id} = $sc;   # undef 表示缺失
			push @{$order{$F[0]}}, $id if (!exists $exonlen{$F[0]}{$id});
		} elsif ($F[2] eq "exon") {
			my ($par) = $F[8] =~ /Parent=([^;]+)/;
			next if (!defined $par);
			for my $p (split(/,/, $par)) {
				$exonlen{$F[0]}{$p} += $F[4] - $F[3] + 1;
			}
		}
	}
	close($fh) if ($f ne "-");
}

# ---- 计算 TPM 的分子（avgcov 已含长度归一，score 需除长度 kb）与全局和 ----
my (%tpm, $sum, $n_noscore, $n_nolen);
$sum = 0; $n_noscore = 0; $n_nolen = 0;
for my $chr (keys %score) {
	for my $id (keys %{$score{$chr}}) {
		my $sc = $score{$chr}{$id};
		if (!defined $sc) { $n_noscore++; next; }
		my $num;
		if ($field eq "avgcov") {
			$num = $sc;
		} else {
			my $len = $exonlen{$chr}{$id} || 0;
			if ($len <= 0) { $n_nolen++; next; }
			$num = $sc / ($len / 1000.0);
		}
		$tpm{$chr}{$id} = $num;
		$sum += $num;
	}
}
die "ERROR: 没有任何可用的 ($field) 信息，检查输入是否为 GeMoSeq 的 GFF3\n" if ($sum <= 0);

for my $chr (keys %tpm) {
	for my $id (keys %{$tpm{$chr}}) {
		$tpm{$chr}{$id} = $tpm{$chr}{$id} / $sum * 1e6;
	}
}

# ---- 第二遍：输出，mRNA 行追加 TPM ----
my $ofh = *STDOUT;
if ($out ne "") { open($ofh, ">", $out) or die "ERROR: 无法写入 $out: $!\n"; }
my $n_tagged = 0;
for my $f (@files) {
	my $fh;
	if ($f eq "-") { $fh = *STDIN; }
	else { open($fh, "<", $f) or die "ERROR: 无法打开 $f: $!\n"; }
	while (my $line = <$fh>) {
		if ($line =~ /^#/ || $line !~ /\S/) {
			print $ofh $line;
			next;
		}
		chomp($line);
		my @F = split(/\t/, $line);
		if (@F >= 9 && ($F[2] eq "mRNA" || $F[2] eq "transcript")) {
			my ($id) = $F[8] =~ /ID=([^;]+)/;
			if (defined $id && exists $tpm{$F[0]}{$id}) {
				$F[8] =~ s/;?$//;
				$F[8] .= ";TPM=" . sprintf("%.3f", $tpm{$F[0]}{$id}) . ";";
				$n_tagged++;
			}
			print $ofh join("\t", @F), "\n";
		} else {
			print $ofh $line, "\n";
		}
	}
	close($fh) if ($f ne "-");
}
close($ofh) if ($out ne "");

my $chk = 0;
$chk += $_ for (map { values %$_ } values %tpm);
printf STDERR "[gemoseq_tpm] 口径=%s; 已标记 %d 个转录本; TPM 总和 = %.1f (应≈1e6); 缺值跳过 %d, 缺 exon 跳过 %d\n",
	$field, $n_tagged, $chk, $n_noscore, $n_nolen;
