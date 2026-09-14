#!/usr/bin/env perl
#=============================================================================
# gemoseq_tpm.pl -- 依据 GeMoSeq 输出的 GFF3 计算并添加 TPM 属性
#
# 用法:
#   perl gemoseq_tpm.pl merged.gff3 > merged.tpm.gff3
#   cat chr01.gff3 chr02.gff3 ... | perl gemoseq_tpm.pl - > all.tpm.gff3
#   perl gemoseq_tpm.pl -o out.gff3 chr01.gff3 chr02.gff3
#
# 原理:
#   TPM = RPK / sum(RPK) * 1e6,  其中 RPK = score / (转录本外显子总长 kb)
#   score 是 GeMoSeq 定量得到的转录本分配 reads 数（GFF3 mRNA 行的 score= 属性）。
#
# 注意:
#   1) score= 是必需属性；缺失时该转录本不加 TPM 并计数警告。
#   2) 分染色体单独跑 GeMoSeq 时，请使用相同的参数；若跑了降采样明显的深度，
#      建议跑 gemoseq 时加 ra=true（丰度按 1/采样概率还原），否则热点区 TPM 偏低。
#   3) 基因/转录本 ID 以 (染色体, ID) 联合键控，分染色体结果 cat 到一起时
#      跨染色体的重复 ID（如 gnwc=false 时的 G1,G2,...）不会串。
#   4) TPM 追加在 mRNA 行的属性列末尾（TPM=<值>），其余行原样输出。
#=============================================================================
use strict;
use warnings;

my $out = "";
my @files;
while (@ARGV) {
	my $a = shift @ARGV;
	if ($a eq "-o") { $out = shift @ARGV; }
	elsif ($a eq "-h" || $a eq "--help") { exec "perldoc $0"; }
	else { push @files, $a; }
}
push @files, "-" if (!@files);

# ---- 第一遍：收集 mRNA 的 (seqid,ID) -> score，exon 归属累长度 ----
my (%score, %exonlen, %order);
my $total_lines = 0;
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
			my ($sc) = $F[8] =~ /score=([\d.eE+-]+)/;
			$score{$F[0]}{$id} = $sc;   # undef 表示缺失
			push @{$order{$F[0]}}, $id if (!exists $exonlen{$F[0]}{$id});
		} elsif ($F[2] eq "exon") {
			my ($par) = $F[8] =~ /Parent=([^;]+)/;
			next if (!defined $par);
			for my $p (split(/,/, $par)) {
				$exonlen{$F[0]}{$p} += $F[4] - $F[3] + 1;
			}
		}
		$total_lines++;
	}
	close($fh) if ($f ne "-");
}

# ---- 计算 RPK 与全局和 ----
my (%tpm, $sum, $n_noscore, $n_nolen);
$sum = 0; $n_noscore = 0; $n_nolen = 0;
for my $chr (keys %score) {
	for my $id (keys %{$score{$chr}}) {
		my $sc = $score{$chr}{$id};
		if (!defined $sc) { $n_noscore++; next; }
		my $len = $exonlen{$chr}{$id} || 0;
		if ($len <= 0) { $n_nolen++; next; }
		my $rpk = $sc / ($len / 1000.0);
		$tpm{$chr}{$id} = $rpk;
		$sum += $rpk;
	}
}
die "ERROR: 没有任何可用 (score, exon) 信息，检查输入是否为 GeMoSeq 的 GFF3\n" if ($sum <= 0);

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
				$F[8] .= ";TPM=" . sprintf("%.3f", $tpm{$F[0]}{$id});
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
printf STDERR "[gemoseq_tpm] 已标记 %d 个转录本; TPM 总和 = %.1f (应≈1e6); 缺 score 跳过 %d, 缺 exon 跳过 %d\n",
	$n_tagged, $chk, $n_noscore, $n_nolen;
