# GeMoSeq-fast

> English | [中文](README_cn.md)

A high-depth optimized fork of [GeMoSeq](https://www.jstacs.de/index.php/GeMoSeq)
(the `projects/gemoseq` sub-project of [Jstacs/Jstacs](https://github.com/Jstacs/Jstacs),
formerly named GeMoRNA), built for **transcript reconstruction from multi-sample merged
BAM files** (transcript evidence production for de novo genome annotation).

The original tool's memory grows with depth × region length on such inputs and gains
nothing beyond 6 threads; this fork makes memory and runtime **essentially independent of
sequencing depth**, and adds the engineering features large-genome work requires:
fragment collapsing, an asynchronous I/O pipeline, sparse EM quantification,
reference-restricted runs, CSI index support (including a critical >512 Mb fix),
prefix-based output naming, abundance rescaling with TPM computation, and reusable
global read statistics.

- Upstream version: GeMoSeq 1.2.3 (2025-11-11)
- This fork's version: `1.2.3-fast`
- License: GPL v3 (inherited from upstream, see [LICENSE](LICENSE))

---

## 1. Motivation

For annotation-oriented transcript calling, merging the alignments of multiple RNA-seq
samples into one BAM before assembly yields clearly better gene completeness (BUSCO)
than assembling each sample separately and merging GTFs afterwards. But merged-BAM
depth is the sum of all samples, and the original GeMoSeq on such input:

1. stores all reads of a genomic region as a linked list of `SAMRecord` objects, capped
   at `maxcov × region length` (default 100 × 750 kb = 75M reads/region), while spike
   coverage (rRNA, chloroplast) bypasses down-sampling entirely;
2. processes every read base-by-base **three times** (graph building, read assignment,
   EM quantification), so runtime scales linearly with depth;
3. ingests BAM (decompress + decode + region building) on a single thread — no gain
   beyond 6 threads, I/O is the bottleneck;
4. has engineering gaps: no CSI index support, no per-chromosome targeting, no control
   over output names, down-sampling stdout spam, whole-genome memory loading, and
   several crashes (including upstream issue #75).

## 2. Features and modifications

### 2.1 Fragment-structure collapsing (ReadGroup) — the key to depth insensitivity

Reads are converted to a compact representation at ingestion (exon-block arrays + gap
types + pre-computed mismatch count) and collapsed by a signature of "full alignment
structure of both mates + strand": **N fragments with identical structure are stored
once with weight N**. Graph counts use `nReads += weight` (exactly equal to per-read
increments); the quantification matrix has one row per group initialized with the group
weight — mathematically equivalent to the original "one row per read, then makeUnique
collapsing". Mates are paired by read name with a region-local buffer, preserving the
original "a pair counts as one observation" semantics.

### 2.2 Down-sampling and memory caps

- the original `mrc` (average region coverage, default 100) down-sampling semantics are
  kept, plus a new **absolute cap** `mrpr` (default 4,000,000 reads/region) which also
  tames spike coverage;
- group-level down-sampling (weighted fragments kept or dropped as a whole,
  expectation-preserving);
- removed the original `#N->M` stdout spam on every down-sampling event.

### 2.3 Asynchronous I/O pipeline and sparse computation

- **background ReadStats**: insert/intron-length statistics overlap with the main pass;
- **async decode + parallel conversion**: a producer thread does BGZF inflate/SAMRecord
  decode, 4 converter threads do cigar/block/mismatch/strand analysis in parallel with
  ordered, back-pressured handoff (order preserved, results deterministic);
- **sparse EM quantification**: the EM loop iterates per-read sparse compatible-candidate
  lists instead of the dense reads×transcripts matrix — O(Σ compat) per iteration,
  identical result;
- `ReadGraph.remove` avoids `LinkedList.removeAll` O(n×m) (now O(n+m));
- `Node` edge maps are created lazily (~112B → ~50B per base position);
- `nextSplit()` rewritten as a single BFS with arrays sized to each component's actual
  span (eliminates O(#components²));
- **lazy genome loading via `.fai`**: chromosomes are read from the fasta only when
  first needed;
- threads 8-12 are the practical sweet spot; for more parallelism, run one process per
  chromosome.

### 2.4 Concurrency and crash fixes

- fixed a multi-threading race: in stranded mode, forward/reverse Regions were read and
  written concurrently by different workers, causing `ConcurrentModificationException`
  (latent in the original as well);
- **fixed the empty-subregion crash** in coverage-based region splitting
  ([upstream issue #75](https://github.com/Jstacs/Jstacs/issues/75)): intervals fully
  inside a long intron receive no reads (every read spans the split boundaries),
  producing empty sub-regions that crashed the worker with a `NullPointerException`;
  empty sub-regions are now dropped and the compute path guards against them;
- unmapped reads (no reference) are skipped safely.
- **fixed a hang at multi-reference rollover**: htsjdk permits only one open query iterator per reader; restricted runs with more than one reference threw `IllegalStateException: Iteration in progress` on the producer thread, and non-daemon worker threads then kept the JVM alive forever (sleep state). Iterators are now closed between references and pool workers are daemon threads (errors exit cleanly). `readstats` also accepts `r=` (single reference) with the same short names as gemoseq.


### 2.5 htsjdk upgrade, CSI and index handling

- the bundled htsjdk was upgraded from 2.5.0-SNAPSHOT (2016, no CSI support) to **2.24.1**;
- automatic index detection: `*.bam.bai` / `*.bam.csi` / `*.bai` / `*.csi` naming
  conventions, `.bai` preferred;
- **a prominent warning is printed when the index is older (mtime) than the BAM** — a
  stale index silently drops data (measured: a stale index returned only 47% of reads);
- with no index at all, falls back to streaming the full BAM and filtering by reference
  name (correct but slower, with a notice);
- **patched `CSIIndex.class`** (source `src/htsjdk/samtools/CSIIndex.java`): upstream
  htsjdk 2.24 computes a `minimumOffset` for whole-reference CSI queries by walking up
  the bin hierarchy; on depth-6 CSIs written by samtools/htslib 1.20+ (records beyond
  512 Mb are placed in level-1 bins), the walk anchored on those tail-holding level-1
  bins and `Chunk.optimizeChunkList` silently **discarded all chunks before
  ~1.07-1.61 Gb** — restricted runs (`r`/`rl`) returned only the chromosome tail. The
  patch skips the walk when `startPos <= 0` (whole-reference queries), restoring full
  coverage.

### 2.6 Restricted runs and outputs

- `r=<chromosome>` / `rl=<list file>`: process only the given references (indexed query);
- `o=<prefix>`: outputs `<prefix>.Transcript_Predictions.gff3` and
  `<prefix>.protocol_gemorna.txt` (inside outdir); the intermediate predictions file is
  also prefixed (`<prefix>.predictions.tmp`);
- `ra=true`: rescales abundances in down-sampled regions by 1/sampling probability
  (restores `score` to the original read-count scale; read the TPM note in §5 before
  using it for quantification).

### 2.7 Statistics reuse: the `readstats` tool and the `rs` parameter

The jar bundles a **`readstats`** tool (peer of gemoseq/predictCDS/GAF/Analyzer/merge):
it computes intron-length statistics (meanSplit/sdSplit/meanReadLen) over the whole BAM
(or a restricted reference set) and writes them to a text file. Passing that file to
GeMoSeq via `rs=<file>` **skips the internal statistics pass** and uses these values
for splice-graph pruning — compute once globally, reuse everywhere, for both speed and
a unified pruning criterion.

## 3. Parameters

All original parameters are unchanged and compatible
(`g/m/s/l/sil/lr/mnor/mfor/mnoir/mfoir/p/mrpg/mrpt/pa/sf/mrl/mrc/mfgl/q/mpl/gp/gnwc`).
New/changed:

| Parameter | short | description | default |
|---|---|---|---|
| Restrict to reference | `r` | process only this reference from the BAM (requires an index) | off |
| Reference list | `rl` | file with one reference per line; process only those into one GFF (requires an index) | off |
| Stream full BAM | `sfb` | when restricting, stream the whole BAM and filter by name instead of using the index (slower but immune to index problems) | false |
| Output prefix | `o` | prefix output and intermediate files | off |
| Collapse identical fragments | `c` | fragment collapsing on/off (for A/B comparison) | true |
| Maximum reads per region | `mrpr` | absolute cap of reads kept per region | 4,000,000 |
| Rescale abundance | `ra` | rescale abundance by 1/down-sampling probability (restore original read-count scale) | false |
| Read statistics | `rs` | precomputed read statistics file (from `readstats`); skips the internal statistics pass | off |

`readstats` tool parameters: `m=<bam>` (required), `o=<output>` (default
readstats.stats), `Shortest intron length` (same as gemoseq's `sil`), `rl=<list>`
(optional restriction). **Output is written as `<outdir>/<o>` and the protocol as
`<outdir>/<o>.protocol`** (same prefix semantics as gemoseq).

## 4. Performance and correctness

### 4.1 Synthetic data (2 Mb genome, 330 genes, with a 100 kb hotspot; reproducible via bench/)

| dataset | original | this fork | output agreement |
|---|---|---|---|
| low (30×/3000×, 1.04M reads) | 18 s / 5.3 GB | 17 s / **2.4 GB** | **coordinates 100% identical** |
| hi (2000×/100,000×, 37.6M reads) | 157 s / 8.8 GB | **76 s / 4.3 GB** | 97.3% exact-coordinate match* |

\* "exact-coordinate match" is a strict metric: an mRNA counts as matching only if its
(chromosome, start, end, strand) tuple is identical on both sides — partial overlaps do
not count. All differences are in the hotspot where both versions apply random
down-sampling (the original itself down-samples randomly with a fixed seed, so any two
implementations fluctuate there); intron-chain agreement 96.9-97.6%.

### 4.2 Real data (plant genome, 1.1 Gb, 15 references; FR_SECOND_STRAND, threads=6)

| run | peak memory | wall time | notes |
|---|---|---|---|
| original, Chr01-only BAM (625 MB) | 11.5 GB | 143 s | 6765 transcripts |
| **this fork**, same input | **8.0 GB** | **83-101 s** | 6710; intron-chain agreement with the original 99.5% |
| this fork, whole-genome BAM (4.8 GB) + `r=Chr01` | 7.5 GB | 108 s | 6710, coordinate-identical to the split-BAM run |
| this fork, whole-genome BAM + `rl=ptg.list o=ptg` | 0.55 GB | 234 s | small-contig list mode |
| this fork, whole-genome BAM, full run + `o=cca` | 11.4 GB | 635 s | 48219 transcripts |

Thread scaling: 6→12 threads gains another 15-20%; **repeated runs on identical input
are byte-identical** (fixed random seed, and all random decisions happen on the single
ingestion thread).

### 4.3 Chromosomes beyond 512 Mb (depth-6 CSI)

Verified locally with a 600 Mb chromosome and a hand-written CSI built with
htslib-1.22 binning rules: before the patch, queries returned only the chromosome tail
(4960/9920); after the patch, `query == stream` (9920=9920) and restricted runs recover
both head and tail genes.

## 5. Typical usage

```bash
# 0) optional but recommended: compute statistics once over the whole BAM
java -jar GeMoSeq-1.2.3-fast.jar readstats m=merged.bam o=genome.stats

# 1) whole genome in one go
java -Xmx16g -jar GeMoSeq-1.2.3-fast.jar gemoseq g=genome.fa m=merged.bam s=FR_SECOND_STRAND threads=6 o=cca

# 2) per-chromosome parallel runs (recommended for large genomes; rs reuses global stats)
java -Xmx16g -jar GeMoSeq-1.2.3-fast.jar gemoseq g=genome.fa m=merged.bam s=FR_SECOND_STRAND threads=6 r=Chr01 o=Chr01 rs=genome.stats

# 3) a set of small contigs in one run (shared statistics, avoids small-sample degeneracy)
java -Xmx16g -jar GeMoSeq-1.2.3-fast.jar gemoseq g=genome.fa m=merged.bam s=FR_SECOND_STRAND threads=6 rl=small_contigs.list o=small rs=genome.stats

# 4) merge per-chromosome results and compute TPM
cat chr*.Transcript_Predictions.gff3 | perl scripts/gemoseq_tpm.pl - > all.tpm.gff3
```

To squeeze memory further: `mrpr=2000000 threads=4`. For unbiased abundances at heavily
down-sampled loci: `ra=true`.

### TPM computation (`scripts/gemoseq_tpm.pl`)

Adds a `TPM=` attribute to every mRNA line; accepts multiple files or `-` (stdin),
`-o <file>` for output.

- **Default mode (`-f avgcov`, recommended)**: TPM is derived from the per-transcript
  `avgCov` attribute (mean exon coverage, already length-normalized):
  `TPM = avgCov / ΣavgCov × 1e6`. Same coverage-based notion as StringTie's TPM.
- **Legacy mode (`-f score`)**: TPM from `score` per exon-kb. Only meaningful on output
  produced with `ra=false`. **Do not use it on `ra=true` output**: hot-spot scores are
  rescaled by up to 1/1e-6, so a handful of transcripts swallow almost all TPM mass and
  everything else is flattened towards 0 (observed on real data: median TPM ≈ 0).
  With `ra=true` data, always use the default avgCov mode (avgCov is built from the
  splice-graph node counts and saturates at hotspots instead of exploding).

## 6. Building from source

Requirements: the official GeMoSeq-1.2.3.jar
([jstacs.de](https://www.jstacs.de/index.php/GeMoSeq)), htsjdk-2.24.1.jar (Maven
Central), JDK 11+, Python 3.

```bash
# Windows uses ';' as classpath separator; Linux ':'
javac -encoding UTF-8 -cp GeMoSeq-1.2.3.jar;htsjdk-2.24.1.jar -d out-fast src/gemoseq/*.java src/gemoma/*.java src/htsjdk/samtools/CSIIndex.java
python3 build/repack_jar.py --base GeMoSeq-1.2.3.jar --htsjdk htsjdk-2.24.1.jar --classes out-fast --out GeMoSeq-1.2.3-fast.jar
```

## 7. Verification

`bench/` contains the synthetic data generator and benchmark scripts (`GenTestData.java`,
`bench.sh`, `bench2.ps1`, `ReproNPE.java`, etc.) reproducing all numbers in section 4.
Acceptance criteria: 100% coordinate identity at non-downsampled loci, intron-chain
agreement ≥ 99% at downsampled loci.

## 8. Known limitations and semantic differences (stated honestly)

1. **Long-read mode** (`lr=true`) has been adapted but not thoroughly tested; long reads
   are structurally almost unique, so collapsing gains little there.
2. Down-sampling is done per group, which differs slightly from the original per-read
   coin flips — abundance values in downsampled regions may vary by a few percent
   (`ra=true` restores the expected scale).
3. With more than two alignments sharing a read name (secondary/supplementary), pair
   collapsing differs slightly from the original idMap semantics; such records are rare
   under the `q=40` MAPQ filter.
4. With `c=false` (collapsing off), output is coordinate-identical to the original
   (verified) — useful as a control.
5. Stale-index detection is based on file mtime; no content-level validation is
   performed (the warning advises re-indexing).
6. **ReadStats scope**: by default, restricted runs (`r`/`rl`) compute statistics from
   the restricted references only (so the output matches a split-BAM run exactly;
   whole-BAM statistics would differ slightly). For a unified criterion across runs —
   recommended for small contigs or when cross-chromosome consistency matters — compute
   whole-genome statistics with `readstats` and feed them via `rs=`.
7. A tiny contig run alone with very few spliced reads can get a degenerate (zero)
   intron-length variance and over-prune junctions — pool contigs via `rl=` or use
   `rs=` with global statistics.
8. `readstats` writes its output and protocol as `<outdir>/<o>` and
   `<outdir>/<o>.protocol`.

## 9. Repository layout

```
GeMoSeq-1.2.3-fast.jar   # ready-to-run fat jar (htsjdk 2.24.1 + CSIIndex patch included)
src/gemoseq/             # modified/added gemoseq sources (override same-named upstream classes)
src/gemoma/ReadStats.java  # ReadStats (with toFile/fromFile and restricted statistics)
src/htsjdk/samtools/CSIIndex.java  # the patched htsjdk class
scripts/gemoseq_tpm.pl   # TPM computation script
build/repack_jar.py      # repackaging script
bench/                   # synthetic data generator and benchmark scripts
LICENSE                  # GPL v3 (inherited from Jstacs)
```

## 10. Acknowledgements and citation

GeMoSeq is the transcript-reconstruction companion of GeMoMa, by Jens Keilwagen et al.;
see [upstream](https://github.com/Jstacs/Jstacs). This repository contains engineering
optimizations only and does not change the core algorithmic design; please cite upstream
as well when citing this work.
