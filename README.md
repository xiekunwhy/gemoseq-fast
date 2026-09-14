# GeMoSeq-fast

> English | [中文](README_cn.md)

A high-depth optimized fork of [GeMoSeq](https://www.jstacs.de/index.php/GeMoSeq)
(the `projects/gemoseq` sub-project of [Jstacs/Jstacs](https://github.com/Jstacs/Jstacs),
formerly named GeMoRNA).

Built for **transcript reconstruction from multi-sample merged BAM files** (transcript evidence
production for de novo genome annotation): the original tool's memory grows with
depth × region length on such inputs, while this fork makes memory and runtime
**essentially independent of sequencing depth**, and adds reference-restricted runs
(chromosome/scaffold/contig level), prefix-based output naming, and CSI index support.

- Upstream version: GeMoSeq 1.2.3 (2025-11-11)
- This fork's version: `1.2.3-fast`
- License: GPL v3 (inherited from upstream, see [LICENSE](LICENSE))

---

## 1. Motivation

For annotation-oriented transcript calling, merging the alignments of multiple RNA-seq samples
into one BAM before assembly yields clearly better gene completeness (BUSCO) than assembling each
sample separately and merging GTFs afterwards. But merged-BAM depth is the sum of all samples,
and the original GeMoSeq on such input:

1. stores all reads of a genomic region as a linked list of `SAMRecord` objects, capped at
   `maxcov × region length` — by default 100 × 750 kb = **75 million reads per region**,
   which explodes memory at high depth;
2. processes every read base-by-base **three times** (graph building, read assignment,
   EM quantification), so runtime scales linearly with depth;
3. `nextSplit()` allocates a full region-length node array and re-filters the edge list for
   every connected component (O(#components²));
4. down-sampling only triggers when the *average* region coverage exceeds the cap — **spike
   coverage (rRNA, chloroplast genes) bypasses it entirely**.

## 2. Core modifications

### 2.1 Fragment-structure collapsing (ReadGroup) — the key to depth insensitivity

Reads are converted to a compact representation at ingestion (exon-block arrays + gap types +
pre-computed mismatch count) and collapsed by a signature of "full alignment structure of both
mates + strand": **N fragments with identical structure are stored once with weight N**.

- Graph building: `nReads += weight` (exactly equal to per-read increments of the original);
- the quantification matrix has one row per group with `readWeights` initialized to the group
  weight — mathematically equivalent to the original "one row per read, then makeUnique
  collapsing with summed weights";
- mates are paired by read name with a region-local buffer, preserving the original
  "a pair counts as one observation" semantics;
- artificial `dummy` gap-fill reads keep the original behavior of sharing one index per region.

### 2.2 Down-sampling and memory caps

- the original `maxcov` (average region coverage) down-sampling semantics are kept, plus a new
  **absolute cap** `mrpr` (default 4,000,000 reads/region) which also tames spike coverage;
- group-level down-sampling (weighted fragments are kept or dropped as a whole);
- removed the original `#N->M` stdout spam on every down-sampling event.

### 2.3 Algorithm/data-structure fixes

- `nextSplit()` rewritten as a single BFS with arrays sized to the component's actual span
  (eliminates the O(#components²) behavior and repeated allocations);
- **fixed the empty-subregion crash** in coverage-based region splitting
  ([upstream issue #75](https://github.com/Jstacs/Jstacs/issues/75)): intervals fully inside a
  long intron receive no reads (every read spans the split boundaries), producing empty
  sub-regions that crashed the worker with a `NullPointerException`; empty sub-regions are now
  dropped and the compute path guards against them;
- `Node.addOutgoing` now uses a single hash lookup instead of containsKey+get+put;
- **lazy genome loading via `.fai`**: chromosomes are read from the fasta only when first needed,
  eliminating the whole-genome resident memory footprint;
- **fixed a multi-threading race**: in stranded mode, forward/reverse Regions were read and
  written concurrently by different workers, causing `ConcurrentModificationException`
  (latent in the original as well); key Region state access is now synchronized.

### 2.4 htsjdk upgrade + index handling

- the bundled htsjdk was upgraded from 2.5.0-SNAPSHOT (2016, no CSI support) to **2.24.1**;
- automatic index detection: `*.bam.bai` / `*.bam.csi` / `*.bai` / `*.csi` naming conventions,
  `.bai` preferred;
- **a prominent warning is printed when the index is older (mtime) than the BAM** — a stale
  index silently drops data;
- with no index at all, falls back to streaming the full BAM and filtering by reference name
  (correct but slower, with a notice);
- `ReadStats` (the insert/intron-length statistics used for splice-graph pruning) is also
  reference-aware: in restricted runs it reads only the requested references via indexed
  queries, so a restricted run matches a split-BAM run in both speed and output.

### 2.5 New parameters

| Parameter | short name | description | default |
|---|---|---|---|
| Restrict to reference | `r` | process only this reference (chromosome/scaffold/contig) from the BAM (requires an index) | off |
| Reference list | `rl` | file with one reference name per line; process only those, combined into one GFF (requires an index) | off |
| Output prefix | `o` | outputs are named `<prefix>.Transcript_Predictions.gff3` and `<prefix>.protocol_gemorna.txt` (inside outdir, default: current directory); the intermediate predictions file is also prefixed (`<prefix>.predictions.tmp`) | off |
| Collapse identical fragments | `c` | fragment collapsing on/off (for A/B comparison) | true |
| Maximum reads per region | `mrpr` | absolute cap of reads kept per region | 4,000,000 |

All original parameters are unchanged and fully compatible.

## 3. Performance and correctness

### 3.1 Synthetic data (2 Mb genome, 330 genes, with a 100 kb hotspot; reproducible via bench/)

| dataset | original | this fork | output agreement |
|---|---|---|---|
| low (30×/3000×, 1.04M reads) | 18 s / 5.3 GB | 17 s / **2.4 GB** | **coordinates 100% identical** |
| hi (2000×/100,000×, 37.6M reads) | 157 s / 8.8 GB | **76 s / 4.3 GB** | 97% coordinate overlap* |

\* "coordinate overlap" is a strict metric: an mRNA counts as matching only if its
(chromosome, start, end, strand) tuple is exactly identical on both sides — partial overlaps do
not count. All differences are in the hotspot where both versions apply random down-sampling
(the original itself down-samples randomly with a fixed seed, so any two implementations
fluctuate there).

### 3.2 Real data (plant genome, 1.1 Gb, 15 references; FR_SECOND_STRAND, threads=6)

| run | peak memory | wall time | notes |
|---|---|---|---|
| original, Chr01-only BAM (625 MB) | 11.5 GB | 143 s | 6765 transcripts |
| **this fork**, same input | **8.0 GB** | **101 s** | 6710 transcripts |
| this fork, whole-genome BAM (4.8 GB) + `r=Chr01` | 7.5 GB | 108 s | 6710, coordinate-identical to the split-BAM run |
| this fork, whole-genome BAM + `rl=ptg.list o=ptg` | 0.55 GB | 234 s | small-contig list mode |
| this fork, whole-genome BAM, full run + `o=cca` | 11.4 GB | 635 s | 48219 transcripts |

On identical data and parameters, intron-chain agreement between this fork and the original jar
is **99.5%** (differences concentrate on terminal exon boundaries at ultra-high-depth loci).

## 4. Usage

```bash
# whole genome
java -Xmx16g -jar GeMoSeq-1.2.3-fast.jar gemoseq g=genome.fa m=merged.bam s=FR_SECOND_STRAND threads=6 o=cca

# Chr01 only (indexed query)
java -Xmx16g -jar GeMoSeq-1.2.3-fast.jar gemoseq g=genome.fa m=merged.bam s=FR_SECOND_STRAND threads=6 r=Chr01 o=Chr01

# a set of small contigs (list file, one name per line, supports # comments and blank lines)
java -Xmx16g -jar GeMoSeq-1.2.3-fast.jar gemoseq g=genome.fa m=merged.bam s=FR_SECOND_STRAND threads=6 rl=small_contigs.list o=small
```

To squeeze memory further: `mrpr=2000000 threads=4`. All other parameters (`mrc`, `mnoir`, ...)
keep their original semantics; see `java -jar GeMoSeq-1.2.3-fast.jar gemoseq`.

## 5. Building from source

Requirements: the official GeMoSeq-1.2.3.jar ([jstacs.de](https://www.jstacs.de/index.php/GeMoSeq)),
htsjdk-2.24.1.jar (Maven Central), JDK 11+, Python 3.

```bash
# Windows uses ';' as classpath separator; Linux ':'
javac -encoding UTF-8 -cp GeMoSeq-1.2.3.jar;htsjdk-2.24.1.jar -d out-fast src/gemoseq/*.java
python3 build/repack_jar.py --base GeMoSeq-1.2.3.jar --htsjdk htsjdk-2.24.1.jar --classes out-fast --out GeMoSeq-1.2.3-fast.jar
```

## 6. Verification

`bench/` contains the synthetic data generator and benchmark scripts:

- `GenTestData.java`: generates a synthetic genome with GT-AG splice sites and a sorted, indexed
  BAM with a controllable depth profile (including a spike hotspot); used for coordinate-level
  diffs between original and fork;
- `bench.sh` / `bench2.ps1`: record wall time, peak RSS, exit code.

See section 3 for the reproduced experiments. Acceptance criteria: 100% coordinate identity at
non-downsampled loci, intron-chain agreement ≥ 99% at downsampled loci.

## 7. Known limitations and semantic differences (stated honestly)

1. **Long-read mode** (`lr=true`) has been adapted but not thoroughly tested; long reads are
   structurally almost unique, so collapsing gains little there.
2. Down-sampling is done per group, which differs slightly from the original per-read coin flips —
   abundance values in downsampled regions may vary by a few percent.
3. With more than two alignments sharing a read name (secondary/supplementary), pair collapsing
   differs slightly from the original idMap semantics; such records are rare under the `q=40`
   MAPQ filter.
4. With `c=false` (collapsing off), output is coordinate-identical to the original (verified) —
   useful as a control.
5. Stale-index detection is based on file mtime; no content-level validation is performed
   (the warning advises re-indexing).
6. In restricted runs (`r`/`rl`), `ReadStats` is computed from the restricted references only,
   so the output matches a split-BAM run exactly; whole-BAM statistics would differ slightly.

## 8. Repository layout

```
GeMoSeq-1.2.3-fast.jar   # ready-to-run fat jar
src/gemoseq/             # the 7 modified/added Java source files (override same-named upstream classes)
build/repack_jar.py      # repackaging script
bench/                   # synthetic data generator and benchmark scripts
LICENSE                  # GPL v3 (inherited from Jstacs)
```

## 9. Acknowledgements and citation

GeMoSeq is the transcript-reconstruction companion of GeMoMa, by Jens Keilwagen et al.; see
[upstream](https://github.com/Jstacs/Jstacs). This repository contains engineering optimizations
only and does not change the core algorithmic design; please cite upstream as well when citing
this work.
