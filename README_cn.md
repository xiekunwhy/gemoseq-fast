# GeMoSeq-fast

> [English](README.md) | 中文

[GeMoSeq](https://www.jstacs.de/index.php/GeMoSeq)（[Jstacs/Jstacs](https://github.com/Jstacs/Jstacs) 的 `projects/gemoseq`，曾用名 GeMoRNA）的高深度优化分叉，
面向**多样本合并 BAM 的转录本重建**（de novo 基因组注释的转录证据生产）。

原版在高深度合并 BAM 上内存随深度×区域长度失控、线程超 6 无收益；本分叉把内存与时间做到
**对测序深度基本不敏感**，并补齐了大型基因组工程化所需的一整套能力：片段折叠、异步 I/O 流水线、
稀疏 EM、按染色体/参考序列定向运行、CSI 索引支持（含 >512Mb 的关键修复）、输出前缀、
丰度还原与 TPM 计算、可复用的全局 read 统计。

- 上游版本：GeMoSeq 1.2.3（2025-11-11）
- 本分叉版本：`1.2.3-fast`
- 许可证：GPL v3（沿用上游，见 [LICENSE](LICENSE)）

---

## 1. 为什么改

注释流程中，多个 RNA 样本合并成一个大 BAM 再 call 转录本，基因完整度（BUSCO）明显优于
"单样本 call 再合并 GTF"。但合并 BAM 深度是各样本之和，原版 GeMoSeq 在这种输入上：

1. 每个基因组区域把全部 reads 以 `SAMRecord` 对象链表存内存，上限 `maxcov × 区域长度`
   （默认 100 × 750 kb = 7500 万条/区域），尖峰深度（rRNA、叶绿体）完全绕开降采样；
2. 构图、read 归属、EM 定量**三遍按 read 逐碱基处理**，时间随深度线性增长；
3. BAM 摄取（解压+解码+区域构建）单线程，线程超 6 无收益，I/O 成瓶颈；
4. 一些工程缺口：不支持 CSI 索引、不能按染色体定向、输出名不可控、降采样刷屏、
   基因组全量入内存、若干崩溃（含上游 issue #75）。

## 2. 功能与修改总览

### 2.1 片段结构折叠（ReadGroup）——对深度不敏感的核心

摄取 read 时即转成紧凑表示（外显子块数组 + gap 类型 + 预算 mismatch 数），按
"双端 mate 的完整比对结构 + 链方向"做签名折叠：**N 个结构完全相同的片段只存一份加权重 N**。
构图计数 `nReads += 权重`（与逐条累加完全相等）；定量矩阵一行对应一个组，初始权重=组权重，
与原版"逐 read 一行再 makeUnique 合并"在数学上等价。mate 按 read name 在区域内缓冲合并，
保持"双端算一个观测"的原版语义。

### 2.2 降采样与内存上限

- 保留原版 `mrc`（区域平均深度，默认 100）降采样语义，新增**绝对上限** `mrpr`
  （默认 4,000,000 reads/区域），尖峰深度也被压住；
- 组级降采样（同权重片段整体取舍，保持期望）；
- 删除原版每次降采样向 stdout 打印 `#N->M` 的刷屏行为。

### 2.3 异步 I/O 流水线与稀疏计算

- **后台 ReadStats**：插入片段/intron 长度统计与主流程重叠计算；
- **异步解码 + 并行转换**：生产者线程做 BGZF 解压/SAMRecord 解码，4 个转换线程并行做
  cigar/block/mismatch/链向分析，按序号带背压交接（顺序不变、结果确定）；
- **稀疏 EM 定量**：EM 循环按每条 read 的稀疏相容候选列表迭代，替代稠密
  reads×transcripts 矩阵，每次迭代 O(Σ相容数)，结果不变；
- `ReadGraph.remove` 规避 `LinkedList.removeAll` 的 O(n×m)（改 O(n+m)）；
- `Node` 边表懒创建（每碱基内存 ~112B → ~50B）；
- `nextSplit()` 重写为单次 BFS + 按组分实际跨度分配数组（消灭 O(组分数²)）；
- **基因组 `.fai` 懒加载**：用到哪条染色体才从 fasta 读哪条；
- 实测 threads 8-12 是甜点；更大的并行建议按染色体多进程。

### 2.4 并发与崩溃修复

- 修复链特异性模式下正反链 Region 被不同 worker 并发读写导致的
  `ConcurrentModificationException`（原版同样潜伏）；
- **修复覆盖度切分产生空子区域导致的 NPE**（[上游 issue #75](https://github.com/Jstacs/Jstacs/issues/75)：
  完全落在长 intron 内的子区间收不到"完全包含"的 read，worker 拆箱空指针；
  空子区域直接丢弃 + 计算路径加保护）；
- 未比对 read（无参考）安全跳过。

### 2.5 htsjdk 升级、CSI 与索引处理

- jar 内 htsjdk 由 2.5.0-SNAPSHOT（2016，不支持 CSI）升级为 **2.24.1**；
- 索引自动探测：`*.bam.bai` / `*.bam.csi` / `*.bai` / `*.csi` 四种命名，`.bai` 优先；
- **索引比 BAM 旧（mtime）时打印醒目警告**——过期索引会静默丢数据（实测旧索引
  只能取回 47% 的 reads）；
- 无索引时回退为流式扫描全 BAM 按参考名过滤（正确但慢，有提示）；
- **打过补丁的 `CSIIndex.class`**（源码 `src/htsjdk/samtools/CSIIndex.java`）：上游
  htsjdk 2.24 在全参考 CSI 查询时沿 bin 层级走查计算 `minimumOffset`，对 samtools
  （htslib 1.20+）写出的 depth-6 csi（>512Mb 记录被放进 level-1 大 bin），走查会锚定
  到装着染色体尾部的 level-1 bin，`Chunk.optimizeChunkList` 随即将 **~1.07-1.61G
  之前的全部 chunks 静默丢弃**（受限运行只剩尾部）。本补丁在 `startPos <= 0` 时跳过
  该走查，恢复完整覆盖。

### 2.6 定向运行与输出

- `r=<染色体>` / `rl=<列表文件>`：只处理指定参考序列（索引查询）；
- `o=<前缀>`：输出 `<前缀>.Transcript_Predictions.gff3` 与 `<前缀>.protocol_gemorna.txt`
  （在 outdir 下），中间临时文件也用 `<前缀>.predictions.tmp`；
- `ra=true`：把降采样区域的丰度按 1/采样概率还原（TPM 类定量需要）。

### 2.7 统计复用：readstats 工具与 `rs` 参数

jar 内并立 **`readstats`** 工具（与 gemoseq/predictCDS/GAF/Analyzer/merge 同级）：
对全基因组（或限定参考序列）计算 intron 长度分布统计（meanSplit/sdSplit/meanReadLen），
写成文本文件。GeMoSeq 用 `rs=<文件>` 读入后**跳过自身统计步骤**，直接用该文件做剪接图剪枝。
一次全局统计，处处复用——速度与口径统一兼得。

## 3. 参数总表

原版参数全部兼容（`g/m/s/l/sil/lr/mnor/mfor/mnoir/mfoir/p/mrpg/mrpt/pa/sf/mrl/mrc/mfgl/q/mpl/gp/gnwc`），
新增/改变如下：

| 参数 | 短名 | 说明 | 默认 |
|---|---|---|---|
| Restrict to reference | `r` | 只处理 BAM 中该条参考序列（需索引） | 关 |
| Reference list | `rl` | 每行一个参考序列名的文件，只处理这些并输出到一个 GFF（需索引） | 关 |
| Stream full BAM | `sfb` | 限定时不走索引，流式扫全 BAM 按名过滤（慢但免疫索引问题） | false |
| Output prefix | `o` | 输出与中间文件加前缀 | 关 |
| Collapse identical fragments | `c` | 片段折叠开关（对照用） | true |
| Maximum reads per region | `mrpr` | 区域 reads 绝对上限 | 4,000,000 |
| Rescale abundance | `ra` | 丰度按 1/降采样概率还原（TPM 类定量用） | false |
| Read statistics | `rs` | 预计算 read 统计文件（readstats 工具生成），给后跳过自身统计 | 关 |

`readstats` 工具参数：`m=<bam>`（必需）、`o=<输出文件>`（默认 readstats.stats）、
`Shortest intron length`（同 gemoseq 的 sil）、`rl=<列表>`（可选限定范围）。
**输出命名为 `<outdir>/<o>`，protocol 为 `<outdir>/<o>.protocol`**（与 gemoseq 的前缀语义一致）。

## 4. 性能与正确性

### 4.1 合成数据（2 Mb 基因组、330 基因、含 100 kb 热点区；bench/ 可复现）

| 数据 | 原版 | 修改版 | 输出一致性 |
|---|---|---|---|
| low（30×/3000×，104 万 reads） | 18 s / 5.3 GB | 17 s / **2.4 GB** | **坐标 100% 一致** |
| hi（2000×/100,000×，3760 万 reads） | 157 s / 8.8 GB | **76 s / 4.3 GB** | 坐标 97.3% 精确一致* |

\* "精确一致"为严格口径：(染色体, start, end, 链向) 四元组完全相等才算，部分重叠不计。
差异全部位于双方都触发随机降采样的热点区（原版自身也是固定种子随机降采样，
两次不同实现必然有涨落）；内含子链口径 96.9-97.6%。

### 4.2 真实数据（植物基因组，1.1 Gb、15 条参考；FR_SECOND_STRAND，threads=6）

| 运行 | 峰值内存 | 耗时 | 说明 |
|---|---|---|---|
| 原版，Chr01 单染色体 BAM（625 MB） | 11.5 GB | 143 s | 6765 转录本 |
| **修改版**，同左 | **8.0 GB** | **83-101 s** | 6710，内含子链与原版一致率 99.5% |
| 修改版，全基因组 BAM（4.8 GB）+ `r=Chr01` | 7.5 GB | 108 s | 6710，与拆分 BAM **逐坐标一致** |
| 修改版，全基因组 BAM + `rl=ptg.list o=ptg` | 0.55 GB | 234 s | 小 contig 列表模式 |
| 修改版，全基因组 BAM 全量 + `o=cca` | 11.4 GB | 635 s | 48219 转录本 |

线程扩展：6→12 线程约再快 15-20%；**同一输入重复运行输出逐字节一致**（随机种子固定，
且随机决策全部发生在主摄取线程）。

### 4.3 超 512Mb 染色体（depth-6 csi）

本地以 htslib 1.22 规则构造的 600 Mb 染色体 + 手写 csi 验证：打补丁前查询只剩尾部
（4960/9920），打补丁后 `query == stream`（9920=9920），受限运行头尾基因齐全。

## 5. 典型用法

```bash
# 0) 可选但推荐：先算一次全局统计（大 BAM 只扫一遍）
java -jar GeMoSeq-1.2.3-fast.jar readstats m=merged.bam o=genome.stats

# 1) 全基因组一次跑
java -Xmx16g -jar GeMoSeq-1.2.3-fast.jar gemoseq g=genome.fa m=merged.bam s=FR_SECOND_STRAND threads=6 o=cca

# 2) 大染色体逐条并行（推荐的多进程姿势；rs 复用全局统计）
java -Xmx16g -jar GeMoSeq-1.2.3-fast.jar gemoseq g=genome.fa m=merged.bam s=FR_SECOND_STRAND threads=6 r=Chr01 o=Chr01 rs=genome.stats

# 3) 一批小 contig 合并跑（共享统计，避免小样本退化）
java -Xmx16g -jar GeMoSeq-1.2.3-fast.jar gemoseq g=genome.fa m=merged.bam s=FR_SECOND_STRAND threads=6 rl=small_contigs.list o=small rs=genome.stats

# 4) 分染色体结果合并后统一算 TPM
cat chr*.Transcript_Predictions.gff3 | perl scripts/gemoseq_tpm.pl - > all.tpm.gff3
```

内存还想再压：`mrpr=2000000 threads=4`。降采样明显的区域想要丰度无偏：`ra=true`。

## 6. 从源码构建

依赖：官方 GeMoSeq-1.2.3.jar（[jstacs.de](https://www.jstacs.de/index.php/GeMoSeq)）、
htsjdk-2.24.1.jar（Maven Central）、JDK 11+、Python 3。

```bash
# Windows 下 classpath 用分号；Linux 用冒号
javac -encoding UTF-8 -cp GeMoSeq-1.2.3.jar;htsjdk-2.24.1.jar -d out-fast src/gemoseq/*.java src/gemoma/*.java src/htsjdk/samtools/CSIIndex.java
python3 build/repack_jar.py --base GeMoSeq-1.2.3.jar --htsjdk htsjdk-2.24.1.jar --classes out-fast --out GeMoSeq-1.2.3-fast.jar
```

## 7. 验证方法

`bench/` 内含合成数据生成器与基准脚本（`GenTestData.java`、`bench.sh`、`bench2.ps1`、
`ReproNPE.java` 等），可复现第 4 节全部数字。验收标准：非降采样位点坐标 100% 一致，
降采样位点内含子链一致率 ≥ 99%。

## 8. 已知限制与语义差异（如实说明）

1. **长读模式**（`lr=true`）已做适配但无充分测试；长读结构几乎不重复，折叠收益小。
2. 降采样在组层面进行，与原版逐 read 抛硬币的分布略有差异——降采样区域的丰度值会有
   百分之几的合理涨落（可用 `ra=true` 还原期望尺度）。
3. 同名比对记录 >2 条（secondary/supplementary）时配对折叠与原版 idMap 语义略有出入；
   `q=40`（MAPQ 过滤）下这类记录很少。
4. `c=false` 关闭折叠时输出与原版逐坐标一致（已验证），可用于对照。
5. 索引过期检测基于文件 mtime，内容级校验不做（交给警告提示重建）。
6. **ReadStats 统计口径**：默认受限运行（`r`/`rl`）只用受限参考序列做统计（因此与拆分
   BAM 运行逐坐标一致；全基因组口径会略有不同）。如需全局统一判据（推荐用于小 contig
   或要求跨染色体一致时），用 `readstats` 先算全基因组统计、再以 `rs=` 喂给各次运行。
7. 小 contig 单独运行时若可剪接 read 极少，intron 长度方差可能退化为 0 导致过度剪枝——
   用 `rl=` 打包或 `rs=` 全局统计规避。
8. `readstats` 的输出与 protocol 按 `o=` 命名为 `<outdir>/<o>` 与 `<outdir>/<o>.protocol`。

## 9. 目录结构

```
GeMoSeq-1.2.3-fast.jar   # 可直接运行的 fat jar（含 htsjdk 2.24.1 + CSIIndex 补丁）
src/gemoseq/             # 修改/新增的 gemoseq 源文件（覆盖上游同名类）
src/gemoma/ReadStats.java  # ReadStats（含 toFile/fromFile 与限定统计）
src/htsjdk/samtools/CSIIndex.java  # 打过补丁的 htsjdk 类
scripts/gemoseq_tpm.pl   # TPM 计算脚本
build/repack_jar.py      # 重新打包脚本
bench/                   # 合成数据生成器与基准脚本
LICENSE                  # GPL v3（沿用 Jstacs）
```

## 10. 致谢与引用

GeMoSeq 是 GeMoMa 的配套转录本重建工具，原作者 Jens Keilwagen 等，见
[upstream](https://github.com/Jstacs/Jstacs)。本仓库仅做工程优化，不改变其核心算法设计；
引用时请同时引用上游。
