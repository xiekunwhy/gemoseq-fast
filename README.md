# GeMoSeq-fast

[GeMoSeq](https://www.jstacs.de/index.php/GeMoSeq)（[Jstacs/Jstacs](https://github.com/Jstacs/Jstacs) 的 `projects/gemoseq`，曾用名 GeMoRNA）的高深度优化分叉。

面向**多样本合并 BAM 的转录本重建**（de novo 基因组注释的转录证据生产）：原版在高深度合并
BAM 上内存随深度×区域长度失控，本分叉把内存与时间做到**对深度基本不敏感**，并增加了
按染色体/参考序列定向运行、输出前缀命名、CSI 索引支持等工程特性。

- 上游版本：GeMoSeq 1.2.3（2025-11-11）
- 本分叉版本号：`1.2.3-fast`
- 许可证：GPL v3（沿用上游，见 [LICENSE](LICENSE)）

---

## 1. 为什么改

转录本注释流程中，把多个 RNA 样本的比对结果合并成一个大 BAM 再 call 转录本，基因完整度
（BUSCO）明显优于"单样本 call 再合并 GTF"。但合并 BAM 的深度是各样本之和，原版 GeMoSeq 在
这种输入上：

1. 每个基因组区域（region）把所有 reads 以 `SAMRecord` 对象链表存内存，上限是
   `maxcov × 区域长度`——默认 100 × 750 kb = **7500 万条 read/区域**，高深度下直接爆内存；
2. 构图、read 归属、EM 定量**三遍按 read 逐碱基处理**，时间随深度线性增长；
3. `nextSplit()` 每切一个连通组分就分配一次全区长度的节点数组并重过滤边表
   （O(组分数²)）；
4. 平均深度不超过 100× 时降采样不触发——**尖峰深度（rRNA、叶绿体基因）完全绕开降采样**。

## 2. 核心修改

### 2.1 片段结构折叠（ReadGroup）——对深度不敏感的关键

摄取 read 时即转成紧凑表示（外显子块数组 + gap 类型 + 预算 mismatch 数），按
"双端 mate 的完整比对结构 + 链方向"做签名折叠：**N 个结构完全相同的片段只存一份加权重 N**。

- 构图：`nReads += 权重`（与原版逐条累加完全相等）；
- 定量矩阵一行对应一个组，`readWeights` 初始化为组权重——与原版"逐 read 一行、再
  makeUnique 合并求和"在数学上等价；
- mate 配对按 read name 在区域内缓冲合并，保持"双端算一个观测"的原版语义；
- `dummy` 填隙 read 维持原版"全区域共用一个索引"的行为。

### 2.2 降采样与内存上限

- 保留原版 `maxcov`（区域平均深度）降采样语义，新增**绝对上限** `mrpr`（默认 4,000,000
  reads/区域），尖峰深度也会被压住；
- 组级降采样（同权重片段整体取舍）；
- 删除原版每次降采样向 stdout 打印 `#N->M` 的刷屏行为。

### 2.3 算法/数据结构修补

- `nextSplit()` 重写为单次 BFS + 按组分实际跨度分配数组（消灭 O(组分数²) 与重复分配）；
- `Node.addOutgoing` 由"containsKey+get+put"三次哈希改为单次 get；
- **基因组 `.fai` 懒加载**：用到哪条染色体才从 fasta 读哪条（有 `.fai` 时），省掉整个基因组的
  常驻内存；
- **修复多线程竞态**：链特异性模式下正反链 Region 被不同 worker 并发读写导致
  `ConcurrentModificationException`（原版同样潜伏此问题）；对 Region 的关键状态访问加了同步。

### 2.4 htsjdk 升级 + 索引处理

- jar 内 htsjdk 由 2.5.0-SNAPSHOT（2016，不支持 CSI）升级为 **2.24.1**；
- 索引自动探测：`*.bam.bai` / `*.bam.csi` / `*.bai` / `*.csi` 四种命名，`.bai` 优先；
- **索引比 BAM 旧（mtime）时打印醒目警告**——过期索引会静默丢数据（见 4.3）；
- 完全没有索引时回退为流式扫描全 BAM 按参考名过滤（正确但慢，有提示）。

### 2.5 新增参数

| 参数 | 短名 | 说明 | 默认 |
|---|---|---|---|
| Restrict to reference | `r` | 只处理 BAM 中该条参考序列（需索引） | 关 |
| Reference list | `rl` | 文件里每行一个参考序列名，只处理这些并输出到一个 GFF（需索引） | 关 |
| Output prefix | `o` | 输出命名为 `<前缀>.Transcript_Predictions.gff3` 与 `<前缀>.protocol_gemorna.txt`（在 outdir 下，默认当前目录）；中间临时文件也用 `<前缀>.predictions.tmp` | 关 |
| Collapse identical fragments | `c` | 片段折叠开关（对照用） | true |
| Maximum reads per region | `mrpr` | 区域 reads 绝对上限 | 4,000,000 |

原版参数全部兼容不变。

## 3. 性能与正确性

### 3.1 合成数据（2 Mb 基因组、330 基因、含 100 kb 热点区；bench/ 可复现）

| 数据 | 原版 | 修改版 | 输出一致性 |
|---|---|---|---|
| low（30×/3000×，104 万 reads） | 18 s / 5.3 GB | 17 s / **2.4 GB** | **坐标 100% 一致** |
| hi（2000×/100,000×，3760 万 reads） | 157 s / 8.8 GB | **76 s / 4.3 GB** | 坐标 97% 重合* |

\* 差异全部位于双方都触发随机降采样的热点区（原版自身也是固定种子随机降采样，两次不同
实现必然有涨落）。

### 3.2 真实数据（植物基因组，1.1 Gb、15 条参考序列；FR_SECOND_STRAND，threads=6）

| 运行 | 峰值内存 | 耗时 | 说明 |
|---|---|---|---|
| 原版，Chr01 单染色体 BAM（625 MB） | 11.5 GB | 143 s | 6765 转录本 |
| **修改版**，同左 | **8.0 GB** | **101 s** | 6710 转录本 |
| 修改版，全基因组 BAM（4.8 GB）+ `r=Chr01` | 7.5 GB | 205 s | 6713，索引查询 |
| 修改版，全基因组 BAM + `rl=ptg.list o=ptg` | 0.55 GB | 234 s | 小 contig 列表模式 |
| 修改版，全基因组 BAM 全量 + `o=cca` | 11.4 GB | 635 s | 48219 转录本 |

同数据同参数下，修改版与原版 jar 的**内含子链一致率 99.5%**（差异集中在超高深度位点的
末端外显子边界）。

### 3.3 真实数据验证时发现的上游数据坑（重要）

测试中用"索引查询数 reads"与"流式数 reads"交叉验证时发现：测试用全基因组 BAM 自带的
`.csi` 索引是**旧版 BAM 的遗留**（mtime 早于 BAM，query 只能取回 47% 的 reads）。
用 htsjdk 重建索引后两者完全一致。**换 BAM 不换索引会导致静默丢数据**，本分叉现在会检测
索引 mtime 早于 BAM 的情况并打印警告。

## 4. 用法

```bash
# 全基因组
java -Xmx16g -jar GeMoSeq-1.2.3-fast.jar gemoseq g=genome.fa m=merged.bam s=FR_SECOND_STRAND threads=6 o=cca

# 只做 Chr01（索引查询）
java -Xmx16g -jar GeMoSeq-1.2.3-fast.jar gemoseq g=genome.fa m=merged.bam s=FR_SECOND_STRAND threads=6 r=Chr01 o=Chr01

# 只做一批小 contig（列表文件，一行一个名字，支持 # 注释与空行）
java -Xmx16g -jar GeMoSeq-1.2.3-fast.jar gemoseq g=genome.fa m=merged.bam s=FR_SECOND_STRAND threads=6 rl=small_contigs.list o=small
```

内存还想再压：`mrpr=2000000 threads=4`。其他参数（`mrc`、`mnoir` 等）语义同原版，见
`java -jar GeMoSeq-1.2.3-fast.jar gemoseq`。

## 5. 从源码构建

依赖：官方 GeMoSeq-1.2.3.jar（[jstacs.de](https://www.jstacs.de/index.php/GeMoSeq)）、
htsjdk-2.24.1.jar（Maven Central）、JDK 11+、Python 3。

```bash
# Windows 下 classpath 用分号；Linux 用冒号
javac -encoding UTF-8 -cp GeMoSeq-1.2.3.jar;htsjdk-2.24.1.jar -d out-fast src/gemoseq/*.java
python3 build/repack_jar.py --base GeMoSeq-1.2.3.jar --htsjdk htsjdk-2.24.1.jar --classes out-fast --out GeMoSeq-1.2.3-fast.jar
```

## 6. 验证方法

`bench/` 内含合成数据生成器与基准脚本：

- `GenTestData.java`：生成带 GT-AG 剪接位点的合成基因组 + 可控深度剖面（含尖峰热点）的
  排序索引 BAM；用于"原版 vs 修改版"的坐标级一致性 diff；
- `bench.sh` / `bench2.ps1`：记录 wall time、峰值 RSS、退出码。

复现实验见第 3 节。正确性验收标准：非降采样位点坐标 100% 一致，降采样位点内含子链
一致率 ≥ 99%。

## 7. 已知限制与语义差异（如实说明）

1. **长读模式**（`lr=true`）已做适配但无充分测试；长读结构几乎不重复，折叠收益小。
2. 降采样在组层面进行，与原版逐 read 抛硬币的分布略有差异——降采样区域的丰度值会有
   百分之几的合理涨落。
3. 同名比对记录 >2 条（secondary/supplementary）时配对折叠与原版 idMap 语义略有出入；
   `q=40`（MAPQ 过滤）下这类记录很少。
4. `c=false` 关闭折叠时输出与原版逐坐标一致（已验证），可用于对照。
5. 索引过期检测基于文件 mtime，内容级校验不做（交给警告提示重建）。

## 8. 目录结构

```
GeMoSeq-1.2.3-fast.jar   # 可直接运行的 fat jar
src/gemoseq/             # 修改/新增的 7 个 Java 源文件（覆盖上游同名类）
build/repack_jar.py      # 重新打包脚本
bench/                   # 合成数据生成器与基准脚本
LICENSE                  # GPL v3（沿用 Jstacs）
```

## 9. 致谢与引用

GeMoSeq 是 GeMoMa 的配套转录本重建工具，原作者 Jens Keilwagen 等，见
[upstream](https://github.com/Jstacs/Jstacs)。本仓库仅做工程优化，不改变其核心算法设计；
引用时请同时引用上游。
