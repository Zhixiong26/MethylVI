# MethylVI 5-kb分析流程：IR与NR项目

本目录保存当前10个样本（5个IR、5个NR）的ALLCools与MethylVI可复现流程。分析从MethSCAn 300k QC细胞出发，重建5-kb mCG整数计数，以`sample_id`作为batch covariate训练MethylVI，并输出latent、UMAP、Leiden及标签引导UMAP。

文档更新至2026-08-12。服务器大任务统一通过`dsub`提交，不在登录节点直接运行分析。

## 1. 快速入口与当前完成状态

服务器工作目录：

```bash
cd /share/home/rzli/MethylVI/20260810/scripts
```

统一入口：

```bash
bash 09_run_pipeline.sh prepare
bash 09_run_pipeline.sh verify
bash 09_run_pipeline.sh build
bash 09_run_pipeline.sh train
bash 09_run_pipeline.sh plots
bash 09_run_pipeline.sh supervised
bash 09_run_pipeline.sh test
bash 09_run_pipeline.sh all
```

正式流程现在默认启用GRCh38 blacklist。`prepare`从原始ALLC开始完整复现；已经有历史MCDS时，也可用快捷入口避免重复生成MCDS：

```bash
bash 09_run_pipeline.sh blacklist
```

`blacklist`快捷入口复用历史`mcg_5kb.mcds`，只生成blacklist过滤后的新H5AD和聚类结果；它不属于`all`。只有读取历史未过滤结果时才设置`MVI_USE_BLACKLIST=0`。

`all`只执行`verify → build → train → plots → supervised`，不会自动执行耗时的`prepare`，也不会重复运行测试。`prepare`直接使用`MVI_ALLCOOLS_ENV`；其他阶段使用`MVI_CONDA_ENV`。

| 项目 | 当前状态 |
|---|---|
| 10个样本原始ALLC及TBI | 已验收，58,534个细胞 |
| MethSCAn 300k细胞QC | 已验收，6,199个细胞 |
| ALLCools 5-kb MCDS与聚类 | 正式默认版本已完成blacklist与低频过滤；6,199×230,306 |
| 当前SCANPY注释 | 使用`20260810/Result0810`版本，匹配5,765/6,199 |
| MethylVI H5MU | 历史231,648-bin版本已完成；正式230,306-bin版本待重建 |
| 正式CPU训练 | 历史版本已完成；正式blacklist版本待H5MU重建后重新训练 |
| 校正前后普通图 | 已用当前注释重画 |
| supervised UMAP | 0.2、0.5、0.7、0.9均已完成，共12张PDF |
| 重构后服务器测试 | 任务`164116`通过：8项单元测试和两轮CPU smoke test均成功 |
| 重构后真实输入审计 | 任务`164118`通过：6,199个细胞、10个样本和全部ALLC完全匹配 |
| 100k版本真实输入审计 | 任务`164132`成功：6,199个细胞、100,206个bins、10样本、6,199份ALLC全部匹配，退出码0 |
| 100k版本H5MU构建 | 任务`164133`成功：6,199行全部组装，H5MU约0.49 GiB，退出码0 |
| 100k版本CPU训练 | 任务`164134`成功：early stopping于第80/500个epoch，80条记录，退出码0 |
| 50k版本阈值估计 | 2026-08-13完成：`MVI_HYPO_PERCENT=2.669785`，内部阈值为非零细胞数>165，预计49,947个bins |
| 50k版本blacklist聚类 | 任务`164163`成功：6,199×49,947 H5AD，10个显著LSI成分，退出码0 |
| 50k版本H5MU构建 | 任务`164165`成功：6,199行全部组装，H5MU约0.26 GiB，退出码0 |
| 50k版本CPU训练 | 任务`164166`成功：96 CPU、120 GiB；第69/500个epoch early stopping，69条记录，退出码0 |
| 230k版本H5MU重建 | 任务`164171`成功：96 CPU、120 GiB；6,199行全部组装，H5MU约1.03 GiB，退出码0 |
| 230k版本CPU训练 | 任务`164172`成功：120 CPU、120 GiB；第78/500个epoch early stopping，78条记录，退出码0 |
| 230k版本普通UMAP | 任务`164173`成功：生成校正后 sample、condition、cell type 三张PDF |
| 230k版本supervised UMAP | 任务`164174`成功：生成0.2、0.5、0.7、0.9四组权重，共12张PDF |
| blacklist运行前验收 | 2026-08-12通过：ENCFF356LFX MD5/gzip、独立路径、ALLCools、pybedtools和bedtools均正常 |
| blacklist独立聚类 | 任务`164127`成功：新H5AD为6,199×230,306，8个显著LSI成分，退出码0 |
| blacklist约10万bin重筛选 | 任务`164130`成功：6,199×100,206 H5AD，9个显著LSI成分，退出码0 |

当前设计固定为`IR01–IR05`和`NR01–NR05`。300k表示每个细胞至少覆盖300,000个CpG位点，并非保留300,000个特征。正式默认版本使用230,306个经过ENCODE GRCh38 blacklist、`hypo-score > 0.95`二值化及低频过滤的5-kb bins；H5AD中的聚类分数不作为MethylVI计数。原231,648-bin结果仅作为历史未过滤版本保留。

## 2. 正式执行顺序

| 步骤 | 入口stage | 执行脚本 | 主要输入 | 主要输出 |
|---:|---|---|---|---|
| 1 | `prepare` | `02_prepare_allcools.sh`、`03_cluster_allcools.py` | 原始ALLC、300k白名单、hg38长度、ENCFF356LFX | 6,199个ALLC/TBI链接、MCDS、blacklist统计及聚类H5AD |
| 2 | `verify` | `04_verify_inputs.py` | H5AD、ALLC、样本表、当前注释 | 独立`input_audit.json` |
| 3 | `build` | `05_build_methylvi_input.py` | 6,199个ALLC、230,306个bins | H5MU、逐细胞检查点、build摘要 |
| 4 | `train` | `06_train_methylvi.py` | H5MU、`sample_id` batch | 模型、20维latent、UMAP、Leiden、训练记录 |
| 5 | `plots` | `07_plot_embeddings.py --stage all` | ALLCools及MethylVI坐标、当前注释 | 校正前3张、校正后3张PDF |
| 6 | `supervised` | `08_plot_supervised_umap.py` | 固定MethylVI latent、当前cell type | 4组坐标、12张PDF、H5AD和JSON摘要 |

测试不是正式分析步骤。`test`单独运行`tests/test_mvi_utils.py`和`tests/test_methylvi_smoke.py`；后者使用合成数据进行两轮CPU训练。

大任务提交示例：

```bash
dsub \
  -n methylvi_cpu_formal \
  -R "cpu=32;mem=65536MB" \
  --cwd /share/home/rzli/MethylVI/20260810/scripts \
  -oo scheduler_logs/methylvi_cpu_formal.%J.out \
  -eo scheduler_logs/methylvi_cpu_formal.%J.err \
  env MVI_THREADS=32 MVI_ACCELERATOR=cpu \
  bash 09_run_pipeline.sh train
```

服务器更新时应一次性替换整个`scripts/`目录并移除旧编号脚本，避免新旧入口并存。历史数据、模型和图片继续保留；正式blacklist版本需要重新执行`verify → build → train → plots → supervised`。

## 3. 文件作用

| 文件 | 作用 |
|---|---|
| `00_config.sh` | 输入输出路径、样本结构、QC、batch correction、环境、资源及模型参数 |
| `01_sample_metadata.tsv` | 10个实际样本的`sample_id`与`condition`对应表 |
| `02_prepare_allcools.sh` | 按300k白名单整理原始ALLC，生成MCDS并调用03；`.cov.gz`仅作备用 |
| `03_cluster_allcools.py` | 默认执行GRCh38 blacklist、5-kb低频过滤、筛选摘要、LSI、聚类、t-SNE和UMAP |
| `04_verify_inputs.py` | 审计H5AD、ALLC、区域、样本结构和当前SCANPY注释 |
| `05_build_methylvi_input.py` | 从ALLC聚合整数`mc/cov`层，支持逐细胞检查点复用并写入H5MU |
| `06_train_methylvi.py` | 训练MethylVI并生成latent、邻居图、UMAP、Leiden及运行摘要 |
| `07_plot_embeddings.py` | `--stage before\|after\|all`；动态读取当前注释并生成6张普通图 |
| `08_plot_supervised_umap.py` | 在固定latent上生成四种`target_weight`的标签引导UMAP |
| `09_run_pipeline.sh` | `prepare/blacklist/verify/build/train/plots/supervised/test/all`统一入口 |
| `mvi_utils.py` | ID标准化、元数据/注释合并、区域解析、ALLC聚合、JSON和绘图公共函数 |
| `requirements.lock.txt` | 正式环境的软件版本记录 |
| `hg38.canonical.chrom.sizes` | hg38 canonical 24条染色体长度 |
| `ENCFF356LFX_GRCh38_blacklist.bed.gz` | ENCODE4 GRCh38 exclusion list；固定MD5并随部署目录同步 |
| `tests/test_mvi_utils.py` | 不读取真实数据的公共函数单元测试 |
| `tests/test_methylvi_smoke.py` | 32×64合成计数的两轮CPU MethylVI smoke test |
| `MethylVI_vs_yuanpei_workflow_comparison.md` | 当前流程与`yuanpei`流程的方法与工程比较 |

## 4. 参数配置

### 4.1 路径和样本

| 内容 | 当前值 |
|---|---|
| 脚本目录 | `/share/home/rzli/MethylVI/20260810/scripts` |
| 图像目录 | `/share/home/rzli/MethylVI/20260810/result/blacklist_f0p2` |
| 数据根目录 | `/share/LCZX_Data/data/allcools` |
| 原始样本目录 | `/share/LCZX_Data/data/allcools/25110891_<sample>_Met` |
| ALLCools 300k输出 | `/share/LCZX_Data/data/allcools/methylvi_5kb_300k_blacklist_f0p2` |
| 聚类H5AD | `/share/LCZX_Data/data/allcools/methylvi_5kb_300k_blacklist_f0p2/mcg_5kb.clustered.h5ad` |
| MethylVI输出根目录 | `/share/LCZX_Data/data/allcools/methylVI_results_300k_blacklist_f0p2` |
| MethylVI输入H5MU | `/share/LCZX_Data/data/allcools/methylVI_results_300k_blacklist_f0p2/methylvi_5kbin_input.h5mu` |
| 正式模型结果 | `/share/LCZX_Data/data/allcools/methylVI_results_300k_blacklist_f0p2/results_ir_nr` |
| 当前SCANPY注释 | `/share/home/rzli/SCANPY/20260810/Result0810/annotation/02_cell_annotation_all_cells.csv` |
| MethylVI环境 | `methylvi` |
| ALLCools环境 | `/share/home/rzli/miniconda3/envs/allcools` |

历史未做blacklist的结果保留在以下路径，只在显式设置`MVI_USE_BLACKLIST=0`时读取：

| 内容 | 历史版本路径 |
|---|---|
| ALLCools输出 | `/share/LCZX_Data/data/allcools/methylvi_5kb_300k` |
| 聚类H5AD | `/share/LCZX_Data/data/allcools/methylvi_5kb_300k/mcg_5kb.clustered.h5ad` |
| MethylVI输出 | `/share/LCZX_Data/data/allcools/methylVI_results_300k` |
| 图片 | `/share/home/rzli/MethylVI/20260810/result/01_before_methylvi`及`02_after_methylvi` |

| 样本 | 分组 | 原始ALLC | 300k QC细胞 |
|---|---|---:|---:|
| IR01 | IR | 7,981 | 583 |
| IR02 | IR | 6,070 | 965 |
| IR03 | IR | 7,383 | 914 |
| IR04 | IR | 8,171 | 318 |
| IR05 | IR | 5,392 | 439 |
| NR01 | NR | 4,340 | 306 |
| NR02 | NR | 5,672 | 755 |
| NR03 | NR | 4,285 | 658 |
| NR04 | NR | 7,057 | 722 |
| NR05 | NR | 2,183 | 539 |
| 合计 | 5 IR + 5 NR | 58,534 | 6,199 |

QC后IR为3,219个细胞，NR为2,980个细胞。筛选参数为`min_sites=300000`、`max_sites=10000000`、`min_meth=55`、`max_meth=none`。

### 4.2 ALLCools与H5MU

| 环节 | 参数 | 当前值 |
|---|---|---|
| ALLC来源 | 优先级 | 原始`*_allc.gz`及TBI优先，cov仅备用 |
| MCDS | region/context | `chrom5k 5000`、`CGN` |
| MCDS | quantifier | `count`和`hypo-score cutoff=0.9` |
| score matrix | quantifier | `CGN hypo-score` |
| blacklist过滤 | 正式默认版本 | `ENCFF356LFX`、`f=0.2`；实际删除14,312个唯一bins |
| binarize | cutoff | `0.95` |
| 低频bin过滤 | `filter_regions()` | 显式`hypo_percent=0.5` |
| LSI | algorithm/seed | `arpack` / `0` |
| significant PC | p cutoff | `0.1`（本次得到8个成分） |
| 邻居/初始Leiden | 参数 | `25` / resolution `1.0` |
| t-SNE | 参数 | Euclidean、perplexity 30、exaggeration -1 |
| ConsensusClustering | 参数 | neighbors 25、resolution 0.5、repeats 500、min cluster 10、rate 0.5、seed 0 |
| H5MU | 正式目标形状 | 6,199×230,306（待重建） |
| H5MU | 层/dtype | `mc`和`cov`，预计继续使用`uint16`，以build验收为准 |
| H5MU build | workers | `32`，逐细胞检查点可复用；特征集合变化后旧检查点不得直接误复用 |

### 4.3 5-kb bins筛选与blacklist依据

Blacklist（exclusion list）是参考基因组中容易因重复序列、多重比对、组装困难或异常高信号产生技术伪影的区域集合。Blacklist过滤按基因组坐标排除不可靠区域；`filter_regions()`则按当前数据中显著hypomethylation信号的出现频率删除低频bins，两者解决的问题不同，不能相互替代。

ALLCools官方mCG 5-kb教程采用以下顺序：

```text
细胞QC
→ 删除与参考基因组blacklist重叠的5-kb bins
→ 提取CGN hypo-score矩阵
→ 以0.95阈值二值化
→ 删除只在极少数细胞中非零的bins
→ TF-IDF、LSI、邻居图和聚类
```

官方API中`remove_black_list_region(..., f=0.2)`的默认`f=0.2`表示当一个5-kb bin至少20%（约1,000 bp）与blacklist区域重叠时将其排除；`filter_regions(..., hypo_percent=0.5)`要求区域在至少约0.5%的细胞中有二值化后的非零hypo信号。对6,199个细胞，0.5%约为31个细胞。具体保留判断和取整以服务器ALLCools 1.1.1源码为准。

历史结果与正式默认流程必须分开记录：

| 项目 | 历史未过滤版本 | 正式默认版本 |
|---|---|---|
| 参考基因组 | hg38 canonical 24条染色体 | 不变 |
| blacklist | 未使用 | ENCODE GRCh38 exclusion list `ENCFF356LFX` |
| blacklist overlap | 不适用 | `f=0.2` |
| hypo-score二值化 | `cutoff=0.95` | 不变 |
| 低频bin过滤 | 软件默认`hypo_percent=0.5` | 显式`hypo_percent=0.5`并写入摘要 |
| 最终bins | 231,648 | 230,306（任务`164127`已验收） |
| MethylVI | 已按231,648个bins训练 | 需按230,306个bins重建H5MU并重新训练 |

正式流程记录blacklist文件来源、assembly、accession、MD5、重叠阈值、blacklist删除数量、低频删除数量及最终bin清单。完整重现时从`prepare`开始生成新MCDS；本次迁移利用`blacklist`快捷入口复用已有MCDS。最终H5AD特征集合改变后，必须重新执行`verify → build → train → plots → supervised`。历史231,648-bin版本保留但不再作为默认输入。

上述方案已经整合为默认profile。参考文件为`ENCFF356LFX_GRCh38_blacklist.bed.gz`，期望MD5为`393688b4f06c9ce26165d47433dd8c37`；`03_cluster_allcools.py`会在运行前强制复核。运行后自动生成：

- `feature_filter_summary.json`：初始bin数、blacklist删除数、低频删除数和最终bin数；
- `retained_5kb_bins.tsv.gz`：最终保留bin的坐标清单；
- `mcg_5kb.clustered.h5ad`：新版本聚类结果。

2026-08-12服务器预检已经确认：文件MD5正确、gzip完整，ALLCools 1.1.1、pybedtools 0.10.0及`/share/home/rzli/miniconda3/envs/allcools/bin/bedtools`均可用；新旧profile路径隔离正确。

任务`164127`的实际筛选及聚类统计如下：

| 指标 | 数量 |
|---|---:|
| 初始5-kb bins | 617,665 |
| blacklist去重后实际删除 | 14,312 |
| blacklist后剩余 | 603,353 |
| 随后由`hypo_percent=0.5`删除 | 373,047 |
| blacklist版本最终保留（f0.2，hypo=0.5） | 230,306 |
| 约10万bin版本最终保留（f0.2，hypo=1.169543） | 100,206 |
| 约10万bin版本显著LSI成分 | 9 |
| 230,306-bin版本显著LSI成分 | 8 |
| 50,000-bin版本最终保留（f0.2，hypo=2.669785） | 49,947 |
| 50,000-bin版本显著LSI成分 | 10 |

约10万bin正式筛选版本（任务`164130`）使用`MVI_HYPO_PERCENT=1.169543`，ALLCools内部阈值为非零细胞数`>72`，最终保留`100,206`个bins；blacklist后删除14,312个唯一bins，低频过滤再删除503,147个bins，显著LSI成分为9。该版本是新的正式候选输入，必须单独重建H5MU和MethylVI模型。

50k版本（任务`164163`）使用`MVI_HYPO_PERCENT=2.669785`，由100,206-bin H5AD的非零细胞频数反推，ALLCools内部阈值为非零细胞数`>165`，实际保留`49,947`个bins；blacklist后删除14,312个唯一bins，低频过滤再删除553,406个bins，显著LSI成分为10。该版本使用独立的`blacklist_f0p2_50k`输出profile，不覆盖100,206-bin版本。

为将最终bin数调整到约100,000，已在6,199个细胞的blacklist版本H5AD上按ALLCools 1.1.1源码反推：`MVI_HYPO_PERCENT=1.169543`使内部`n_cell=int(6199*percent/100)=72`，筛选条件为非零细胞数`>72`，预计保留`100,206`个bins。该值需要重新运行blacklist聚类生成H5AD；随后必须重建H5MU并重新训练MethylVI，不能继续使用230,306-bin版本的H5MU或模型。

ALLCools控制台打印了`14318 chrom5k features removed`，而MCDS维度实际从617,665降为603,353，即实际删除14,312个唯一bins。根据ALLCools 1.1.1的实现，这是因为bedtools相交结果中同一5-kb bin可与多个blacklist区间相交，打印值按相交记录计数；矩阵筛选通过`isin`按唯一bin删除。项目统计和后续形状必须采用`initial_5kb_bins - bins_after_blacklist = 14,312`，不能采用14,318。

服务器正式运行前先做轻量检查：

```bash
cd /share/home/rzli/MethylVI/20260810/scripts

md5sum ENCFF356LFX_GRCh38_blacklist.bed.gz
/share/home/rzli/miniconda3/envs/allcools/bin/python - <<'PY'
import shutil
import pybedtools
print("pybedtools:", pybedtools.__version__)
print("bedtools:", shutil.which("bedtools"))
assert shutil.which("bedtools"), "ALLCools blacklist过滤需要bedtools可执行文件"
PY

bash -c '
source 00_config.sh
printf "source MCDS=%s\nnew H5AD=%s\nnew MethylVI root=%s\n" \
  "$MVI_SOURCE_MCDS" "$MVI_H5AD" "$MVI_ROOT"
'
```

确认检查通过后由计算节点运行；首次只运行到新H5AD和统计文件，不立即重建MethylVI：

```bash
mkdir -p scheduler_logs

dsub \
  -n allcools_blacklist_f0p2 \
  -R "cpu=32;mem=65536MB" \
  --cwd /share/home/rzli/MethylVI/20260810/scripts \
  -oo scheduler_logs/allcools_blacklist_f0p2.%J.out \
  -eo scheduler_logs/allcools_blacklist_f0p2.%J.err \
  env MVI_THREADS=32 \
  bash 09_run_pipeline.sh blacklist
```

任务成功后先读取`feature_filter_summary.json`，确认实际删除量，再决定是否依次提交：

```bash
bash 09_run_pipeline.sh verify
bash 09_run_pipeline.sh build
MVI_ACCELERATOR=cpu bash 09_run_pipeline.sh train
bash 09_run_pipeline.sh plots
bash 09_run_pipeline.sh supervised
```

以上大任务仍应分别通过`dsub`提交，不能直接在登录节点运行。默认已经启用blacklist；只有检查历史结果时才临时设置`MVI_USE_BLACKLIST=0`。

证据边界如下：

- ENCODE blacklist由Amemiya、Kundaje和Boyle系统定义，是功能基因组测序中成熟的技术异常区域QC，不是2026年新提出的算法。
- ALLCools官方5-kb教程在提取hypo-score矩阵前明确执行blacklist过滤，因此对本项目具有最直接的软件流程依据。
- 2023年Nature成年小鼠脑单细胞DNA甲基化图谱的方法明确包含coverage和ENCODE blacklist基础特征过滤，说明该策略已用于大规模单细胞甲基化研究。
- 近期Amethyst单细胞DNA甲基化研究在hg38中使用`black_list_fraction=0.2`，并说明采用ALLCools开发者建议，构成人类数据的直接实践证据。
- 2025年的“Beyond Blacklists”研究指出，不同生成策略、比对器和参考版本会影响exclusion sets。因此blacklist仍可作为常规QC，但必须固定文件版本和assembly，不能把它解释为绝对的生物学无效区域。

### 4.4 MethylVI正式参数

| 参数 | 当前值 |
|---|---|
| batch key | `sample_id` |
| modality/context | `mCG` |
| count layers | `mc` / `cov` |
| likelihood / dispersion | `betabinomial` / `region` |
| latent / hidden / layers | `20` / `128` / `1` |
| batch size | `32` |
| 最大epochs | `500`，early stopping开启 |
| accelerator / devices | 正式任务`cpu` / `1` |
| seed | `0` |
| latent neighbors | `15`，representation=`X_methylVI` |
| Leiden | resolution `1.0`、seed `0` |
| CPU线程 | 正式任务`32` |

学习率、优化器、训练/验证划分和early-stopping patience未显式覆盖，使用锁定的scvi-tools 1.3.3默认值。`sample_id`与IR/NR分组绑定，因此样本校正可能同时削弱真实condition差异；解释时必须保留校正前图。如果以后有独立技术批次，应改用`technical_batch`。

### 4.5 普通图与supervised UMAP

普通图散点`size=4`、`alpha=0.75`、无边线、随机绘制seed 0；类别不超过20时使用`tab20`，否则使用`gist_ncar`。`07_plot_embeddings.py`只重新着色既有坐标，不重训模型。

| supervised参数 | 当前值 |
|---|---|
| target key | `cell_type` |
| target weights | `0.2 0.5 0.7 0.9` |
| neighbors | `15` |
| metric / target metric | `euclidean` / `categorical` |
| min_dist / spread | `0.5` / `1.0` |
| seed | `0` |

未注释或`exclude_from_main_analysis=True`的细胞使用`-1`作为未标记目标：仍显示在图中，但标签不会拉近它们。weight越大，cell type标签对二维布局影响越强；这些图是标签引导可视化，不能代替无监督latent证据。

## 5. 已验收结果和服务器任务记录

### 5.1 关键结果

- 正式默认ALLCools H5AD为6,199×230,306：从617,665个初始bins删除14,312个唯一blacklist bins，再执行二值化及低频过滤；显著LSI成分为8。
- 历史MCDS约3.8 GiB，历史未做blacklist的H5AD为6,199×231,648，仅作为旧结果保留。
- H5AD cell ID唯一，包含`X_pca`、`X_tsne`和`X_umap`；L1四簇为c0 3,600、c1 1,398、c2 665、c3 536，十折准确率0.964。
- 当前SCANPY表56,746行，匹配5,765个，未匹配434个，排除17个；cell type可更新，样本/分组以`01_sample_metadata.tsv`为准。
- 历史H5MU为6,199×231,648，最大`mc=7172`、最大`cov=8287`，满足`0≤mc≤cov`；正式230,306-bin H5MU尚待重建验收。
- 正式训练89条记录（epoch 0–88），由early stopping停止；`train_loss_epoch`从92,240.95降至76,516.16。
- 正式结果包含461 MiB模型、485 KiB latent、约3.0 MiB embedding H5AD、UMAP表、训练历史和运行摘要。
- 环境已验收scvi-tools 1.3.3、Torch 2.12.1、AnnData 0.11.4、MuData 0.3.9、Scanpy 1.11.5、igraph 1.0.0和leidenalg 0.12.0。完整记录见`requirements.lock.txt`。

### 5.2 dsub任务历史

| 任务 | 阶段与资源 | 结果、时间、内存或修复 |
|---:|---|---|
| `163899` | 全部58,534细胞stage-only；4 CPU、16 GiB | 成功；约400秒、峰值321 MiB；验证10样本、5 IR、5 NR |
| `163900` | 300k白名单stage-only；4 CPU、16 GiB | 成功；134秒、峰值135 MiB；确认6,199个ALLC和TBI |
| `163901` | 完整ALLCools；32 CPU、194,560 MiB | 成功；约5小时、峰值约20.3 GiB；生成3.8 GiB MCDS与聚类H5AD |
| `164085` | 首次H5MU build；32 CPU、194,560 MiB | 逐细胞计数完成后失败：AnnData不能直接写`numpy.memmap`；已改为`np.asarray()`视图 |
| `164090` | H5MU build重试；32 CPU、65,536 MiB | 成功；复用全部6,199个检查点，生成约1.03 GiB H5MU |
| `164094` | CPU一轮基准；32 CPU、65,536 MiB | 训练1轮56.44秒、峰值9,579 MiB；随后因缺少igraph失败；安装igraph/leidenalg并验收Leiden |
| `164095` | 正式CPU训练；32 CPU、65,536 MiB | 成功；4,819秒、峰值9,358 MiB、平均CPU约18.38核；89轮early stopping |
| `164096` | 校正后普通图；4 CPU、16,384 MiB | 成功；9秒、峰值82 MiB；生成3张有效单页PDF |
| `164097` | 首次supervised UMAP；2 CPU、16 GiB | 失败；17秒、峰值约500 MiB；Categorical列`fillna(False)`报错，未产生结果 |
| `164099` | 当前注释审计；2 CPU、8 GiB | 成功；10秒、峰值437 MiB；匹配5,765/6,199 |
| `164100` | 当前注释校正前重画；2 CPU、8 GiB | 成功；18秒、峰值884 MiB |
| `164101` | 当前注释校正后重画；2 CPU、8 GiB | 成功；4秒、峰值95 MiB |
| `164102` | supervised UMAP重试；2 CPU、16 GiB | 成功；69秒、峰值709 MiB；生成4组坐标与12张PDF |
| `164116` | 重构后完整测试；2 CPU、16 GiB | 成功；20秒、峰值544 MiB；8项单元测试及32×64合成数据两轮CPU训练全部通过 |
| `164118` | 重构后真实输入审计；2 CPU、8 GiB | 成功；13秒、峰值460 MiB；6,199个细胞、231,648个bins、10样本及6,199份ALLC全部通过 |
| `164127` | GRCh38 blacklist独立聚类；32 CPU、64 GiB | 成功；248秒、峰值15,451 MiB；删除14,312个唯一blacklist bins，生成6,199×230,306 H5AD及8个显著LSI成分 |
| `164129` | GRCh38 blacklist约10万bin重筛选；32 CPU、64 GiB | 失败；18秒、峰值621 MiB；pybedtools找不到`intersectBed`，旧230,306-bin文件未被视为本次结果 |
| `164130` | GRCh38 blacklist约10万bin重筛选重试；32 CPU、64 GiB | 成功；退出码0；生成6,199×100,206 H5AD、9个显著LSI成分；显式传入ALLCools PATH和`MVI_HYPO_PERCENT=1.169543` |
| `164132` | 100k版本真实输入审计；2 CPU、8 GiB | 成功；7秒、峰值418 MiB；6,199个细胞、100,206个bins、10样本及6,199份ALLC全部通过；预计稠密`mc/cov`约2.31 GiB |
| `164133` | 100k版本H5MU构建；32 CPU、64 GiB | 成功；6,199/6,199行组装完成，H5MU约0.49 GiB；输出根目录`methylVI_results_300k_blacklist_f0p2_100k` |
| `164134` | 100k版本正式CPU训练；64 CPU、122,880 MB | 成功；1,861秒、峰值4,532 MiB；第80/500个epoch early stopping，80条记录，最佳`elbo_validation=46075.895`，20维latent及下游结果已生成 |
| `164163` | 50k版本GRCh38 blacklist重筛选；32 CPU、64 GiB | 成功；退出码0；生成6,199×49,947 H5AD、10个显著LSI成分；`MVI_HYPO_PERCENT=2.669785` |
| `164165` | 50k版本H5MU构建；32 CPU、64 GiB | 成功；3,160秒、峰值2,836 MiB；6,199/6,199行组装完成；H5MU约0.26 GiB |
| `164166` | 50k版本正式CPU训练；96 CPU、122,880 MB | 成功；987秒、峰值2,560 MiB；第69/500个epoch early stopping，最佳`elbo_validation=27898.900`；20维latent及下游结果已生成 |
| `164171` | 230k版本H5MU重建；96 CPU、122,880 MB | 成功；1,529秒、峰值11,602 MiB；6,199/6,199行组装完成；H5MU约1.03 GiB |
| `164172` | 230k版本正式CPU训练；120 CPU、122,880 MB | 成功；3,509秒、峰值9,432 MiB；第78/500个epoch early stopping，78条记录；退出码0 |
| `164173` | 230k版本普通UMAP；4 CPU、16,384 MB | 成功；生成`02_after_methylvi`下3张PDF；退出码0 |
| `164174` | 230k版本supervised UMAP；4 CPU、16,384 MB | 成功；生成4个target weight目录、12张PDF；退出码0 |

`164097`的修复是先将`exclude_from_main_analysis`转为Pandas字符串，再填充缺失并解析布尔值。`164094`的模型训练本身成功，错误发生在训练后的Leiden阶段。

`164127`日志中的两条Scanpy `No data for colormapping provided via 'c'`为绘图参数提示；任务正常写出结果并以退出码0结束，不影响H5AD、LSI或聚类结果。

任务`164129`失败原因是dsub计算节点的作业环境没有把`/share/home/rzli/miniconda3/envs/allcools/bin`加入`PATH`，导致pybedtools找不到`intersectBed`；登录节点的预检并不能代表作业节点环境。入口脚本现已显式加入`MVI_ALLCOOLS_ENV/bin`，重新上传后再提交。任务留下的旧`feature_filter_summary.json`和230,306-bin H5AD来自前一次成功运行，不能作为`164129`的结果。

监控`164133`时出现的`djob` token expired只影响登录节点的任务查询，不影响作业执行；日志中的`wrote ...h5mu`和`EXIT_CODE: 0`表明构建已成功完成。

## 6. 输出目录

```text
/share/home/rzli/MethylVI/20260810/result/blacklist_f0p2/
├── 01_before_methylvi/
│   ├── allcools_original_embedding_cell_type.pdf
│   ├── allcools_original_embedding_sample_id.pdf
│   └── allcools_original_embedding_condition.pdf
├── 02_after_methylvi/
│   ├── methylvi_umap_cell_type.pdf
│   ├── methylvi_umap_sample_id.pdf
│   └── methylvi_umap_condition.pdf
└── 03_supervised_umap/
    ├── target_weight_0p2/
    ├── target_weight_0p5/
    ├── target_weight_0p7/
    └── target_weight_0p9/

/share/LCZX_Data/data/allcools/methylVI_results_300k_blacklist_f0p2/
├── methylvi_5kbin_input.h5mu
├── count_rows/
├── build_summary.json
└── results_ir_nr/
    ├── sample_by_condition.csv
    ├── training_history.csv
    ├── model/
    ├── latent_representation.npy
    ├── methylvi_embedding.h5ad
    ├── cell_annotations_umap.tsv.gz
    ├── supervised_umap/
    │   ├── methylvi_supervised_umap.h5ad
    │   ├── supervised_umap_summary.json
    │   └── target_weight_*_coordinates.tsv.gz
    └── run_summary.json
```

正式图片写入项目`result/blacklist_f0p2/`；数据、模型和坐标写入LCZX数据盘的独立blacklist目录。更新cell type只需重新运行`verify`、`plots`和`supervised`，不需要重建H5MU或重训MethylVI。

## 7. 失败修复与已知问题

- 登录Shell默认未初始化Conda；入口会加载`MVI_CONDA_INIT`。只在已经激活正确环境时使用`MVI_SKIP_CONDA=1`。
- 服务器无GPU，正式提交必须显式设置`MVI_ACCELERATOR=cpu`。`pip check`中的`nvidia-cusparselt-cu13`平台提示不影响已验收CPU运行。
- 多进程阶段把OpenBLAS、OMP、MKL、NumExpr及Numba内部线程设为1，避免进程×线程过度并行；训练脚本单独设置PyTorch线程。
- MuData、Lightning、DataLoader worker和Scanpy后端提示属于FutureWarning或性能提示，不等于失败。
- Scanpy未来可能更改Leiden默认后端；本次结果按当前`leidenalg`实现验收，不在分析中途改变后端。
- H5MU逐细胞检查点位于`count_rows/*.npz`；manifest防止不同输入误复用。不要同时启动两个写入相同`MVI_ROOT`或`MVI_RESULTS`的任务。
- `model.save(..., overwrite=True)`会覆盖同名模型目录，重新训练前必须确认结果目录或另设`MVI_RESULTS`。
- 交互Shell保护条件不要直接使用`exit`，否则会结束SSH会话。文档中的`<JOB_ID>`是占位符，不能原样输入Shell。
- 当前暂不做batch correction定量评估。普通UMAP只能用于定性观察，不能证明condition信号保留，也不能排除过校正。
- 正式默认版本为230,306个blacklist-filtered bins。历史231,648-bin H5MU/模型仍可通过`MVI_USE_BLACKLIST=0`读取，但不得与正式版本混用；正式版本必须使用独立H5MU和重新训练的模型。

本地MethSCAn上游资料路径：

```text
/Users/luozhixiong/Library/Mobile Documents/com~apple~CloudDocs/Documents/PHD/脚本/Methscan/01_Upstream
```

## 8. 参考资料

- [MethylVI paper](https://www.nature.com/articles/s42256-026-01225-9)
- [MethylVI model documentation](https://docs.scvi-tools.org/en/latest/user_guide/models/methylvi.html)
- [MethylVI integration tutorial](https://docs.scvi-tools.org/en/latest/tutorials/notebooks/scbs/MethylVI_batch.html)
- [MethylVI reproducibility repository](https://github.com/suinleelab/methylVI-reproducibility)
- [ALLCools documentation](https://lhqing.github.io/ALLCools/intro.html)
- [ALLCools official mCG 5-kb clustering tutorial](https://lhqing.github.io/ALLCools/cell_level/basic/mcg_5kb_basic.html)
- [ALLCools `remove_black_list_region` and `filter_regions` API](https://lhqing.github.io/ALLCools/api/ALLCools/clustering/mcad/index.html)
- [Amemiya, Kundaje & Boyle. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Scientific Reports (2019)](https://www.nature.com/articles/s41598-019-45839-z)
- [ENCODE GRCh38 exclusion list `ENCFF356LFX`](https://www.encodeproject.org/files/ENCFF356LFX/)
- [Liu et al. Single-cell DNA methylome and 3D multi-omic atlas of the adult mouse brain. Nature (2023)](https://www.nature.com/articles/s41586-023-06805-y)
- [Amethyst: single-cell DNA methylation analysis with hg38 ENCODE blacklist filtering (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12521561/)
- [Beyond Blacklists: a critical assessment of exclusion-set generation strategies (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11839099/)
