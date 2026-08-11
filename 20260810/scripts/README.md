# MethylVI 5-kb workflow（2026-08-10）

本目录是针对当前项目的 10 个样本（5 个 IR、5 个 NR）整理的 ALLCools/MethylVI 流程。所有新脚本都在本目录内；父目录中的原始资料不被覆盖。

## 设计要点

- `sample_id` 是样本身份，`condition` 是生物学分组，取值严格为 `IR` 或 `NR`。
- 最终分析只使用通过 MethSCAn 300k 细胞 QC 的细胞。这里的 300k 是“每个细胞至少覆盖 300,000 个 CpG 位点”，不是从全基因组挑选 300,000 个特征；现有 MethSCAn QC 同时设置 `max_sites=10,000,000` 和 `min_meth=55`。
- MethylVI 默认按 `sample_id` 做 batch correction（`MVI_BATCH_KEY=sample_id`），适用于本项目以去除样本间批次差异为主要目的的分析。
- 当前 5 个 IR 和 5 个 NR 与 10 个样本完全绑定，因此 sample batch 与 condition 不是统计上独立的因素；校正后的 latent/UMAP 可能削弱真实的 IR/NR 信号。建议同时检查未校正结果，或在有独立技术批次列时改用 `MVI_BATCH_KEY=technical_batch`。
- `MVI_SAMPLE_ID_REGEX` 已按 MethSCAn 的实际命名设置，可从 `25110891_IR01_Met__barcode`、`IR01__barcode` 或 `IR01_cell123` 提取 `IR01`。如果实际 cell ID 不同，请修改正则，或在注释文件中提供逐细胞 `sample_id`。

## 文件作用

| 顺序 | 文件 | 类型 | 作用 |
|---:|---|---|---|
| 00 | `mvi_00_config.sh` | 配置 | 输入/输出路径、10 样本期望数量、IR/NR 数量和模型参数；入口脚本自动加载 |
| 01 | `mvi_01_sample_metadata.tsv` | 输入 | 当前 MethSCAn 项目的实际 10 个样本表（IR01–IR05、NR01–NR05） |
| 02 | `mvi_02_allcools_prepare_5kb_counts.sh` | 可选上游 | dsub 任务执行脚本；优先直接链接原始 ALLC，生成 MCDS 和 5-kb 聚类结果；`.cov.gz` 转换仅作备用 |
| 03 | `mvi_03_allcools_cluster_5kb.py` | 可选上游 | 执行 ALLCools 5-kb 聚类、t-SNE 和 UMAP |
| 04 | `mvi_04_run_pipeline.sh` | 主入口 | 统一调度 `smoke`、`verify`、`original-sample`、`build`、`train` 和 `plots` |
| 05 | `mvi_05_smoke_test.py` | 主流程 | 用合成数据进行两轮快速 API/环境检查 |
| 06 | `mvi_06_verify_inputs.py` | 主流程 | 审计 H5AD、ALLC、5-kb 坐标、样本元数据和 IR/NR 数量；生成 `mvi_06_input_audit.json` |
| 07 | `mvi_07_plot_original_embedding.py` | 主流程 | 在原始 ALLCools UMAP/t-SNE 上绘制 SCANPY cell type、sample 与 condition |
| 08 | `mvi_08_build_input.py` | 主流程 | 从 ALLC 文件重建 MethylVI 所需的整数 `mc/cov` 层，生成 H5MU |
| 09 | `mvi_09_train_model.py` | 主流程 | 训练 MethylVI，计算 latent、邻居图、UMAP 和 Leiden |
| 10 | `mvi_10_plot_celltype_sample.py` | 主流程 | 绘制训练后 UMAP 的 cell type 与 sample |
| 11 | `mvi_11_plot_condition.py` | 主流程 | 绘制训练后 UMAP 的 IR/NR condition |
| 12 | `mvi_12_run_pipeline.sh` | 集群入口 | 作为 dsub 任务的执行脚本，运行 `mvi_04_run_pipeline.sh all` |
| 13 | `mvi_13_utils_pipeline.py` | 支撑模块 | ID 标准化、样本元数据合并、区域解析、ALLC 聚合和 JSON 输出公共函数 |
| 14 | `mvi_14_utils_plot.py` | 支撑模块 | 分类嵌入绘图公共函数 |
| 15 | `mvi_15_requirements.lock.txt` | 环境记录 | 核心软件版本记录 |
| 16 | `mvi_16_test_utils_pipeline.py` | 测试 | 不接触真实数据的轻量单元测试 |
| 参考 | `hg38.canonical.chrom.sizes` | 参考数据 | hg38 canonical 染色体长度；供 ALLCools 生成 5-kb MCDS 使用 |

## 输入配置

编辑 `mvi_00_config.sh`，至少设置：

```bash
export MVI_DATA_ROOT=/path/to/your/project
export MVI_METHSCAN_UPSTREAM="/Users/luozhixiong/Library/Mobile Documents/com~apple~CloudDocs/Documents/PHD/脚本/Methscan/01_Upstream"
export MVI_ALLCOOLS_OUTPUT=$MVI_DATA_ROOT/methylvi_5kb_300k
export MVI_H5AD=$MVI_ALLCOOLS_OUTPUT/mcg_5kb.clustered.h5ad
export MVI_CHROM_SIZES=/path/to/20260810/scripts/hg38.canonical.chrom.sizes
export MVI_ALLC_DIR=$MVI_ALLCOOLS_OUTPUT/input_allc
export MVI_SAMPLE_METADATA=/path/to/20260810/scripts/mvi_01_sample_metadata.tsv
export MVI_ROOT=$MVI_DATA_ROOT/methylVI_results_300k
```

`mvi_01_sample_metadata.tsv` 必须有唯一的两列：

```text
sample_id  condition
IR01       IR
...
IR05       IR
NR01       NR
...
NR05       NR
```

实际文件可用 TSV 或 CSV。当前 `MVI_ANNOTATION` 默认读取 SCANPY 的全细胞注释表；公共读取器会将 SCANPY 的 `sample`、`group`、`cell_type_integrated` 映射为 MethylVI 内部的 `sample_id`、`condition`、`cell_type`，并保留原始列以便追溯。若同时提供 sample metadata，`condition` 以 sample metadata 为准。

所有新生成的 PDF/PNG 统一保存到 `${MVI_FIGURES_DIR}`。默认为当前项目的 `result` 目录：服务器上是 `/share/home/rzli/MethylVI/20260810/result`，本地则是 `20260810/result`。其下只设两个编号子目录：`${MVI_FIGURES_BEFORE_DIR}` 即 `01_before_methylvi/` 存放ALLCools和校正前图，`${MVI_FIGURES_AFTER_DIR}` 即 `02_after_methylvi/` 存放MethylVI校正后图。模型、H5MU、latent、审计和训练表继续保存在各自的数据结果目录，不与图像混放。

## 服务器验证记录

截至 2026-08-10，服务器上已确认的目录如下：

| 内容 | 服务器绝对路径 |
|---|---|
| MethylVI 脚本工作目录 | `/share/home/rzli/MethylVI/20260810/scripts` |
| MethSCAn 数据根目录 | `/share/LCZX_Data/data/allcools` |
| 原始 ALLCools 样本目录 | `/share/LCZX_Data/data/allcools/25110891_<IR/NR>_Met/allcools` |
| 10 个样本的 `.cov.gz` 目录 | `/share/LCZX_Data/data/allcools/25110891_{IR01..IR05,NR01..NR05}_Met/cov_dedup_probability` |
| hg38 canonical 染色体长度文件 | `/share/home/rzli/MethylVI/20260810/scripts/hg38.canonical.chrom.sizes` |
| 58,534 个原始细胞的输入审计目录 | `/share/LCZX_Data/data/allcools/methylvi_5kb` |
| 最终 300k QC 细胞的 ALLCools 5-kb 输出目录 | `/share/LCZX_Data/data/allcools/methylvi_5kb_300k` |

已确认 10 个样本目录全部存在，`cov_dedup_probability` 中约有 5.85 万个细胞文件，IR 和 NR 两组均包含 5 个样本。

服务器端此前已确认文件编号从 `mvi_00` 连续到 `mvi_16`；重复的 sample metadata example 已删除，实际运行只读取 `mvi_01_sample_metadata.tsv`。本地现已进一步把两个集群脚本统一为普通 `.sh`：`mvi_02_allcools_prepare_5kb_counts.sh` 和 `mvi_12_run_pipeline.sh`，这两个新文件及本 README 需要重新同步到服务器。

服务器上的 `hg38.canonical.chrom.sizes` 已确认存在，共 24 行，包含 canonical chromosomes。登录 shell 默认未初始化 Conda，系统 `/usr/bin/python` 不包含 ALLCools；加载 Conda 后检查了 `methylvi`、`scDNAm`、`scDNAm_v2`、`scWGBS` 和 `scanpy310`，均未安装 ALLCools。`scDNAm` 仅提供 `bgzip` 和 `tabix`，其余候选环境也没有完整工具链。因此 ALLCools 上游需要建立独立的 `allcools` Conda 环境，避免修改现有 MethylVI 和 MethSCAn 环境。

2026-08-10 已在 `linux-aarch64` 服务器上成功创建 `/share/home/rzli/miniconda3/envs/allcools`，指定 Python 3.8，并通过 `conda-forge`、`bioconda` 和 `defaults` 安装 ALLCools 的基础依赖。回传日志已显示 `Executing transaction: done` 并返回 shell 提示符，因此 Conda 环境创建已确认完成。

ALLCools 本体安装后，`python`、`allcools`、`bgzip` 和 `tabix` 均已从 `/share/home/rzli/miniconda3/envs/allcools/bin/` 解析到正确命令。首次 Python 导入验收时，OpenBLAS 根据登录节点的 128 个 CPU 尝试创建 128 个线程，超过账号当前的 `RLIMIT_NPROC` 余量，导致 `pthread_create failed` 并在用户中断后出现 `KeyboardInterrupt`。这不是缺少 Python 包；登录节点验收时应先把 `OPENBLAS_NUM_THREADS`、`OMP_NUM_THREADS`、`MKL_NUM_THREADS` 和 `NUMEXPR_NUM_THREADS` 设为 1。ALLCools 与 MethylVI 的 dsub 执行脚本也将每个工作进程的这些变量固定为 1，并通过 ALLCools 的 `--cpu`/`n_jobs`、Python 进程池或 PyTorch 显式使用 `MVI_THREADS=32` 做任务级并行，避免嵌套超额并行。

限制登录节点线程后，`python -m pip check` 已通过，`allcools --help` 也能正常显示 `generate-dataset` 等子命令。首次验收发现 `scanpy 1.10.2` 使用了 Python 3.8 无法解析的语法，在 `scanpy/readwrite.py` 报 `SyntaxError`。将 Scanpy 固定为明确支持 Python 3.8 的 `1.9.8` 后，导入验收已全部通过。当前已确认的 ALLCools 上游环境为：ALLCools 1.1.1、AnnData 0.9.2、Scanpy 1.9.8、NumPy 1.24.4 和 Pandas 1.5.3；`python -m pip check` 返回 `No broken requirements found`。

ALLCools 1.1.1 API 预检已通过。`MCDS.open`、`MCDS.get_score_adata`、`binarize_matrix`、`filter_regions`、`lsi`、`significant_pc_test`、`tsne` 和 `ConsensusClustering` 均可正常导入，当前 `mvi_03_allcools_cluster_5kb.py` 使用的参数均存在。`allcools generate-dataset --help` 也确认 `--allc_table`、`--output_path`、`--chrom_size_path`、`--obs_dim`、`--cpu`、`--chunk_size`、`--regions` 和可重复的 `--quantifiers` 语法与 `mvi_02_allcools_prepare_5kb_counts.sh` 一致。

2026-08-10 在 `node-4` 尝试用 `sbatch` 提交 `MVI_STAGE_ONLY=1` 预检时，shell 返回 `sbatch: command not found`。失败发生在02脚本启动之前，因此未创建预检输出。用户随后确认该服务器的大任务使用 `dsub` 提交，不使用 Slurm。根据现有 MethSCAn 成功任务，标准申请格式为 `-R "cpu=32;mem=194560MB"`，通过 `--cwd`、`-oo` 和 `-eo` 设置工作目录与日志，并用 `djob` 监控。stage-only 虽然不计算 MCDS，但会遍历并链接 58,534 个细胞文件，也统一使用较小资源的 `dsub` 任务；全量 MCDS 和 MethylVI 使用 32 CPU、190 GiB 的 `dsub` 任务。

2026-08-10 17:26，stage-only 预检已通过 `dsub` 成功提交，任务 ID 为 `163899`，任务名为 `methylvi_stage_allc`，申请资源为 4 CPU、16 GiB。`djob` 显示任务状态为 `RUNNING`，执行节点为 `node-11`，说明当前 `dsub` 参数、工作目录及任务启动方式有效。此记录只表示任务已成功启动；仍需等待任务结束，并从标准输出确认输入模式为 `original_allc`、细胞数为 58,534、样本数为 10、IR=5、NR=5，以及出现 `MVI_STAGE_ONLY=1` 完成信息，之后才能提交完整 MCDS 任务。

任务运行期间再次检查：标准输出已写入 `host=node-11`、输入 `/share/LCZX_Data/data/allcools`、输出 `/share/LCZX_Data/data/allcools/methylvi_5kb` 和 `threads=4`；错误日志仍为 0 字节。02脚本在扫描原始 ALLC 并逐一建立 ALLC/TBI 软链接期间不会逐文件打印进度，因此这一阶段主日志只有启动行属于预期行为，可通过 `source_allc_manifest.tsv.tmp.*` 的行数和 `input_allc/` 中软链接数观察进度。

后续进度检查显示独立 ALLCools 环境记录完整：Python 3.8.20、ALLCools 1.1.1、AnnData 0.9.2、NumPy 1.24.4、Pandas 1.5.3、Scanpy 1.9.8、SciPy 1.10.1 和 scikit-learn 1.3.2。任务 `163899` 仍在 `node-11` 运行，临时 source manifest 已写入 22,736 行，平铺 ALLC 软链接已建立 22,872 个，错误日志为空。运行中链接数略高于 manifest 行数是因为脚本先建立 ALLC/TBI 链接、再追加 manifest 记录，属于正常的瞬时进度差；任务完成时二者应统一为 58,534。

随后正式 `source_allc_manifest.tsv` 已生成且为 58,534 行，`input_allc/` 中的 `*.allc.tsv.gz` 软链接也为 58,534 个。标准输出明确记录 `using original indexed ALLC files: 58534 cells`，证明脚本使用的是原始带 tabix 索引的 ALLC，没有进入 cov fallback；错误日志仍为空。此时任务仍为 `RUNNING`，正在逐细胞核对 ALLC/TBI 链接并生成 `selected_cells.allc.tsv`，必须继续等待最终的 `verified mode=original_allc cells=58534 samples=10 IR=5 NR=5` 和 stage-only 完成信息。

stage-only 任务 `163899` 最终成功结束，`djob` 已无活动匹配，调度器报告 `EXIT_CODE: 0` 和 `Job execution succeeded`。正式 `selected_cells.allc.tsv` 为 58,534 行，ALLC 软链接为 58,534 个，TBI 软链接也为 58,534 个，因此原始输入整理与索引完整性检查通过。业务日志在 17:33:24 明确记录 `verified mode=original_allc cells=58534 samples=10 IR=5 NR=5` 和 `MVI_STAGE_ONLY=1; ALLC staging and verification completed, skipping MCDS generation`，且没有 `ERROR` 或 `Traceback`。该任务实际运行 400 秒，申请 4 CPU/16 GiB，实测内存峰值 321 MiB；错误日志为空。检查还确认当前不存在 `mcds.COMPLETE` 或 `mcg_5kb.clustered.h5ad`，并且没有其他活动任务，因此可以进入完整 MCDS 生成和 5-kb 聚类步骤。

在准备提交完整 MCDS 前，进一步确认最终分析需要沿用 MethSCAn 的 300k 细胞 QC。MethSCAn 工作流定义为 `min_sites=300000`、`max_sites=10000000`、`min_meth=55`、`max_meth=none`，每个样本通过细胞记录在 `<sample>/qc_minmeth55_maxmethnone_maxsites10000000_covdedupprob/filtered_data_single_300k/column_header.txt`。服务器逐样本核查结果为：IR01 583、IR02 965、IR03 914、IR04 318、IR05 439、NR01 306、NR02 755、NR03 658、NR04 722、NR05 539；IR 合计 3,219、NR 合计 2,980，总计 6,199。10个 `filter_provenance.tsv` 的四项阈值全部一致，抽查 header 为每行一个原始 barcode。上述 58,534 细胞 stage-only 因此只作为原始 ALLC 完整性审计，不能直接用于最终 MCDS。02脚本已改为强制按这10个 header 白名单选择原始 ALLC，期望细胞数改为 6,199，并使用独立的 `methylvi_5kb_300k` 输出目录。

修改后的 `mvi_00_config.sh` 和 `mvi_02_allcools_prepare_5kb_counts.sh` 已重新同步服务器并通过 `bash -n`。服务器 grep 已确认默认输出为 `methylvi_5kb_300k`、`MVI_EXPECTED_CELLS=6199`、`MVI_USE_FILTERED_CELLS=1`，02脚本包含白名单加载、原始 ALLC/cov fallback 白名单过滤和最终细胞数硬检查。新的 `/share/LCZX_Data/data/allcools/methylvi_5kb_300k` 尚不存在，说明不会覆盖或混入此前 58,534 细胞的审计目录，可以提交300k stage-only 预检。

2026-08-10 17:45，300k stage-only 预检已通过 `dsub` 成功提交，任务 ID 为 `163900`，任务名为 `methylvi_stage_300k`，申请4 CPU/16 GiB，`djob` 显示在 `node-11` 运行。该任务使用独立输出 `/share/LCZX_Data/data/allcools/methylvi_5kb_300k`；需等待日志确认白名单6,199个、原始 ALLC 6,199个、样本10个及 IR/NR=5/5后，才能提交完整 MCDS。

任务 `163900` 的早期日志已确认 `loaded MethSCAn 300k whitelist: cells=6199 samples=10`，错误日志为空。正式 `filtered_cell_whitelist.tsv` 为6,199行，`filtered_qc_summary.tsv` 包含 IR01 583、IR02 965、IR03 914、IR04 318、IR05 439、NR01 306、NR02 755、NR03 658、NR04 722、NR05 539，且每行均记录实际 header 和 provenance 绝对路径。白名单加载与逐样本数量检查通过，任务继续扫描全部原始 ALLC并匹配这6,199个 barcode。

任务 `163900` 最终成功完成：`source_allc_manifest.tsv`、ALLC软链接和TBI软链接均为6,199；业务日志依次记录 `loaded MethSCAn 300k whitelist: cells=6199 samples=10`、`using original indexed ALLC files: 6199 cells`、`verified mode=original_allc cells=6199 samples=10 IR=5 NR=5` 和 stage-only 完成信息。调度器报告 `EXIT_CODE: 0`、`Job execution succeeded`，运行134秒，内存峰值135 MiB。标准错误文件末尾出现的内容仅为调度器追加的成功状态与资源摘要，不是程序错误。300k白名单、原始ALLC和索引已全部验收，可以正式生成MCDS及5-kb聚类结果。

2026-08-10 17:48，完整 ALLCools 任务已通过 `dsub` 成功提交，任务 ID 为 `163901`，任务名为 `methylvi_allcools_5kb_300k`，申请32 CPU、190 GiB，`djob` 显示在 `node-11` 运行。任务输出仍为独立的 `/share/LCZX_Data/data/allcools/methylvi_5kb_300k`；流程会先复核6,199个白名单细胞，再运行 `allcools generate-dataset` 和03脚本的5-kb聚类、t-SNE及UMAP。

任务 `163901` 的早期检查显示白名单已再次成功加载为6,199个细胞，错误日志为空；此时任务仍在扫描全部58,534个原始 ALLC进行白名单匹配，尚未创建 `mcg_5kb.mcds`。上一轮同一匹配阶段约需2分钟，因此当前状态正常，需等日志出现最终 `verified mode` 后再判断是否进入 `generate-dataset`。

后续检查时任务仍为 `RUNNING`，日志尚未新增 `using original` 或 `verified mode`，输出目录也尚未出现 `mcg_5kb.mcds`；标准错误中没有 `ERROR` 或 `Traceback`。02脚本只有完成整轮原始文件扫描和软链接核对后才打印下一条业务日志，因此需结合本轮 `source_allc_manifest.tsv.tmp.*` 的行数与修改时间判断内部进度，不能仅凭主日志暂时不增长认定任务卡死。

2026-08-11 09:02 复查时，任务 `163901` 已不在 `djob` 活动列表。业务日志确认6,199个300k白名单细胞再次通过原始 ALLC、样本数和IR/NR核验，并在日志第635行记录 `generate-dataset finished.`。输出目录已存在3.8 GiB的 `mcg_5kb.mcds` 及 `mcds.COMPLETE`，且没有残留的临时 source manifest，因此MCDS生成阶段已完成。调度日志更新时间为2026-08-10 22:58；尚需核查任务最终退出码、日志末尾及 `mcg_5kb.clustered.h5ad`，才能判断后续03聚类、t-SNE和UMAP是否也成功，当前暂不进入MethylVI训练。

完整日志与结果验收确认任务 `163901` 全部成功，不需要重跑。任务输出 `ALLCools mCG-5kb clustering completed`，调度器报告 `EXIT_CODE: 0` 和 `Job execution succeeded`。初始矩阵为6,199细胞×617,665个5-kb bins，区域过滤后为6,199×231,648，显著LSI成分为8个；500次Leiden得到9个原始簇和21个临时outlier，监督合并在4个非outlier标签时停止，最终将所有outlier重新分配，十折交叉验证准确率为0.964。已生成77 MiB的 `mcg_5kb.clustered.h5ad`、61 KiB的 `cell_clusters.csv.gz`、1.3 MiB的 `tsne.L1.png` 和1.7 MiB的 `umap.L1.png`。任务运行18,035秒（约5小时），内存峰值20,745 MiB，平均CPU利用率约28.77核。stderr中的Numba FNV hashing和Scanpy colormap信息均为 `UserWarning`，不影响结果。下一步先验收H5AD内部维度、cell ID、聚类列和embedding，再运行MethylVI smoke/verify。

2026-08-11 首次在新的登录shell中运行H5AD只读验收时，未预先限制BLAS线程，OpenBLAS按节点128个CPU尝试创建128线程并超过账号 `RLIMIT_NPROC`，在导入NumPy/h5py阶段出现 `pthread_create failed`，随后由用户中断并显示 `KeyboardInterrupt`。此次失败发生在 `anndata.read_h5ad` 之前，H5AD未被打开或修改，不属于文件损坏。每次新登录会话运行Python轻量检查前都需先导出 `OPENBLAS_NUM_THREADS=1`、`OMP_NUM_THREADS=1`、`MKL_NUM_THREADS=1`、`NUMEXPR_NUM_THREADS=1` 和 `NUMBA_NUM_THREADS=1`；dsub执行脚本本身已包含相应限制。

设置上述5个线程变量为1后，H5AD只读验收全部通过：shape为6,199×231,648，obs名称唯一，样本计数与10份300k白名单逐一一致；obs包含 `leiden`、`L1`、`L1_proba`，obsm包含 `X_pca`、`X_tsne`、`X_umap`，t-SNE和UMAP均为6,199×2。最终 `L1` 包含4簇：c0 3,600、c1 1,398、c2 665、c3 536。所有断言通过并输出 `H5AD 内容验收通过`，可以进入MethylVI环境及smoke test验收。

MethylVI环境首次验收发现两个问题：`python -m pip check` 报告 `nvidia-cusparselt-cu13 0.8.1 is not supported on this platform`，说明当前linux-aarch64环境残留一个不适配平台的NVIDIA依赖；核心导入则因本地05/09脚本使用 `from methyl_vi.model import MethylVI` 而报 `ModuleNotFoundError: methyl_vi`。根据MethylVI官方复现仓库和scvi-tools 1.3.3文档，参考实现位于 `scvi.external.METHYLVI`，不是必须另装的独立 `methyl_vi` 包；官方 `setup_mudata` 参数为 `methylation_contexts` 和 `modalities`。因此当前不能直接安装所谓 `methyl-vi==0.1.0`，必须先读取服务器已装scvi-tools的真实类与函数签名，再修正05、09和版本锁文件并重新运行smoke test。

服务器真实API检查确认架构为aarch64、scvi-tools 1.3.3、Torch运行时2.12.1+cu130，当前节点 `CUDA available=False`、GPU数为0；`scvi.external.METHYLVI` 可正常导入。官方1.3.3源码确认 `setup_mudata` 使用 `methylation_contexts=["mCG"]` 和 `modalities={"batch_key": "mCG"}`，底层模型的 `likelihood` 可选 `betabinomial/binomial`，`dispersion` 可选 `region/region-cell`。本地05/09脚本已按此修正，移除错误的独立包导入和 `gene/gene-cell` 参数；15版本文件也改为记录scvi-tools内置实现及服务器实际Torch版本。Python编译检查和6个单元测试均通过，需重新同步05、09、15和README后运行CPU smoke test。当前CUDA依赖警告先保留观察，不在API与smoke验证前卸载Torch依赖。

修正后的05/09/15已同步服务器，CPU smoke test成功完成。合成数据为32细胞×64特征，使用 `sample_id` 作为batch，官方METHYLVI模型在CPU训练2轮后正常停止，latent维度为32×4且数值有效，最终输出 `MethylVI smoke test passed: latent shape (32, 4)`。MuData 0.4行为变化、Lightning pytree、DataLoader worker数和logging interval信息均为FutureWarning/性能提示，不是失败。由于当前节点无GPU，后续正式训练需先评估CPU可行性或寻找GPU队列；在此之前先运行真实输入verify及计数构建审计。

在运行真实输入verify前，06脚本进一步加固：除原有样本数、IR/NR=5/5、cell ID、5-kb坐标、缺失ALLC和未知分组检查外，现在拒绝ALLC目录中任何不属于H5AD的额外细胞，并在审计完成后显式关闭backed H5AD句柄。H5AD细胞数不固定写死为6,199；6,199只是当前300k QC数据的实测值，通用审计以H5AD与当前ALLC白名单一一对应为准，便于以后更换数据集或QC阈值。修正后Python编译检查和6个单元测试通过，需重新同步06和README后再运行verify。

真实输入verify已成功通过并生成 `mvi_06_input_audit.json`：H5AD为6,199细胞×231,648个保留5-kb区域，全部区域坐标解析成功；10个样本均被识别，无缺失样本，IR为3,219细胞、NR为2,980细胞；ALLC目录正好6,199个文件，缺失和额外ALLC均为0，未知sample和condition均为0。H5AD不含计数layer属于预期，因为其中X只作为ALLCools聚类分数，MethylVI所需整数 `mc/cov` 会从原始ALLC重建。预计两层稠密uint16矩阵占5.35 GiB。下一步先生成校正前ALLCools嵌入对照图，再通过dsub执行计数build。

已在服务器确认 SCANPY 正式逐细胞注释文件为 `/share/home/rzli/SCANPY/20260714/result/annotation/02_cell_annotation_all_cells.csv`，文件为3.4 MiB，共56,747行（包含表头，即56,746个细胞）。字段为 `cell_id,sample,group,batch,leiden_integrated,cell_type_integrated,exclude_from_main_analysis,analysis_status`；例如 `IR01_GAGGTGTATTTGGTGAG` 注释为 `CD14_Monocytes`。读取器会将 SCANPY 的 `IR01_<barcode>` 与 ALLCools 的 `IR01__<barcode>` 标准化后匹配，不再把 ALLCools `L1` c0–c3 当作真实 cell type。07脚本已补充 cell type、sample、condition 三类校正前嵌入图，并显式关闭backed H5AD句柄。

2026-08-11 实际交集核验显示：当前6,199个300k甲基化QC细胞中，5,765个（93.0%）能唯一匹配 SCANPY 注释，434个（7.0%）不在 SCANPY 注释对象中，标准化后注释ID重复数为0。各样本匹配数为 IR01 562、IR02 880、IR03 852、IR04 300、IR05 408、NR01 283、NR02 699、NR03 630、NR04 642、NR05 509。匹配细胞涵盖16个 cell type；其中5,748个为 `Keep`，17个 `Platelet_erythroid_contamination` 被标记为 `Exclude`。MethylVI建模保留全部6,199个通过甲基化QC的细胞；434个未匹配细胞记为 `Unannotated`，避免用RNA/SCANPY是否保留反向改变甲基化建模队列。17个 `Exclude` 细胞的原始标记保留在H5MU元数据中，细胞类型下游统计时可明确排除，不影响batch correction模型使用所有合格甲基化细胞。

同日首次在服务器检查新注释代码时，`mvi_00_config.sh` 仍显示空的 `MVI_ANNOTATION`，`mvi_13_utils_pipeline.py` 也未检出 `cell_type_integrated`、`annotation_match_rate` 或 `scanpy_excluded` 字段，说明服务器尚未同步本地新版文件。当时运行的6个旧版单元测试全部通过，但不包含新增的 SCANPY/ALLCools cell ID 兼容和列别名测试；新版应显示8个测试。因此必须先同步00、07、08、13、16和README，再重新测试与verify。

新版文件同步后，服务器已能正确检出 SCANPY 路径、列别名和新增的8个测试；其中7个通过，`test_scanpy_annotation_column_aliases` 在将 `condition` 转大写时失败。原因是服务器 NumPy 对 `object` dtype 数组调用 `np.char.upper()` 会报 `TypeError: string operation on non-string array`，与数据文件无关。公共读取器已改为先在 Pandas 字符串列上调用 `.str.upper()`，再显式转为字符串数组；需重新同步13和README并确认8个测试全部通过。

修正后服务器8个单元测试已全部通过，真实输入 `verify` 也成功完成。审计确认H5AD为6,199细胞×231,648个保留5-kb bins，6,199个ALLC全部匹配，缺失和额外ALLC均为0，10个样本和IR/NR分组均无未知值。SCANPY注释表共56,746个细胞，当前队列匹配5,765个、未匹配434个，匹配率为0.9299887078561059；5,748个匹配细胞为SCANPY `Keep`，17个为`Exclude`。审计报告已列出16个正式cell type和434个`Unannotated`，预计稠密uint16 `mc/cov` 两层占5.35 GiB。输入已满足生成校正前嵌入图和后续build的条件。

校正前嵌入绘图已成功完成，所有图像已统一保存在 `/share/home/rzli/MethylVI/20260810/result`。当前共5个文件：`allcools_original_embedding_cell_type.pdf` 122 KiB、`allcools_original_embedding_condition.pdf` 110 KiB、`allcools_original_embedding_sample_id.pdf` 116 KiB、`allcools_5kb_tsne_L1.png` 1.3 MiB和 `allcools_5kb_umap_L1.png` 1.7 MiB。后两张ALLCools图已从数据输出目录迁移至项目 `result`，没有保留重复副本；后续MethylVI UMAP PDF也使用同一目录。校正前对照阶段已完成，下一步为通过dsub重建6,199细胞×231,648区域的整数 `mc/cov` 计数并生成H5MU。

build提交前复核发现，旧的主入口会将 `NUMBA_NUM_THREADS` 设为 `MVI_THREADS`，而08脚本同时使用 `MVI_THREADS` 个多进程聚合ALLC，存在每个进程再创建32个Numba线程的过度订阅风险。04脚本已改为：build、verify、smoke和plots阶段的 `NUMBA_NUM_THREADS=1`；仅train阶段启动单个Python进程时，再临时使用 `MVI_THREADS`。这不改变计数结果，只防止build在计算节点上过度创建线程。

图片目录已进一步整理为两个编号子目录。服务器现有的5张图已全部迁移到 `/share/home/rzli/MethylVI/20260810/result/01_before_methylvi`：3张原始嵌入PDF和2张ALLCools L1 PNG；`result/02_after_methylvi` 已创建且当前为空，用于后续MethylVI校正后图。`result` 根目录不再直接存放图片，未留下重复副本。

H5MU build提交前检查已通过。当前 `MVI_ROOT=/share/LCZX_Data/data/allcools/methylVI_results_300k`，目标为其下 `methylvi_5kbin_input.h5mu`；目标H5MU尚不存在，`count_rows` 检查点目录也尚未创建，因此这是无旧结果冲突的首次build。主入口已确认build默认 `NUMBA_NUM_THREADS=1`，仅train临时使用 `MVI_THREADS`。`/share/LCZX_Data` 文件系统容量50 TiB，已用25 TiB，可用26 TiB；`/share/home` 可用573 TiB；当前MethylVI输出根仅360 KiB，空间充足。

2026-08-11 10:22:23，首次H5MU build已通过dsub成功提交，任务ID为 `164085`，任务名为 `methylvi_build_300k`，申请32 CPU和194,560 MiB内存，当前在 `node-11` 上处于 `RUNNING` 状态。调度器标准输出和错误日志分别为 `scheduler_logs/methylvi_build_300k.164085.out` 和 `.err`，脚本内部进度日志为 `logs/mvi_08_build_input.log`。

任务启动日志已输出 `aggregating 6,199 cells x 231,648 retained regions with 32 workers`，stderr为空，djob仍显示在node-11运行。首次监控时 `count_rows` 中已出现数十个检查点；`du` 扫描期间有一个 `.tmp.npz` 消失，是工作进程完成写入后将临时文件原子重命名为正式 `.npz` 导致的正常并发竞态，不是文件损坏。后续监控将分开统计正式检查点和 `.tmp.npz`。

后续监控时已生成242个正式 `.npz` 检查点，抽查时临时 `.tmp.npz` 为0；内部日志已连续记录 `count rows 50/6,199`、100、150和200，均为 `built`、`reused=0`。调度器stderr仍为空，任务 `164085` 在node-11保持 `RUNNING`，说明32进程聚合和检查点写入均正常。

临时30秒刷新监控显示任务已达到658/6,199个正式检查点（10.61%），日志已连续到 `count rows 650/6,199 (built=650, reused=0)`，与文件计数一致。任务仍在node-11运行，当前阶段仍为逐细胞ALLC聚合，未进入稠密矩阵组装。

服务器上尚未发现 H5AD 或 H5MU，但原始逐细胞 ALLC 已确认存在。每个样本在 `/share/LCZX_Data/data/allcools/25110891_<sample>_Met/allcools` 下包含多个按 barcode 前缀拆分的 `*_merged_fr_bam_allcools` 子目录，逐细胞文件命名为 `<barcode>_allc.gz`。各样本 ALLC 数量为 IR01 7,981、IR02 6,070、IR03 7,383、IR04 8,171、IR05 5,392、NR01 4,340、NR02 5,672、NR03 4,285、NR04 7,057、NR05 2,183，总计 58,534，与 `cov_dedup_probability` 的细胞数完全一致。IR01 样本目录下同时存在 `total` 和 `allcools` 目录。之前检索只匹配 `*.allc.tsv.gz`，因此漏掉了实际的 `*_allc.gz` 命名和更深的目录层级。

上游主路线已确定直接从 58,534 个原始 ALLC 中按 MethSCAn 300k 白名单选择 6,199 个细胞，不从 `.cov.gz` 反向重建。IR01 的 `total` 已确认为包含 7,981 个 ALLC 的软链接汇总目录；抽样 ALLC 为标准 7 列格式（chrom、1-based position、strand、context、mc、cov、call），`gzip -t` 和 tabix 查询均通过。10 个样本的 ALLC 和 `.tbi` 数量逐一相等，合计均为 58,534。`mvi_02_allcools_prepare_5kb_counts.sh` 会验证每个样本的 filter provenance，加载 6,199 个 `sample__barcode` 白名单，扫描嵌套的原始 `<barcode>_allc.gz`，仅为白名单细胞建立 ALLC/TBI 软链接，并硬性检查细胞总数 6,199、样本数 10 及 IR/NR=5/5。在上游生成 H5AD 之前，仍不能运行 `verify`、`build` 或 `train`。

服务器运行结果、核查结论、使用的绝对路径和后续修正应持续补充到本节。

## 与当前 MethSCAn 上游数据的对应关系

已检查 `/Users/luozhixiong/Library/Mobile Documents/com~apple~CloudDocs/Documents/PHD/脚本/Methscan/01_Upstream` 中的工作流说明和 `Report.md`。其中明确的样本目录名为：

```text
25110891_IR01_Met  ...  25110891_IR05_Met
25110891_NR01_Met  ...  25110891_NR05_Met
```

每个样本的数据层级为：

```text
<sample>/cov
  → <sample>/cov_dedup_probability
  → <sample>/compact_data_dedup_probability
  → <sample>/qc_minmeth55_maxmethnone_maxsites10000000_covdedupprob/filtered_data_single_300k
```

当前 MethylVI/ALLCools 上游优先使用每个样本 `allcools` 子目录中的原始 ALLC；只有原始 ALLC 不可用时，才使用 `cov_dedup_probability` 作为备用输入，不使用未去重的原始 `cov`。`mvi_01_sample_metadata.tsv` 已按 `IR01–IR05` 和 `NR01–NR05` 建好；上游脚本会将 cell ID 统一为 `IR01__barcode`/`NR01__barcode`，`mvi_00_config.sh` 中的正则可提取对应 sample ID。服务器数据根目录为 `/share/LCZX_Data/data/allcools`。

已核对 MethSCAn 历史转换逻辑：原始 ALLC 的正链 CpG 使用原 `pos`，负链 CpG 先使用 `pos - 1`，再写成 `start=end` 的 cov 坐标。这说明 cov 已丢失原始链方向信息；在原始 ALLC 可用时，MethylVI/ALLCools 上游应优先直接使用原始 ALLC。`.cov.gz` 转回 ALLC 的逻辑仅作为无原始 ALLC 时的备用方案。

当前服务器应直接把包含 10 个样本目录的数据根目录作为第一个参数；脚本会优先发现原始 ALLC。首次先提交 stage-only 预检：

```bash
cd /share/home/rzli/MethylVI/20260810/scripts
mkdir -p scheduler_logs

dsub \
  -n methylvi_stage_300k \
  -R "cpu=4;mem=16384MB" \
  --cwd /share/home/rzli/MethylVI/20260810/scripts \
  -oo scheduler_logs/methylvi_stage_300k.%J.out \
  -eo scheduler_logs/methylvi_stage_300k.%J.err \
  env MVI_STAGE_ONLY=1 MVI_THREADS=4 \
  bash mvi_02_allcools_prepare_5kb_counts.sh \
  /share/LCZX_Data/data/allcools \
  /share/LCZX_Data/data/allcools/methylvi_5kb_300k \
  /share/home/rzli/MethylVI/20260810/scripts/hg38.canonical.chrom.sizes
```

如果数据根目录下没有原始 ALLC，脚本才会自动回退到嵌套的 `*_Met/cov_dedup_probability/*.cov.gz`，在输出目录建立 `staged_cov/` 并转换为 ALLC。

首次服务器预检时应设置 `MVI_STAGE_ONLY=1`。此模式验证300k QC provenance和白名单，建立6,199个 ALLC/索引软链接，生成 manifest 和 ALLC table，并核验 10 个样本及 5 IR + 5 NR，不会启动耗时的 MCDS 生成和聚类。预检通过后，再用下面的 `dsub` 命令运行完整上游。

可选的 ALLCools 上游步骤使用参数或环境变量提供路径。当前服务器的实际命令为：

```bash
cd /share/home/rzli/MethylVI/20260810/scripts
mkdir -p scheduler_logs

dsub \
  -n methylvi_allcools_5kb_300k \
  -R "cpu=32;mem=194560MB" \
  --cwd /share/home/rzli/MethylVI/20260810/scripts \
  -oo scheduler_logs/methylvi_allcools_5kb_300k.%J.out \
  -eo scheduler_logs/methylvi_allcools_5kb_300k.%J.err \
  env MVI_THREADS=32 MVI_STAGE_ONLY=0 \
  bash mvi_02_allcools_prepare_5kb_counts.sh \
  /share/LCZX_Data/data/allcools \
  /share/LCZX_Data/data/allcools/methylvi_5kb_300k \
  /share/home/rzli/MethylVI/20260810/scripts/hg38.canonical.chrom.sizes
```

等价环境变量为 `MVI_DATA_ROOT`、`MVI_ALLCOOLS_OUTPUT`、`MVI_CHROM_SIZES` 和 `MVI_ALLCOOLS_ENV`；软件命令可用 `ALLCOOLS_EXE`、`PYTHON_BIN`、`BGZIP_EXE` 和 `TABIX_EXE` 覆盖。

## 数据处理流程与脚本对应关系

整个流程分为“可选的 ALLCools 上游准备”和“MethylVI 主流程”两部分。每一步的脚本、输入和输出如下。

| 步骤 | 处理内容 | 使用脚本 | 主要输入 | 主要输出 |
|---|---|---|---|---|
| 0. 环境与参数 | 设置数据目录、H5AD、ALLC 目录、样本元数据和模型参数 | `mvi_00_config.sh` | 配置文件 | 环境变量 |
| 1. ALLCools 上游（可选） | 验证 MethSCAn 300k QC 白名单；从原始 ALLC 选择6,199个细胞；生成 5-kb MCDS；按 mCG hypo-score 做聚类和 UMAP | `mvi_02_allcools_prepare_5kb_counts.sh`、`mvi_03_allcools_cluster_5kb.py` | 58,534 个原始 ALLC、10个 filtered header、染色体长度文件；去重 cov 仅作备用 | `filtered_cell_whitelist.tsv`、`input_allc/`、MCDS、`mcg_5kb.clustered.h5ad`；聚类图写入 `${MVI_FIGURES_DIR}` |
| 2. 依赖与 API 检查 | 用小型合成数据训练两轮，确认 Python、MethylVI、scvi-tools 和 PyTorch 环境可用 | `mvi_05_smoke_test.py` | 无真实数据 | 终端通过/失败状态 |
| 3. 输入审计 | 检查 H5AD 细胞和 5-kb bins、ALLC 文件匹配、SCANPY cell type、样本 ID、IR/NR=5/5、坐标和元数据 | `mvi_06_verify_inputs.py`；公共函数来自 `mvi_13_utils_pipeline.py` | H5AD、ALLC、`mvi_01_sample_metadata.tsv`、SCANPY全细胞注释表 | `${MVI_AUDIT}` |
| 4. 原始嵌入对照 | 在 ALLCools 已有 UMAP/t-SNE 上绘制 SCANPY cell type、sample 和 IR/NR，作为 MethylVI 前的对照 | `mvi_07_plot_original_embedding.py`、`mvi_14_utils_plot.py` | H5AD、SCANPY注释、样本元数据 | `allcools_original_embedding_cell_type.pdf`、`allcools_original_embedding_sample_id.pdf`、`allcools_original_embedding_condition.pdf` |
| 5. 计数重建 | 从原始 ALLC 在保留的 5-kb bins 内重新聚合 mCG 甲基化数和覆盖数；不使用 H5AD 的聚类分数作为 MethylVI 计数 | `mvi_08_build_input.py`；公共函数来自 `mvi_13_utils_pipeline.py` | ALLC、H5AD 中的 bins | `${MVI_INPUT}`、`count_rows/`、`build_summary.json` |
| 6. MethylVI 训练 | 读取 H5MU 的 `mc/cov`，默认按 `sample_id` 做 batch correction；计算 latent、邻居图、UMAP 和 Leiden | `mvi_09_train_model.py` | `${MVI_INPUT}`、`MVI_SAMPLE_KEY`、`MVI_CONDITION_KEY`、`MVI_BATCH_KEY` | 模型、latent、embedding、训练历史、`sample_by_condition.csv` |
| 7. 结果可视化 | 按 cell type、sample 和 condition 绘制训练后 UMAP | `mvi_10_plot_celltype_sample.py`、`mvi_11_plot_condition.py`、`mvi_14_utils_plot.py` | `methylvi_embedding.h5ad` 或 `cell_annotations_umap.tsv.gz` | `methylvi_umap_cell_type.pdf`、`methylvi_umap_sample_id.pdf`、`methylvi_umap_condition.pdf` |

其中第 1 步只有在当前项目还没有合格的 ALLCools 5-kb 结果时才执行；如果已有 `mcg_5kb.clustered.h5ad` 和对应的逐细胞 ALLC 文件，可直接从第 2 步开始。第 3 步必须在第 5、6 步之前通过。第 4 步是校正前对照，不改变输入数据；第 5 步才是真正生成 MethylVI 的 `mc/cov` 计数矩阵。

## 运行顺序

```bash
cd /path/to/20260810/scripts
bash mvi_04_run_pipeline.sh verify
bash mvi_04_run_pipeline.sh smoke
bash mvi_04_run_pipeline.sh original-sample
bash mvi_04_run_pipeline.sh build
bash mvi_04_run_pipeline.sh train
bash mvi_04_run_pipeline.sh plots
```

也可以运行 `bash mvi_04_run_pipeline.sh all`；大任务使用下面的命令提交完整流程。首次运行建议先单独执行 `verify` 和 `smoke`。

```bash
cd /share/home/rzli/MethylVI/20260810/scripts
mkdir -p scheduler_logs

dsub \
  -n methylvi_ir_nr \
  -R "cpu=32;mem=194560MB" \
  --cwd /share/home/rzli/MethylVI/20260810/scripts \
  -oo scheduler_logs/methylvi_ir_nr.%J.out \
  -eo scheduler_logs/methylvi_ir_nr.%J.err \
  env MVI_THREADS=32 \
  bash mvi_12_run_pipeline.sh
```

提交后使用 `djob` 查看任务状态，标准输出和错误日志位于 `scheduler_logs/`。本服务器不使用 `sbatch` 或 `bjobs`。

输入审计会拒绝：样本数不是 10、IR/NR 不是 5/5、10 个元数据样本没有全部出现在选中细胞中、无法从 cell ID 或注释得到 sample、condition 不是 IR/NR、缺少 ALLC 文件、5-kb 坐标无法解析，或坐标与 ALLC 不一致。

## 主要输出

- `${MVI_INPUT}`：包含 `mCG.layers['mc']` 和 `mCG.layers['cov']` 的 H5MU；
- `${MVI_RESULTS}/model/`：训练后的 MethylVI 模型；
- `${MVI_RESULTS}/latent_representation.npy`：latent 表示；
- `${MVI_RESULTS}/methylvi_embedding.h5ad`：邻居图、UMAP 和 Leiden；
- `${MVI_RESULTS}/cell_annotations_umap.tsv.gz`：UMAP、sample、condition、cell type 和 latent；
- `${MVI_RESULTS}/sample_by_condition.csv`：样本与 IR/NR 的细胞数交叉表；
- `${MVI_FIGURES_BEFORE_DIR}`：ALLCools聚类图和MethylVI校正前的cell type/sample/condition对照图；
- `${MVI_FIGURES_AFTER_DIR}`：MethylVI校正后的cell type/sample/condition UMAP；
- `${MVI_RESULTS}/run_summary.json`、`${MVI_RESULTS}/training_history.csv`：参数、软件版本和训练记录；
- `${MVI_AUDIT}`：训练前输入审计报告。

## 断点与安全

`mvi_08_build_input.py` 按细胞写入 `count_rows/*.npz` 检查点，并用 manifest 防止不同输入被错误复用。已有形状和图层正确的 H5MU 会被验证后复用；需要强制重组时才使用 `--force-assemble`。不要让两个 build/train 作业同时写入相同的 `MVI_ROOT`。

## 参考资料

- [MethylVI paper](https://www.nature.com/articles/s42256-026-01225-9)
- [MethylVI model documentation](https://docs.scvi-tools.org/en/latest/user_guide/models/methylvi.html)
- [MethylVI integration tutorial](https://docs.scvi-tools.org/en/latest/tutorials/notebooks/scbs/MethylVI_batch.html)
- [MethylVI reproducibility repository](https://github.com/suinleelab/methylVI-reproducibility)
- [ALLCools documentation](https://lhqing.github.io/ALLCools/intro.html)
