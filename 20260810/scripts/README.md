# MethylVI 5-kb workflow（2026-08-10）

本目录是针对当前项目的 10 个样本（5 个 IR、5 个 NR）整理的 ALLCools/MethylVI 流程。所有新脚本都在本目录内；父目录中的原始资料不被覆盖。

## 设计要点

- `sample_id` 是样本身份，`condition` 是生物学分组，取值严格为 `IR` 或 `NR`。
- MethylVI 默认按 `sample_id` 做 batch correction（`MVI_BATCH_KEY=sample_id`），适用于本项目以去除样本间批次差异为主要目的的分析。
- 当前 5 个 IR 和 5 个 NR 与 10 个样本完全绑定，因此 sample batch 与 condition 不是统计上独立的因素；校正后的 latent/UMAP 可能削弱真实的 IR/NR 信号。建议同时检查未校正结果，或在有独立技术批次列时改用 `MVI_BATCH_KEY=technical_batch`。
- `MVI_SAMPLE_ID_REGEX` 已按 MethSCAn 的实际命名设置，可从 `25110891_IR01_Met__barcode`、`IR01__barcode` 或 `IR01_cell123` 提取 `IR01`。如果实际 cell ID 不同，请修改正则，或在注释文件中提供逐细胞 `sample_id`。

## 文件作用

| 顺序 | 文件 | 类型 | 作用 |
|---:|---|---|---|
| 00 | `mvi_00_config.sh` | 配置 | 输入/输出路径、10 样本期望数量、IR/NR 数量和模型参数；入口脚本自动加载 |
| 01 | `mvi_01_sample_metadata.tsv` | 输入 | 当前 MethSCAn 项目的实际 10 个样本表（IR01–IR05、NR01–NR05） |
| 02 | `mvi_02_sample_metadata.tsv.example` | 模板 | 5 IR + 5 NR 的元数据模板；样本命名变化时据此修改 |
| 03 | `mvi_03_allcools_prepare_5kb_counts.sbatch` | 可选上游 | 从 `.cov`/`.cov.gz` 生成 ALLC、MCDS 和 5-kb 聚类结果；路径通过参数或环境变量提供 |
| 04 | `mvi_04_allcools_cluster_5kb.py` | 可选上游 | 执行 ALLCools 5-kb 聚类、t-SNE 和 UMAP |
| 05 | `mvi_05_run_pipeline.sh` | 主入口 | 统一调度 `smoke`、`verify`、`original-sample`、`build`、`train` 和 `plots` |
| 06 | `mvi_06_smoke_test.py` | 主流程 | 用合成数据进行两轮快速 API/环境检查 |
| 07 | `mvi_07_verify_inputs.py` | 主流程 | 审计 H5AD、ALLC、5-kb 坐标、样本元数据和 IR/NR 数量 |
| 08 | `mvi_08_plot_original_embedding.py` | 主流程 | 在原始 ALLCools UMAP/t-SNE 上绘制 sample 与 condition |
| 09 | `mvi_09_build_input.py` | 主流程 | 从 ALLC 文件重建 MethylVI 所需的整数 `mc/cov` 层，生成 H5MU |
| 10 | `mvi_10_train_model.py` | 主流程 | 训练 MethylVI，计算 latent、邻居图、UMAP 和 Leiden |
| 11 | `mvi_11_plot_celltype_sample.py` | 主流程 | 绘制训练后 UMAP 的 cell type 与 sample |
| 12 | `mvi_12_plot_condition.py` | 主流程 | 绘制训练后 UMAP 的 IR/NR condition |
| 13 | `mvi_13_submit_cpu.sbatch` | 集群入口 | 在 Slurm 中提交并运行 `mvi_05_run_pipeline.sh all` |
| 14 | `mvi_14_utils_pipeline.py` | 支撑模块 | ID 标准化、样本元数据合并、区域解析、ALLC 聚合和 JSON 输出公共函数 |
| 15 | `mvi_15_utils_plot.py` | 支撑模块 | 分类嵌入绘图公共函数 |
| 16 | `mvi_16_requirements.lock.txt` | 环境记录 | 核心软件版本记录 |
| 17 | `mvi_17_test_utils_pipeline.py` | 测试 | 不接触真实数据的轻量单元测试 |

## 输入配置

编辑 `mvi_00_config.sh`，至少设置：

```bash
export MVI_DATA_ROOT=/path/to/your/project
export MVI_METHSCAN_UPSTREAM="/Users/luozhixiong/Library/Mobile Documents/com~apple~CloudDocs/Documents/PHD/脚本/Methscan/01_Upstream"
export MVI_H5AD=$MVI_DATA_ROOT/mcg_5kb.clustered.h5ad
export MVI_ALLC_DIR=$MVI_DATA_ROOT/input_allc
export MVI_SAMPLE_METADATA=/path/to/20260810/scripts/mvi_01_sample_metadata.tsv
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

实际文件可用 TSV 或 CSV。注释文件 `MVI_ANNOTATION` 可留空；如使用，至少需要 `cell_id`，并可增加 `sample_id`、`condition`、`cell_type` 及真实技术批次列。若同时提供 sample metadata，`condition` 以 sample metadata 为准。

## 服务器验证记录

截至 2026-08-10，服务器上的脚本工作目录已经确认。服务器上的 MethSCAn 数据根目录也已确认，但具体路径不在公开 README 中记录。

已确认 10 个样本目录全部存在，`cov_dedup_probability` 中约有 5.85 万个细胞文件，IR 和 NR 两组均包含 5 个样本。

当前服务器上尚未发现 H5AD、H5MU 或 ALLC 文件，只有按样本分目录保存的 `.cov.gz`。因此不能直接运行 `verify`、`build` 或 `train`。`mvi_03_allcools_prepare_5kb_counts.sbatch` 已更新为同时支持平铺 cov 目录和嵌套的 `*_Met/cov_dedup_probability` 目录，并会自动为嵌套输入添加样本前缀、检查 cell ID 冲突，再生成 ALLC、MCDS 和 5-kb H5AD。

服务器运行结果、脱敏后的核查结论和后续修正应持续补充到本节；具体路径仅在服务器上的配置文件中维护。

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

MethylVI 应优先使用 `cov_dedup_probability`，不要直接使用原始 `cov`。`mvi_01_sample_metadata.tsv` 已按 `IR01–IR05` 和 `NR01–NR05` 建好；`mvi_00_config.sh` 中的正则可以从合并 cell ID（如 `25110891_IR01_Met__barcode`）提取对应的 `IR01`。当前检查发现本地 `01_Upstream` 目录保存的是脚本、报告和说明文件，未发现 `.cov/.cov.gz`、H5AD 或 H5MU 数据文件；实际数据根目录由 MethSCAn 脚本中的 `BASE_DIR` 指定，若运行环境不是 `/share/LCZX_Data/data/allcools`，需要在运行前覆盖相应路径。

如果使用 MethSCAn 的 `.cov.gz` 作为 ALLCools 上游输入，`mvi_03_allcools_prepare_5kb_counts.sbatch` 已支持 `.cov.gz`，例如：

```bash
sbatch mvi_03_allcools_prepare_5kb_counts.sbatch \
  /path/to/merged_or_staged_cov_dedup_probability \
  /path/to/allcools_5kb_ir_nr \
  /path/to/hg38.chrom.sizes
```

如果 cov 按 10 个样本分别放在 `*_Met/cov_dedup_probability` 下，也可以直接把包含这些样本目录的数据根目录作为第一个参数；脚本会在输出目录下建立 `staged_cov/` 软链接，并把 cell ID 规范为 `IR01__barcode` 或 `NR01__barcode`。

可选的 ALLCools 上游步骤使用参数或环境变量提供路径，不含任何旧项目固定路径：

```bash
sbatch mvi_03_allcools_prepare_5kb_counts.sbatch /path/to/cov /path/to/allcools_output /path/to/chrom.sizes
```

等价环境变量为 `MVI_COV_DIR`、`MVI_ALLCOOLS_OUTPUT`、`MVI_CHROM_SIZES`；软件命令可用 `ALLCOOLS_EXE` 和 `PYTHON_BIN` 覆盖。

## 数据处理流程与脚本对应关系

整个流程分为“可选的 ALLCools 上游准备”和“MethylVI 主流程”两部分。每一步的脚本、输入和输出如下。

| 步骤 | 处理内容 | 使用脚本 | 主要输入 | 主要输出 |
|---|---|---|---|---|
| 0. 环境与参数 | 设置数据目录、H5AD、ALLC 目录、样本元数据和模型参数 | `mvi_00_config.sh` | 配置文件 | 环境变量 |
| 1. ALLCools 上游（可选） | 将每个 `.cov`/`.cov.gz` 转为压缩 ALLC；生成 5-kb MCDS；按 mCG hypo-score 做聚类和 UMAP | `mvi_03_allcools_prepare_5kb_counts.sbatch`、`mvi_04_allcools_cluster_5kb.py` | MethSCAn `cov_dedup_probability`、染色体长度文件 | `input_allc/`、MCDS、`mcg_5kb.clustered.h5ad` |
| 2. 依赖与 API 检查 | 用小型合成数据训练两轮，确认 Python、MethylVI、scvi-tools 和 PyTorch 环境可用 | `mvi_06_smoke_test.py` | 无真实数据 | 终端通过/失败状态 |
| 3. 输入审计 | 检查 H5AD 细胞和 5-kb bins、ALLC 文件匹配、样本 ID、IR/NR=5/5、坐标和元数据 | `mvi_07_verify_inputs.py`；公共函数来自 `mvi_14_utils_pipeline.py` | H5AD、ALLC、`mvi_01_sample_metadata.tsv`、可选注释表 | `${MVI_AUDIT}` |
| 4. 原始嵌入对照 | 在 ALLCools 已有 UMAP/t-SNE 上绘制 sample 和 IR/NR，作为 MethylVI 前的对照 | `mvi_08_plot_original_embedding.py`、`mvi_15_utils_plot.py` | H5AD、样本元数据 | `allcools_original_embedding_sample.pdf`、`allcools_original_embedding_condition.pdf` |
| 5. 计数重建 | 从原始 ALLC 在保留的 5-kb bins 内重新聚合 mCG 甲基化数和覆盖数；不使用 H5AD 的聚类分数作为 MethylVI 计数 | `mvi_09_build_input.py`；公共函数来自 `mvi_14_utils_pipeline.py` | ALLC、H5AD 中的 bins | `${MVI_INPUT}`、`count_rows/`、`build_summary.json` |
| 6. MethylVI 训练 | 读取 H5MU 的 `mc/cov`，默认按 `sample_id` 做 batch correction；计算 latent、邻居图、UMAP 和 Leiden | `mvi_10_train_model.py` | `${MVI_INPUT}`、`MVI_SAMPLE_KEY`、`MVI_CONDITION_KEY`、`MVI_BATCH_KEY` | 模型、latent、embedding、训练历史、`sample_by_condition.csv` |
| 7. 结果可视化 | 按 cell type、sample 和 condition 绘制训练后 UMAP | `mvi_11_plot_celltype_sample.py`、`mvi_12_plot_condition.py`、`mvi_15_utils_plot.py` | `methylvi_embedding.h5ad` 或 `cell_annotations_umap.tsv.gz` | `methylvi_umap_cell_type.pdf`、`methylvi_umap_sample_id.pdf`、`methylvi_umap_condition.pdf` |

其中第 1 步只有在当前项目还没有合格的 ALLCools 5-kb 结果时才执行；如果已有 `mcg_5kb.clustered.h5ad` 和对应的逐细胞 ALLC 文件，可直接从第 2 步开始。第 3 步必须在第 5、6 步之前通过。第 4 步是校正前对照，不改变输入数据；第 5 步才是真正生成 MethylVI 的 `mc/cov` 计数矩阵。

## 运行顺序

```bash
cd /path/to/20260810/scripts
bash mvi_05_run_pipeline.sh verify
bash mvi_05_run_pipeline.sh smoke
bash mvi_05_run_pipeline.sh original-sample
bash mvi_05_run_pipeline.sh build
bash mvi_05_run_pipeline.sh train
bash mvi_05_run_pipeline.sh plots
```

也可以运行 `bash mvi_05_run_pipeline.sh all`，或使用 `sbatch mvi_13_submit_cpu.sbatch` 提交完整流程。首次运行建议先单独执行 `verify` 和 `smoke`。

输入审计会拒绝：样本数不是 10、IR/NR 不是 5/5、10 个元数据样本没有全部出现在选中细胞中、无法从 cell ID 或注释得到 sample、condition 不是 IR/NR、缺少 ALLC 文件、5-kb 坐标无法解析，或坐标与 ALLC 不一致。

## 主要输出

- `${MVI_INPUT}`：包含 `mCG.layers['mc']` 和 `mCG.layers['cov']` 的 H5MU；
- `${MVI_RESULTS}/model/`：训练后的 MethylVI 模型；
- `${MVI_RESULTS}/latent_representation.npy`：latent 表示；
- `${MVI_RESULTS}/methylvi_embedding.h5ad`：邻居图、UMAP 和 Leiden；
- `${MVI_RESULTS}/cell_annotations_umap.tsv.gz`：UMAP、sample、condition、cell type 和 latent；
- `${MVI_RESULTS}/sample_by_condition.csv`：样本与 IR/NR 的细胞数交叉表；
- `${MVI_RESULTS}/methylvi_umap_condition.pdf`：IR/NR UMAP；
- `${MVI_RESULTS}/run_summary.json`、`${MVI_RESULTS}/training_history.csv`：参数、软件版本和训练记录；
- `${MVI_AUDIT}`：训练前输入审计报告。

## 断点与安全

`mvi_09_build_input.py` 按细胞写入 `count_rows/*.npz` 检查点，并用 manifest 防止不同输入被错误复用。已有形状和图层正确的 H5MU 会被验证后复用；需要强制重组时才使用 `--force-assemble`。不要让两个 build/train 作业同时写入相同的 `MVI_ROOT`。

## 参考资料

- [MethylVI paper](https://www.nature.com/articles/s42256-026-01225-9)
- [MethylVI model documentation](https://docs.scvi-tools.org/en/latest/user_guide/models/methylvi.html)
- [MethylVI integration tutorial](https://docs.scvi-tools.org/en/latest/tutorials/notebooks/scbs/MethylVI_batch.html)
- [MethylVI reproducibility repository](https://github.com/suinleelab/methylVI-reproducibility)
- [ALLCools documentation](https://lhqing.github.io/ALLCools/intro.html)
