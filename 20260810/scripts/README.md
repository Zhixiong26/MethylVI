# MethylVI 5-kb分析流程：IR与NR项目

本目录保存当前10个样本（5个IR、5个NR）的ALLCools与MethylVI分析脚本。流程从MethSCAn筛选后的细胞出发，重建5-kb mCG整数计数，使用MethylVI按样本进行batch correction，再生成latent、UMAP、Leiden及分组图。

本文档整理至2026-08-11。原有资料不会被覆盖；大任务统一通过`dsub`提交，服务器不使用Slurm。

## 当前状态

| 项目 | 状态 |
|---|---|
| 10个样本原始ALLC及TBI | 已验收，共58,534个细胞 |
| MethSCAn 300k细胞QC | 已验收，共6,199个细胞 |
| ALLCools 5-kb MCDS与聚类 | 已完成并验收 |
| SCANPY cell type注释 | 已切换到20260810新注释；任务`164099`审计通过，匹配5,765/6,199 |
| MethylVI H5MU输入 | 已完成并验收，约1.03 GiB |
| MethylVI环境与Leiden | smoke test和实际运行验收均通过 |
| 正式CPU训练 | 任务`164095`已成功完成 |
| 新注释图片重绘 | 校正前任务`164100`、校正后任务`164101`均成功 |
| supervised UMAP | 首次任务`164097`失败后已修复；重试任务`164102`成功生成4组UMAP和12张PDF |

正式任务`164095`已退出活动队列，调度器报告`EXIT_CODE=0`和`Job execution succeeded`。`run_summary.json`与训练历史确认实际训练89轮（epoch编号0–88），未达到500轮上限，由early stopping正常终止；随后成功生成模型、latent、邻居图、UMAP、Leiden及完整结果表。新注释审计、校正前后重绘和四组supervised UMAP也均已完成。

2026-08-11 16:12:49，校正后绘图任务`164096`通过`dsub`提交，申请4 CPU和16,384 MiB，并在`node-11`成功完成。调度器报告`EXIT_CODE=0`和`Job execution succeeded`，运行9秒，峰值内存82 MiB，没有`Traceback`、`ERROR`或`Exception`。`02_after_methylvi`已生成cell type、condition和sample三张单页PDF，文件大小分别为121 KiB、108 KiB和114 KiB，均通过PDF文件类型检查。至此本项目从输入整理到最终绘图的完整流程已全部完成。

## 项目设计

- 样本为`IR01–IR05`和`NR01–NR05`；`sample_id`表示样本，`condition`只允许`IR`或`NR`。
- 最终使用6,199个通过MethSCAn 300k QC的细胞。这里的300k表示每个细胞至少覆盖300,000个CpG位点，不是筛选300,000个特征。
- MethSCAn筛选参数为`min_sites=300000`、`max_sites=10000000`、`min_meth=55`、`max_meth=none`。
- MethylVI输入特征为ALLCools保留的231,648个5-kb mCG bins；H5AD中的聚类分数不作为MethylVI计数。
- `MVI_BATCH_KEY=sample_id`，即MethylVI按10个样本校正batch。
- 本项目中样本与IR/NR分组绑定，sample batch和condition不是独立因素。校正可能同时削弱真实IR/NR差异，因此必须同时保留校正前图，并结合细胞类型和样本分布解释校正后结果。如以后获得独立技术批次，应改为`MVI_BATCH_KEY=technical_batch`。

## 服务器路径

| 内容 | 绝对路径 |
|---|---|
| 脚本工作目录 | `/share/home/rzli/MethylVI/20260810/scripts` |
| 项目图像目录 | `/share/home/rzli/MethylVI/20260810/result` |
| MethSCAn/ALLCools数据根目录 | `/share/LCZX_Data/data/allcools` |
| 原始样本目录 | `/share/LCZX_Data/data/allcools/25110891_<sample>_Met` |
| 原始ALLC目录 | `/share/LCZX_Data/data/allcools/25110891_<sample>_Met/allcools` |
| MethSCAn 300k白名单 | `<sample>/qc_minmeth55_maxmethnone_maxsites10000000_covdedupprob/filtered_data_single_300k/column_header.txt` |
| 300k ALLCools输出 | `/share/LCZX_Data/data/allcools/methylvi_5kb_300k` |
| 5-kb聚类H5AD | `/share/LCZX_Data/data/allcools/methylvi_5kb_300k/mcg_5kb.clustered.h5ad` |
| MethylVI输出根目录 | `/share/LCZX_Data/data/allcools/methylVI_results_300k` |
| MethylVI输入H5MU | `/share/LCZX_Data/data/allcools/methylVI_results_300k/methylvi_5kbin_input.h5mu` |
| 正式模型结果 | `/share/LCZX_Data/data/allcools/methylVI_results_300k/results_ir_nr` |
| SCANPY注释表 | `/share/home/rzli/SCANPY/20260810/Result0810/annotation/02_cell_annotation_all_cells.csv` |
| hg38染色体长度 | `/share/home/rzli/MethylVI/20260810/scripts/hg38.canonical.chrom.sizes` |

本地MethSCAn上游资料位于：

```text
/Users/luozhixiong/Library/Mobile Documents/com~apple~CloudDocs/Documents/PHD/脚本/Methscan/01_Upstream
```

## 样本与细胞数

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

QC后IR共3,219个细胞，NR共2,980个细胞。10个样本的原始ALLC数量与TBI索引数量逐一相等。

## 文件与运行顺序

| 顺序 | 文件 | 作用 |
|---:|---|---|
| 00 | `mvi_00_config.sh` | 集中定义输入、输出、样本规则、batch correction、环境、资源和模型参数 |
| 01 | `mvi_01_sample_metadata.tsv` | 10个真实样本的`sample_id`与`condition`对应表 |
| 02 | `mvi_02_allcools_prepare_5kb_counts.sh` | 按300k白名单整理原始ALLC，生成MCDS并调用03；`.cov.gz`仅作备用输入 |
| 03 | `mvi_03_allcools_cluster_5kb.py` | ALLCools 5-kb过滤、LSI、聚类、t-SNE和UMAP |
| 04 | `mvi_04_run_pipeline.sh` | 主入口：`smoke`、`verify`、`original-sample`、`build`、`train`、`plots`、`supervised`、`all` |
| 05 | `mvi_05_smoke_test.py` | 用合成数据快速检查MethylVI API和训练环境 |
| 06 | `mvi_06_verify_inputs.py` | 审计H5AD、ALLC、坐标、注释、样本及IR/NR分组 |
| 07 | `mvi_07_plot_original_embedding.py` | 绘制MethylVI校正前的cell type、sample和condition图 |
| 08 | `mvi_08_build_input.py` | 从ALLC重建整数`mc/cov`层并生成H5MU，支持逐细胞断点复用 |
| 09 | `mvi_09_train_model.py` | 训练MethylVI并生成latent、邻居图、UMAP和Leiden |
| 10 | `mvi_10_plot_celltype_sample.py` | 绘制校正后的cell type和sample UMAP |
| 11 | `mvi_11_plot_condition.py` | 绘制校正后的IR/NR UMAP |
| 12 | `mvi_12_run_pipeline.sh` | `dsub`完整流程执行入口 |
| 13 | `mvi_13_utils_pipeline.py` | ID标准化、注释合并、区域解析、ALLC聚合和JSON输出公共函数 |
| 14 | `mvi_14_utils_plot.py` | 嵌入图绘制公共函数 |
| 15 | `mvi_15_requirements.lock.txt` | 服务器实测软件版本 |
| 16 | `mvi_16_test_utils_pipeline.py` | 不读取真实数据的轻量单元测试 |
| 17 | `mvi_17_plot_supervised_umap.py` | 在已训练的MethylVI latent上生成`target_weight=0.2/0.5/0.7/0.9`的cell type标签引导UMAP |
| 文档 | `MethylVI_vs_yuanpei_workflow_comparison.md` | 比较当前流程与`yuanpei/reproducible_methylVI_pipeline`的MethylVI核心、工程差异、证据边界和结果可比性 |
| 参考 | `hg38.canonical.chrom.sizes` | hg38 canonical 24条染色体长度 |

## 数据处理流程

| 步骤 | 脚本 | 输入 | 输出 | 当前状态 |
|---:|---|---|---|---|
| 1. 300k白名单与ALLC整理 | 02 | 58,534个原始ALLC、10个QC header | 6,199个ALLC/TBI软链接、manifest | 已完成 |
| 2. MCDS与5-kb聚类 | 02、03 | 6,199个ALLC、hg38长度 | MCDS、聚类H5AD、t-SNE、UMAP | 已完成 |
| 3. 环境smoke test | 04、05 | 合成`mc/cov` | 2轮CPU训练和latent检查 | 已通过 |
| 4. 真实输入审计 | 04、06、13 | H5AD、ALLC、样本表、SCANPY注释 | `mvi_06_input_audit.json` | 已通过 |
| 5. 校正前绘图 | 04、07、14 | ALLCools embedding和注释 | 3张校正前PDF | 已完成 |
| 6. H5MU构建 | 04、08、13 | 6,199个ALLC、231,648个bins | H5MU、count checkpoints、build summary | 已完成 |
| 7. MethylVI训练 | 04、09 | H5MU、`sample_id` batch | 模型、latent、UMAP、Leiden、训练记录 | 已完成 |
| 8. 校正后绘图 | 04、10、11、14 | MethylVI embedding | 3张校正后PDF | 已完成 |
| 9. supervised UMAP | 04、17 | 已训练的`X_methylVI`和cell type注释 | 4组坐标、12张分类PDF、H5AD和JSON摘要 | 任务`164102`已完成 |

主流程命令为：

```bash
cd /share/home/rzli/MethylVI/20260810/scripts

bash mvi_04_run_pipeline.sh smoke
bash mvi_04_run_pipeline.sh verify
bash mvi_04_run_pipeline.sh original-sample
bash mvi_04_run_pipeline.sh build
bash mvi_04_run_pipeline.sh train
bash mvi_04_run_pipeline.sh plots
bash mvi_04_run_pipeline.sh supervised
```

大任务必须通过`dsub`执行，不要在登录节点直接运行Python分析。使用`djob`查看任务状态；日志统一保存在`scheduler_logs/`。当前正式训练日志为：

```text
scheduler_logs/methylvi_cpu_formal.164095.out
scheduler_logs/methylvi_cpu_formal.164095.err
```

## 核心配置

当前`mvi_00_config.sh`的正式设置为：

| 参数 | 当前值 | 含义 |
|---|---|---|
| `MVI_EXPECTED_SAMPLES` | `10` | 样本数检查 |
| `MVI_EXPECTED_IR` | `5` | IR样本数检查 |
| `MVI_EXPECTED_NR` | `5` | NR样本数检查 |
| `MVI_EXPECTED_CELLS` | `6199` | 02上游白名单和ALLC整理期望数；06审计不强制固定H5AD细胞数 |
| `MVI_BATCH_KEY` | `sample_id` | MethylVI batch correction字段 |
| `MVI_BIN_SIZE` | `5000` | 区域宽度 |
| `MVI_MC_CONTEXT` | `CGN` | 甲基化context |
| `MVI_BATCH_SIZE` | `32` | 训练batch size |
| `MVI_MAX_EPOCHS` | `500` | 最大训练轮数；启用early stopping |
| `MVI_N_LATENT` | `20` | latent维度 |
| `MVI_N_HIDDEN` | `128` | hidden维度 |
| `MVI_N_LAYERS` | `1` | hidden层数 |
| `MVI_NEIGHBORS` | `15` | latent邻居图参数 |
| `MVI_LEIDEN_RESOLUTION` | `1.0` | Leiden resolution |
| `MVI_THREADS` | `32` | 任务级CPU线程数 |
| `MVI_ACCELERATOR` | `cpu`（正式任务显式覆盖） | 服务器无GPU |
| `MVI_SUPERVISED_TARGET_KEY` | `cell_type` | supervised UMAP的标签目标 |
| `MVI_SUPERVISED_TARGET_WEIGHTS` | `0.2 0.5 0.7 0.9` | 依次计算的标签引导强度 |
| `MVI_SUPERVISED_MIN_DIST` | `0.5` | supervised UMAP的`min_dist` |

`MVI_SAMPLE_ID_REGEX`可识别`25110891_IR01_Met__barcode`、`IR01__barcode`和`IR01_barcode`等形式，并统一生成`IR01__barcode`形式的cell ID。

### MethSCAn与ALLC整理参数

| 参数 | 实际值 |
|---|---|
| 使用筛选后细胞 | `MVI_USE_FILTERED_CELLS=1` |
| 最少覆盖CpG位点 | `300000` |
| 最多覆盖CpG位点 | `10000000` |
| 最低甲基化筛选值 | `55` |
| 最高甲基化筛选值 | `none` |
| 预期细胞数 | `6199` |
| 预期样本结构 | 10个样本，5 IR + 5 NR |
| ALLC来源优先级 | 原始`*_allc.gz`优先，`cov_dedup_probability/*.cov.gz`仅作备用 |
| ALLC甲基化context | `CGN` |
| 染色体范围 | hg38 canonical 24条染色体 |

### ALLCools MCDS与聚类参数

| 环节 | 参数 | 实际值 |
|---|---|---|
| generate-dataset | observation维度 | `cell` |
| generate-dataset | 区域 | `chrom5k 5000` |
| generate-dataset | count quantifier | `chrom5k count CGN` |
| generate-dataset | hypo-score quantifier | `chrom5k hypo-score CGN cutoff=0.9` |
| generate-dataset | chunk size | `10` |
| score matrix | `mc_type` / `quant_type` | `CGN` / `hypo-score` |
| binarization | cutoff | `0.95` |
| region filtering | 参数 | 使用ALLCools `filter_regions()`默认值 |
| LSI | algorithm | `arpack` |
| LSI | random seed | `0` |
| significant PC test | `p_cutoff` | `0.1` |
| 邻居图 | `n_neighbors` | `25` |
| 初始Leiden | resolution / seed | `1.0` / `0` |
| t-SNE | metric | `euclidean` |
| t-SNE | perplexity | `30` |
| t-SNE | exaggeration | `-1`（由ALLCools自动策略处理） |
| UMAP | random seed | `0`，其余沿用Scanpy默认值 |
| ConsensusClustering | `n_neighbors` | `25` |
| ConsensusClustering | metric | `euclidean` |
| ConsensusClustering | `min_cluster_size` | `10` |
| ConsensusClustering | `leiden_repeats` | `500` |
| ConsensusClustering | `leiden_resolution` | `0.5` |
| ConsensusClustering | `consensus_rate` | `0.5` |
| ConsensusClustering | `train_frac` / `train_max_n` | `0.5` / `500` |
| ConsensusClustering | `max_iter` / seed | `20` / `0` |

实际运行中显著LSI成分为8个，这是根据`p_cutoff=0.1`计算得到的结果，不是预先固定的参数。

### H5MU计数构建参数

| 参数 | 实际值 |
|---|---|
| 细胞数 | `6199` |
| 特征数 | `231648`个保留5-kb bins |
| bin size | `5000` |
| context | `CGN` |
| 并行worker | `32` |
| 输出dtype策略 | `auto` |
| 实际dtype | `uint16` |
| 强制重组 | `--force-assemble`未使用 |
| 断点复用 | 6,199个`count_rows/*.npz`全部复用 |

### MethylVI正式模型参数

| 环节 | 参数 | 实际值 |
|---|---|---|
| setup_mudata | methylated count layer | `mc` |
| setup_mudata | coverage layer | `cov` |
| setup_mudata | modality/context | `mCG` |
| batch correction | batch key | `sample_id`（10个样本） |
| 网络 | likelihood | `betabinomial` |
| 网络 | dispersion | `region` |
| 网络 | latent维度 | `20` |
| 网络 | hidden维度 | `128` |
| 网络 | hidden层数 | `1` |
| 训练 | 最大epochs | `500` |
| 训练 | early stopping | `True` |
| 训练 | batch size | `32` |
| 训练 | accelerator | `cpu` |
| 训练 | devices | `1` |
| 训练 | PyTorch/SCVI线程 | `32` |
| 训练 | random seed | `0` |
| latent邻居图 | `n_neighbors` | `15` |
| latent邻居图 | representation | `X_methylVI` |
| latent UMAP | random seed | `0`，其余沿用Scanpy默认值 |
| latent Leiden | resolution / seed | `1.0` / `0` |

学习率、优化器、训练/验证划分比例、early-stopping patience等没有在本项目脚本中显式覆盖，使用已锁定的scvi-tools 1.3.3内部默认值。若后续需要严格记录这些内部默认值，应从本次保存的模型和训练配置中导出，不应凭其他版本文档推断。

### 计算资源与测试参数

| 任务 | 参数或资源 | 实际值 |
|---|---|---|
| 300k stage-only | dsub | 4 CPU、16 GiB，`MVI_STAGE_ONLY=1` |
| 完整ALLCools | dsub | 32 CPU、194,560 MiB |
| 首次H5MU build | dsub | 32 CPU、194,560 MiB |
| H5MU重试 | dsub | 32 CPU、65,536 MiB |
| CPU 1轮基准 | dsub | 32 CPU、65,536 MiB，1 epoch，batch size 32 |
| 正式MethylVI | dsub | 32 CPU、65,536 MiB，最多500 epochs |
| smoke test | 合成矩阵 | 32细胞×64特征，coverage随机范围0–7 |
| smoke test | 模型 | latent 4、hidden 16、1层、2 epochs、batch size 16、CPU |

### 绘图参数

| 参数 | 实际值 |
|---|---|
| 散点大小 | `s=4` |
| 透明度 | `alpha=0.75` |
| 边线宽度 | `0` |
| 画布 | 类别≤12时`8×7`，否则`11×7`英寸 |
| 调色板 | 类别≤20时`tab20`，否则`gist_ncar` |
| 图例字体/点缩放 | `fontsize=7`、`markerscale=2.5` |
| 绘制顺序随机种子 | `0` |
| ALLCools PNG分辨率 | `300 dpi` |

### supervised UMAP参数与含义

17脚本不重新训练MethylVI，而是在已保存的20维`X_methylVI`上重新计算UMAP。

| 参数 | 值 |
|---|---|
| 目标标签 | `cell_type` |
| target weights | `0.2`、`0.5`、`0.7`、`0.9` |
| neighbors | `15` |
| metric / target metric | `euclidean` / `categorical` |
| min_dist / spread | `0.5` / `1.0` |
| seed | `0` |

新注释表中未注释或`exclude_from_main_analysis=True`的细胞使用`-1`作为未标记目标：它们仍会出现在UMAP中，但不用cell type标签拉近。具体数量由17脚本在运行时根据新注释表计算，不沿用旧注释的434个未匹配统计。`target_weight`越大，cell type标签对二维布局的影响越强。这些图属于标签引导可视化，不能当作独立证据证明MethylVI无监督latent自然形成了同样的分类边界。

## 已验收的关键结果

### ALLCools 5-kb结果

- 完整任务：`163901`，32 CPU、190 GiB，运行约5小时，峰值内存约20.3 GiB，成功退出。
- MCDS约3.8 GiB。
- 初始矩阵为6,199×617,665；过滤后H5AD为6,199×231,648。
- H5AD包含`X_pca`、`X_tsne`和`X_umap`，cell ID唯一。
- 最终`L1`为4簇：c0 3,600、c1 1,398、c2 665、c3 536。
- 十折交叉验证准确率为0.964。

### SCANPY注释

- 当前默认路径为`/share/home/rzli/SCANPY/20260810/Result0810/annotation/02_cell_annotation_all_cells.csv`。
- 新注释审计任务`164099`成功，`EXIT_CODE=0`，运行10秒，峰值内存437 MiB。
- 新表共56,746行；与6,199个MethylVI细胞匹配5,765个，未匹配434个，匹配率0.9299887；排除17个，保留5,748个匹配细胞。
- 新表与旧`20260714`表的总行数、匹配数和排除数恰好相同，但cell type标签内容可以发生变化，因此仍已全部重绘相关图片。
- `cell_type_integrated`映射为内部`cell_type`；样本和分组最终以`mvi_01_sample_metadata.tsv`为准。
- 07、10和17脚本运行时都会重新读取当前`MVI_ANNOTATION`，因此更新cell type不需要重建H5MU或重训MethylVI。
- 8个公共函数单元测试全部通过。

### MethylVI H5MU

- 构建任务`164090`成功，复用了6,199个逐细胞计数检查点。
- H5MU为6,199×231,648，`mc/cov`均为`uint16`整数层。
- 最大`mc=7172`、最大`cov=8287`，所有抽查满足`0≤mc≤cov`。
- 预计稠密两层内存约5.35 GiB，最终H5MU约1.03 GiB。
- 最初写入失败是AnnData不能直接序列化`numpy.memmap`；08脚本已使用`np.asarray()`标准ndarray视图修复，小型写入和真实构建均通过。

### MethylVI环境与训练

- 服务器为linux-aarch64，无GPU，正式训练使用CPU。
- MethylVI来自`scvi.external.METHYLVI`，不是独立的`methyl_vi`包。
- CPU smoke test使用32×64合成数据训练2轮并成功获得32×4 latent。
- 真实数据1轮基准耗时56.44秒，峰值内存约9.36 GiB，平均使用约8.27个CPU核。
- `igraph 1.0.0`、`leidenalg 0.12.0`和`texttable 1.7.0`已安装；100细胞Scanpy Leiden运行验收通过。
- 正式任务`164095`申请32 CPU和64 GiB，使用batch size 32、最多500轮及early stopping。
- 正式任务最终`EXIT_CODE=0`；调度器总运行4,819秒（约80分19秒），其中模型报告训练运行4,745.02秒（约79分05秒）。内存峰值9,358 MiB、平均9,157 MiB，平均CPU利用率1,838%（约18.38核）。
- 训练历史共89条记录（epoch 0–88），说明early stopping在89轮后生效。`train_loss_epoch`由92,240.95降至76,516.16；最后三轮仍轻微下降，但`validation_loss`由81,154.08升至81,281.05，符合验证指标无持续改善后的停止行为。
- 正式输入为6,199细胞×231,648特征，IR 3,219、NR 2,980；batch key为`sample_id`，latent 20维，hidden 128×1层，batch size 32，likelihood为`betabinomial`，dispersion为`region`。
- 已生成461 MiB模型、485 KiB latent、3.0 MiB embedding H5AD、788 KiB细胞UMAP表、14 KiB训练历史和879字节运行摘要。
- 正式运行版本为scvi-tools 1.3.3、Torch 2.12.1、AnnData 0.11.4、MuData 0.3.9和Scanpy 1.11.5；版本锁文件已按正式任务的`run_summary.json`将MuData从旧记录0.3.2更正为0.3.9。

完整软件版本见`mvi_15_requirements.lock.txt`。

## 输出目录

图像与数据结果分开保存：

```text
/share/home/rzli/MethylVI/20260810/result/
├── 00_previous_annotation_20260714/  # 更新前校正前后图的安全归档
├── 01_before_methylvi/       # ALLCools及MethylVI校正前图
│   └── legacy_results_ir_nr_20260811/  # 保留的另一批旧版PDF
├── 02_after_methylvi/        # MethylVI校正后图
└── 03_supervised_umap/     # 第17步4组target_weight的标签引导UMAP

/share/LCZX_Data/data/allcools/methylVI_results_300k/
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

训练完成后，10和11脚本把下列图保存到`02_after_methylvi/`：

- `methylvi_umap_cell_type.pdf`
- `methylvi_umap_sample_id.pdf`
- `methylvi_umap_condition.pdf`

视觉检查显示：校正后10个样本仍在主要区域内充分混合，没有单一样本独占的大型孤立结构；B、CD14、CD4以及T/NK相关结构较校正前更清楚，未见明显cell type结构坍缩。IR/NR在校正前后均高度混合，因此UMAP只能支持“未见明显样本残留batch”的定性结论，不能证明condition信号被保留，也不能排除过校正。本阶段按当前需求暂不进行batch correction定量评估，17脚本改为target-weighted supervised UMAP。

17脚本已在`03_supervised_umap/target_weight_0p2`、`target_weight_0p5`、`target_weight_0p7`和`target_weight_0p9`四个子目录中，分别生成按cell type、sample和condition着色的PDF，共12张。

### supervised UMAP运行记录

- 首次任务`164097`于2026-08-11在node-11提交，申请2 CPU和16 GiB内存。
- 任务在17秒后代码退出，调度器`EXIT_CODE=10001`，峰值内存约500 MiB。
- 失败发生在UMAP计算前：`exclude_from_main_analysis`是Pandas `Categorical`列，旧代码`fillna(False)`尝试写入一个新category，引发`TypeError`。
- 17脚本已改为先转换到Pandas字符串类型，再填充缺失值和识别true/false。首次任务未产生PDF、坐标文件或摘要，不需要清理中间结果。
- 重试任务`164102`于2026-08-11在node-11成功完成，`EXIT_CODE=0`，运行69秒，峰值内存709 MiB，平均CPU利用率98%。
- 日志确认依次完成`target_weight=0.2`、`0.5`、`0.7`和`0.9`，最终生成4组坐标和12张PDF。
- 新注释下的普通校正前绘图任务`164100`成功，运行18秒，峰值内存884 MiB；校正后绘图任务`164101`成功，运行4秒，峰值内存95 MiB。

## 环境与已知问题

- 登录Shell默认没有初始化Conda。使用MethylVI前运行：

```bash
source /share/home/rzli/miniconda3/etc/profile.d/conda.sh
conda activate methylvi
```

- 登录节点运行轻量Python检查前，将以下变量设为1，避免OpenBLAS按整机CPU数创建线程并触发`RLIMIT_NPROC`：

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
```

- `python -m pip check`仍会提示`nvidia-cusparselt-cu13 0.8.1 is not supported on this platform`。服务器无GPU，而CPU smoke test和真实训练均已运行，因此当前不卸载该Torch依赖。
- MuData 0.4行为变化、Lightning pytree、DataLoader worker数、Numba及Scanpy colormap信息属于FutureWarning或性能提示，不等于任务失败。
- Scanpy提示Leiden未来默认后端将从`leidenalg`切换到`igraph`；当前版本和运行结果已按`leidenalg`验收，不在本次正式分析中途改变算法后端。
- 交互式Shell临时命令不要使用`exit`作为条件保护，否则会退出整个SSH会话；应只输出错误信息，或在子Shell/脚本中运行。
- 尖括号形式的`<JOB_ID>`只是文档占位符，不能原样输入Shell，否则会被解释为重定向。

## 断点与安全

- 08脚本逐细胞写入`count_rows/*.npz`，中断后可复用；manifest用于阻止不同输入误复用。
- H5MU存在时会先检查形状和图层；只有明确需要重组时才使用`--force-assemble`。
- 不要同时运行两个写入相同`MVI_ROOT`或`MVI_RESULTS`的build/train任务。
- 09脚本在读取大型H5MU及训练前检查`igraph`和`leidenalg`，避免训练结束后才因聚类依赖缺失失败。
- `model.save(..., overwrite=True)`会覆盖同名模型目录，因此正式提交前必须确认`MVI_RESULTS`为空或使用新的输出目录。
- 不删除内容不同的同名文件；旧图已无覆盖迁移到项目result下的legacy目录。

## 参考资料

- [MethylVI paper](https://www.nature.com/articles/s42256-026-01225-9)
- [MethylVI model documentation](https://docs.scvi-tools.org/en/latest/user_guide/models/methylvi.html)
- [MethylVI integration tutorial](https://docs.scvi-tools.org/en/latest/tutorials/notebooks/scbs/MethylVI_batch.html)
- [MethylVI reproducibility repository](https://github.com/suinleelab/methylVI-reproducibility)
- [ALLCools documentation](https://lhqing.github.io/ALLCools/intro.html)
