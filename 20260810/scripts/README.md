# MethylVI 5-kb分析流程：IR与NR项目

本目录保存当前10个样本（5个IR、5个NR）的ALLCools与MethylVI分析脚本。流程从MethSCAn筛选后的细胞出发，重建5-kb mCG整数计数，使用MethylVI按样本进行batch correction，再生成latent、UMAP、Leiden及分组图。

本文档整理至2026-08-11。原有资料不会被覆盖；大任务统一通过`dsub`提交，服务器不使用Slurm。

## 当前状态

| 项目 | 状态 |
|---|---|
| 10个样本原始ALLC及TBI | 已验收，共58,534个细胞 |
| MethSCAn 300k细胞QC | 已验收，共6,199个细胞 |
| ALLCools 5-kb MCDS与聚类 | 已完成并验收 |
| SCANPY cell type注释合并 | 已完成，匹配5,765个细胞 |
| MethylVI H5MU输入 | 已完成并验收，约1.03 GiB |
| MethylVI环境与Leiden | smoke test和实际运行验收均通过 |
| 正式CPU训练 | 任务`164095`正在运行 |

截至2026-08-11 14:41，正式任务已进入`Epoch 10/500`，平均约52.59秒/epoch，`train_loss_epoch`从约8.44万下降至约8.21万，没有`Traceback`或`ERROR`。500轮是上限，训练启用early stopping。

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
| SCANPY注释表 | `/share/home/rzli/SCANPY/20260714/result/annotation/02_cell_annotation_all_cells.csv` |
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
| 04 | `mvi_04_run_pipeline.sh` | 主入口：`smoke`、`verify`、`original-sample`、`build`、`train`、`plots`、`all` |
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
| 7. MethylVI训练 | 04、09 | H5MU、`sample_id` batch | 模型、latent、UMAP、Leiden、训练记录 | 运行中 |
| 8. 校正后绘图 | 04、10、11、14 | MethylVI embedding | 3张校正后PDF | 待训练完成 |

主流程命令为：

```bash
cd /share/home/rzli/MethylVI/20260810/scripts

bash mvi_04_run_pipeline.sh smoke
bash mvi_04_run_pipeline.sh verify
bash mvi_04_run_pipeline.sh original-sample
bash mvi_04_run_pipeline.sh build
bash mvi_04_run_pipeline.sh train
bash mvi_04_run_pipeline.sh plots
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

`MVI_SAMPLE_ID_REGEX`可识别`25110891_IR01_Met__barcode`、`IR01__barcode`和`IR01_barcode`等形式，并统一生成`IR01__barcode`形式的cell ID。

## 已验收的关键结果

### ALLCools 5-kb结果

- 完整任务：`163901`，32 CPU、190 GiB，运行约5小时，峰值内存约20.3 GiB，成功退出。
- MCDS约3.8 GiB。
- 初始矩阵为6,199×617,665；过滤后H5AD为6,199×231,648。
- H5AD包含`X_pca`、`X_tsne`和`X_umap`，cell ID唯一。
- 最终`L1`为4簇：c0 3,600、c1 1,398、c2 665、c3 536。
- 十折交叉验证准确率为0.964。

### SCANPY注释

- 注释表含56,746个细胞。
- 与当前6,199个细胞成功匹配5,765个，匹配率为0.92999；434个记为`Unannotated`。
- 5,748个匹配细胞为`Keep`，17个为`Exclude`。
- `cell_type_integrated`映射为内部`cell_type`；样本和分组最终以`mvi_01_sample_metadata.tsv`为准。
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

完整软件版本见`mvi_15_requirements.lock.txt`。

## 输出目录

图像与数据结果分开保存：

```text
/share/home/rzli/MethylVI/20260810/result/
├── 01_before_methylvi/       # ALLCools及MethylVI校正前图
│   └── legacy_results_ir_nr_20260811/  # 保留的另一批旧版PDF
└── 02_after_methylvi/        # MethylVI校正后图

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
    └── run_summary.json
```

训练完成后，10和11脚本把下列图保存到`02_after_methylvi/`：

- `methylvi_umap_cell_type.pdf`
- `methylvi_umap_sample_id.pdf`
- `methylvi_umap_condition.pdf`

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
