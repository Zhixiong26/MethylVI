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

`all`只执行`verify → build → train → plots → supervised`，不会自动执行耗时的`prepare`，也不会重复运行测试。`prepare`直接使用`MVI_ALLCOOLS_ENV`；其他阶段使用`MVI_CONDA_ENV`。

| 项目 | 当前状态 |
|---|---|
| 10个样本原始ALLC及TBI | 已验收，58,534个细胞 |
| MethSCAn 300k细胞QC | 已验收，6,199个细胞 |
| ALLCools 5-kb MCDS与聚类 | 已完成并验收 |
| 当前SCANPY注释 | 使用`20260810/Result0810`版本，匹配5,765/6,199 |
| MethylVI H5MU | 已完成，约1.03 GiB |
| 正式CPU训练 | 已完成，early stopping后共89条epoch记录 |
| 校正前后普通图 | 已用当前注释重画 |
| supervised UMAP | 0.2、0.5、0.7、0.9均已完成，共12张PDF |

当前设计固定为`IR01–IR05`和`NR01–NR05`。300k表示每个细胞至少覆盖300,000个CpG位点，并非保留300,000个特征。当前MethylVI使用231,648个ALLCools保留的5-kb bins；H5AD中的聚类分数不作为MethylVI计数。

## 2. 正式执行顺序

| 步骤 | 入口stage | 执行脚本 | 主要输入 | 主要输出 |
|---:|---|---|---|---|
| 1 | `prepare` | `02_prepare_allcools.sh`、`03_cluster_allcools.py` | 原始ALLC、300k白名单、hg38长度 | 6,199个ALLC/TBI链接、MCDS、聚类H5AD |
| 2 | `verify` | `04_verify_inputs.py` | H5AD、ALLC、样本表、当前注释 | `mvi_06_input_audit.json`（保留既有文件名） |
| 3 | `build` | `05_build_methylvi_input.py` | 6,199个ALLC、231,648个bins | H5MU、逐细胞检查点、build摘要 |
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

服务器更新时应一次性替换整个`scripts/`目录并移除旧编号脚本，避免新旧入口并存。已生成的数据、模型、坐标和图片不需要重建。

## 3. 文件作用

| 文件 | 作用 |
|---|---|
| `00_config.sh` | 输入输出路径、样本结构、QC、batch correction、环境、资源及模型参数 |
| `01_sample_metadata.tsv` | 10个实际样本的`sample_id`与`condition`对应表 |
| `02_prepare_allcools.sh` | 按300k白名单整理原始ALLC，生成MCDS并调用03；`.cov.gz`仅作备用 |
| `03_cluster_allcools.py` | ALLCools 5-kb过滤、LSI、聚类、t-SNE和UMAP |
| `04_verify_inputs.py` | 审计H5AD、ALLC、区域、样本结构和当前SCANPY注释 |
| `05_build_methylvi_input.py` | 从ALLC聚合整数`mc/cov`层，支持逐细胞检查点复用并写入H5MU |
| `06_train_methylvi.py` | 训练MethylVI并生成latent、邻居图、UMAP、Leiden及运行摘要 |
| `07_plot_embeddings.py` | `--stage before\|after\|all`；动态读取当前注释并生成6张普通图 |
| `08_plot_supervised_umap.py` | 在固定latent上生成四种`target_weight`的标签引导UMAP |
| `09_run_pipeline.sh` | `prepare/verify/build/train/plots/supervised/test/all`统一入口 |
| `mvi_utils.py` | ID标准化、元数据/注释合并、区域解析、ALLC聚合、JSON和绘图公共函数 |
| `requirements.lock.txt` | 正式环境的软件版本记录 |
| `hg38.canonical.chrom.sizes` | hg38 canonical 24条染色体长度 |
| `tests/test_mvi_utils.py` | 不读取真实数据的公共函数单元测试 |
| `tests/test_methylvi_smoke.py` | 32×64合成计数的两轮CPU MethylVI smoke test |
| `MethylVI_vs_yuanpei_workflow_comparison.md` | 当前流程与`yuanpei`流程的方法与工程比较 |

## 4. 参数配置

### 4.1 路径和样本

| 内容 | 当前值 |
|---|---|
| 脚本目录 | `/share/home/rzli/MethylVI/20260810/scripts` |
| 图像目录 | `/share/home/rzli/MethylVI/20260810/result` |
| 数据根目录 | `/share/LCZX_Data/data/allcools` |
| 原始样本目录 | `/share/LCZX_Data/data/allcools/25110891_<sample>_Met` |
| ALLCools 300k输出 | `/share/LCZX_Data/data/allcools/methylvi_5kb_300k` |
| 聚类H5AD | `/share/LCZX_Data/data/allcools/methylvi_5kb_300k/mcg_5kb.clustered.h5ad` |
| MethylVI输出根目录 | `/share/LCZX_Data/data/allcools/methylVI_results_300k` |
| MethylVI输入H5MU | `/share/LCZX_Data/data/allcools/methylVI_results_300k/methylvi_5kbin_input.h5mu` |
| 正式模型结果 | `/share/LCZX_Data/data/allcools/methylVI_results_300k/results_ir_nr` |
| 当前SCANPY注释 | `/share/home/rzli/SCANPY/20260810/Result0810/annotation/02_cell_annotation_all_cells.csv` |
| MethylVI环境 | `methylvi` |
| ALLCools环境 | `/share/home/rzli/miniconda3/envs/allcools` |

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
| binarize | cutoff | `0.95` |
| LSI | algorithm/seed | `arpack` / `0` |
| significant PC | p cutoff | `0.1`（本次得到8个成分） |
| 邻居/初始Leiden | 参数 | `25` / resolution `1.0` |
| t-SNE | 参数 | Euclidean、perplexity 30、exaggeration -1 |
| ConsensusClustering | 参数 | neighbors 25、resolution 0.5、repeats 500、min cluster 10、rate 0.5、seed 0 |
| H5MU | 形状 | 6,199×231,648 |
| H5MU | 层/dtype | `mc`和`cov`，实际`uint16` |
| H5MU build | workers | `32`，逐细胞检查点可复用 |

### 4.3 MethylVI正式参数

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

### 4.4 普通图与supervised UMAP

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

- ALLCools完整MCDS约3.8 GiB；初始矩阵6,199×617,665，过滤后H5AD为6,199×231,648。
- H5AD cell ID唯一，包含`X_pca`、`X_tsne`和`X_umap`；L1四簇为c0 3,600、c1 1,398、c2 665、c3 536，十折准确率0.964。
- 当前SCANPY表56,746行，匹配5,765个，未匹配434个，排除17个；cell type可更新，样本/分组以`01_sample_metadata.tsv`为准。
- H5MU为6,199×231,648，最大`mc=7172`、最大`cov=8287`，满足`0≤mc≤cov`；稠密内存展开约5.35 GiB，磁盘文件约1.03 GiB。
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

`164097`的修复是先将`exclude_from_main_analysis`转为Pandas字符串，再填充缺失并解析布尔值。`164094`的模型训练本身成功，错误发生在训练后的Leiden阶段。

## 6. 输出目录

```text
/share/home/rzli/MethylVI/20260810/result/
├── 00_previous_annotation_20260714/
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

图片只写入项目`result/`；数据、模型和坐标仍写入LCZX数据盘。更新cell type只需重新运行`verify`、`plots`和`supervised`，不需要重建H5MU或重训MethylVI。

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
