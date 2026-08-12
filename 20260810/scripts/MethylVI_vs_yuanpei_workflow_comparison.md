# `20260810/scripts`与`yuanpei` MethylVI流程的本质和技术实现比较

## 1. 修正后的结论

`yuanpei/reproducible_methylVI_pipeline`中确实存在完整的MethylVI训练流程。综合该目录的README、配置、入口脚本、输入审计和软件版本后，结论是：

> 两套流程的MethylVI核心统计路线本质相同，但输入数据、元数据字段、上游QC、工程实现、运行环境和附加分析并不完全相同。

更具体地说：

- 两者都从ALLCools聚类H5AD保留的5-kb bins出发。
- 两者都明确禁止把H5AD的`X`直接当作MethylVI计数。
- 两者都重新遍历逐细胞ALLC，构建整数`mc/cov`。
- 两者都把个体身份作为batch covariate：`yuanpei`使用`donor`，当前流程使用`sample_id`。
- 两者都训练20维MethylVI latent，并以15近邻、UMAP和Leiden resolution 1.0进行下游分析。
- 两者都设置最大500 epochs、batch size 32、early stopping和seed 0。

因此，当前流程并不是与`yuanpei`完全不同的方法，而是将同一套MethylVI核心路线适配到当前5 IR + 5 NR项目，并增加了更严格的上游审计、断点保护和服务器适配。当前阶段将可选的supervised UMAP设为0.2、0.5、0.7和0.9四个权重。

## 2. 本次比较实际读取到的文件

### 2.1 `yuanpei`的ALLCools上游

- `yuanpei/run_5kbin_from_cov.sbatch`
- `yuanpei/cluster_5kbin.py`

### 2.2 `yuanpei`的MethylVI复现入口

- `yuanpei/reproducible_methylVI_pipeline/README.md`
- `yuanpei/reproducible_methylVI_pipeline/config.sh`
- `yuanpei/reproducible_methylVI_pipeline/run_pipeline.sh`
- `yuanpei/reproducible_methylVI_pipeline/run_core_cu03.sbatch`
- `yuanpei/reproducible_methylVI_pipeline/00_verify_inputs.py`
- `yuanpei/reproducible_methylVI_pipeline/input_audit.json`
- `yuanpei/reproducible_methylVI_pipeline/requirements.lock.txt`

ZIP包中的文件与已解压目录一致。

### 2.3 当前目录

当前流程为`20260810/scripts/00_config.sh`至`09_run_pipeline.sh`，其中：

- 02–03：ALLCools上游；
- 04：输入审计；
- 05：重建MethylVI `mc/cov`；
- 06：训练MethylVI；
- 07：统一绘制校正前和校正后普通图；
- 08：`target_weight=0.2/0.5/0.7/0.9`的cell type标签引导UMAP；
- 09：统一执行入口。

## 3. 一个必须说明的证据边界

`yuanpei/reproducible_methylVI_pipeline/run_pipeline.sh`引用以下核心脚本：

```text
01_plot_existing_umap_donor.py
02_build_methylvi_input.py
03_train_methylvi_donor_batch.py
04_redraw_batch_corrected_celltype_donor.py
05_plot_batch_corrected_disease.py
results_donor_batch_corrected/method_1/run_supervised_umap.py
```

但这些文件不在当前`yuanpei`目录，也不在`reproducible_methylVI_pipeline.zip`中。因此：

- 可以根据README、参数、入口调用和已完成审计确认总体算法路线；
- 可以确认build、train、plots和supervised等阶段确实存在；
- 不能对`yuanpei`的02/03训练实现做逐行代码比较；
- 不能在缺少源码时断言其likelihood、dispersion、内部API调用和当前`06_train_methylvi.py`逐行一致。

本文会把“已由本地文件直接确认”和“根据入口/README确认”区分开，不把缺失源码部分写成已逐行验证。

## 4. 两条完整技术路线

### 4.1 `yuanpei`流程

```text
逐细胞cov
    ↓
转换为ALLC并建立tabix索引
    ↓
ALLCools 5-kb CGN hypo-score MCDS
    ↓
二值化、区域过滤、LSI、ConsensusClustering
    ↓
10,488个细胞 × 272,521个保留5-kb bins的H5AD
    ↓
从10,488个ALLC重新聚合整数mc/cov
    ↓
构建约14 GB H5MU
    ↓
以donor为batch covariate训练MethylVI
    ↓
20维donor-corrected latent
    ↓
15近邻、UMAP、Leiden resolution 1.0
    ↓
cell type、donor、disease绘图
    ↓
可选：使用标签引导的supervised UMAP
```

### 4.2 当前流程

```text
10个样本的58,534个原始ALLC
    ↓
验证MethSCAn 300k QC provenance
    ↓
筛选5 IR + 5 NR的6,199个细胞
    ↓
优先使用原始ALLC和TBI；cov只作备用
    ↓
ALLCools count + hypo-score 5-kb MCDS
    ↓
与yuanpei相同参数的二值化、LSI和ConsensusClustering
    ↓
6,199个细胞 × 231,648个保留5-kb bins的H5AD
    ↓
从6,199个原始ALLC重新聚合整数mc/cov
    ↓
合并sample、IR/NR和SCANPY cell type注释
    ↓
构建约1.03 GiB H5MU
    ↓
以sample_id为batch covariate训练MethylVI
    ↓
20维sample-corrected latent
    ↓
15近邻、UMAP、Leiden resolution 1.0
    ↓
cell type、sample、condition绘图
    ↓
可选：生成4个target_weight的cell type标签引导UMAP
```

## 5. 哪些步骤不同（修正版）

| 比较项 | `yuanpei` | `20260810/scripts` | 影响 |
|---|---|---|---|
| 完整目标 | ALLCools聚类 + MethylVI donor批次校正 | ALLCools聚类 + MethylVI sample批次校正 | 核心分析目标相同，batch字段语义略有不同 |
| 原始输入布局 | 上游从单个平铺`*.cov`目录转换ALLC；MethylVI阶段使用10,488个ALLC | 10个嵌套样本目录，优先直接使用原始带TBI的ALLC | 输入组织和上游来源不同，MethylVI都最终读取ALLC |
| 调度器 | Slurm `sbatch` | 集群`dsub` | 运行环境不同，不改变统计模型 |
| 线程与内存 | 默认50线程，记录为250 GB | 正式训练32线程、64 GiB | 主要影响速度和调度，不代表模型不同 |
| 细胞选择/QC | 对10,488个选中细胞、ALLC和donor赋值做输入审计；当前副本未记录MethSCAn 300k provenance | 核验MethSCAn 300k阈值、10份provenance和6,199细胞白名单 | 当前流程的QC来源和可追溯性更明确 |
| 样本/分组检查 | 验证donor/disease注释和匹配状态 | 强制10样本、5 IR、5 NR | 当前流程额外防止样本遗漏或混入 |
| cell ID | donor前缀并规范化`-/_` | `sample_id__barcode` | 当前命名可避免跨样本barcode冲突 |
| ALLC来源审计 | 有`input_audit.json`和ALLC映射检查 | 保存source manifest、selected-cell table、QC summary和校验信息 | 两者都审计，当前流程记录更细 |
| MCDS quantifier | 上游脚本只生成`hypo-score CGN` | 同时生成`count CGN`与`hypo-score CGN` | 当前MCDS保留的上游信息更完整 |
| ALLCools聚类 | LSI + ConsensusClustering | 同一算法和核心参数 | 该阶段本质相同 |
| 校正前UMAP | ALLCools聚类阶段根据`X_pca`计算；其01脚本只重新着色 | `03_cluster_allcools.py`根据`X_pca`计算；`07_plot_embeddings.py --stage before`只重新着色 | 时机和用途相同 |
| cell type注释 | 合并cell type、donor和disease注释 | 合并SCANPY cell type、sample和IR/NR注释 | 都用于解释和绘图，不作为核心MethylVI监督标签 |
| MethylVI输入 | 从ALLC重建整数`mc/cov`并写入H5MU，约14 GB | 从ALLC重建整数`mc/cov`并写入H5MU，约1.03 GiB | 构建原则相同，数据规模不同 |
| batch key | `donor` | `sample_id` | 都去除个体/样本差异；当前`sample_id`与IR/NR完全绑定，需警惕过校正 |
| MethylVI训练 | **有**；最大500 epochs、batch size 32、hidden 128×1、early stopping、seed 0 | **有**；同核心参数，明确使用beta-binomial和region dispersion | 两者都训练MethylVI；元培的03源码缺失，不能确认所有内部参数逐项相同 |
| 实际停止 | epoch 73 | 89条训练记录（epoch 0–88） | 数据和收敛过程不同，不是方法变化 |
| 校正后latent | 20维donor-corrected latent | 20维`X_methylVI` | 两者都产生20维MethylVI latent |
| 校正后UMAP | 03训练阶段在MethylVI latent上计算；04/05重新着色 | 09训练阶段在`X_methylVI`上计算；10/11重新着色 | 时机和主要实现相同 |
| 训练产物 | README记录model、latent、embedding和训练输出 | model、latent、embedding H5AD、history和run summary | 两者都保存训练结果；当前流程的产物命名和验收记录更明确 |
| supervised UMAP | 可选，标签权重0.2、0.5、0.7、0.9、1.0，不属于核心校正embedding | 可选，本次设置0.2、0.5、0.7、0.9，同样不属于核心校正embedding | 两者目标相同；当前流程不运行1.0 |
| 校正评估 | 有校正前后绘图；README未记录统一的定量前后指标 | 当前按分析需求暂不进行定量batch评估 | 两者当前都不以统一定量评估作为本阶段输出 |

这张表的关键修正是：`yuanpei`的“MethylVI训练”、“校正后latent”和“训练产物”都不是“无”。它与当前流程的核心差异不在于是否训练MethylVI，而在于数据集、batch字段、QC来源、软件实现证据和附加验收步骤。

## 6. MethylVI核心参数对照

| 参数 | `yuanpei` | 当前流程 | 是否相同 |
|---|---:|---:|---|
| cells | 10,488 | 6,199 | 否，数据集不同 |
| retained 5-kb bins | 272,521 | 231,648 | 否，特征集不同 |
| context | mCG | mCG/CGN | 本质相同 |
| batch key | `donor` | `sample_id` | 语义相同，字段名不同 |
| latent维度 | 20 | 20 | 是 |
| hidden维度 | 128 | 128 | 是 |
| hidden层数 | 1 | 1 | 是 |
| 最大epochs | 500 | 500 | 是 |
| early stopping | 开启 | 开启 | 是 |
| 实际停止 | epoch 73 | 89条训练记录，epoch 0–88 | 否，数据和收敛不同 |
| batch size | 32 | 32 | 是 |
| neighbors | 15 | 15 | 是 |
| Leiden resolution | 1.0 | 1.0 | 是 |
| seed | 0 | 0 | 是 |
| CPU threads | 50 | 32（正式训练） | 否，资源不同 |

从参数上看，当前流程明显沿用了`yuanpei`的MethylVI核心设置。不同的停止轮数不是方法变化，而是不同细胞、bins、数据分布和运行实现产生的收敛差异。

## 7. 数据构建本质是否相同

### 7.1 相同原则

两者都明确指出：

```text
ALLCools H5AD中的X是经过处理的hypo-score，不能作为MethylVI原始计数。
```

两者都使用H5AD来确定：

1. 进入分析的细胞；
2. ALLCools最终保留的5-kb bins。

然后重新读取对应ALLC，在这些bins内聚合：

```text
mc  = 甲基化计数
cov = 总覆盖计数
```

最后构建包含`mc/cov`的H5MU。这一原则是两套MethylVI流程最核心的共同点。

### 7.2 当前流程能逐行确认的实现

当前`05_build_methylvi_input.py`会：

1. 解析H5AD中的5-kb坐标；
2. 为每个细胞匹配唯一ALLC；
3. 并行聚合CGN位点到保留bins；
4. 按细胞保存压缩`.npz`检查点；
5. 验证`mc ≤ cov`；
6. 根据最大coverage自动选择dtype；
7. 本次实际选择`uint16`；
8. 组装6,199×231,648的`mc/cov`层；
9. 合并样本、condition和cell type注释；
10. 写出并回读验证H5MU。

### 7.3 `yuanpei`能确认但不能逐行检查的部分

`yuanpei` README明确说build会从10,488个ALLC重新聚合`mc/cov`，并通过`methylVI/count_rows/`断点复用。因此设计原则与当前08一致。

但`02_build_methylvi_input.py`源码缺失，不能确认：

- dtype选择是否完全相同；
- ALLC区域边界处理是否完全相同；
- context匹配是否完全相同；
- memmap和H5MU写入方法是否完全相同；
- checkpoint manifest和溢出保护是否与当前08一致。

所以应写成“数据构建思路相同”，而不是“实现逐行相同”。

## 8. batch correction本质是否相同

### 8.1 batch变量语义相同

`yuanpei`：

```text
batch_key = donor
```

当前流程：

```text
batch_key = sample_id
```

两者都试图去除不同个体/样本之间的系统差异。在本项目中每个`sample_id`对应一个生物样本，因此其统计角色与`yuanpei`的`donor`基本一致。

### 8.2 不能确认的内部差异

当前09源码明确使用：

```text
scvi.external.METHYLVI
likelihood = betabinomial
dispersion = region
```

`yuanpei`版本锁同时记录：

```text
methyl-vi==0.1.0
scvi-tools==1.3.3
```

但其`03_train_methylvi_donor_batch.py`源码缺失，所以无法从当前文件确认：

- 它导入独立`methyl-vi`包还是`scvi.external.METHYLVI`；
- likelihood是否为`betabinomial`；
- dispersion是否为`region`；
- early-stopping patience和其他内部训练参数是否完全相同。

因此最稳妥的结论是：

> batch correction目标、输入层、batch变量类型和主要网络参数相同；具体Python类和未列出的内部参数尚不能证明完全一致。

若要完成严格的逐行比较，需要补充服务器上的：

```text
/home/lijia/jiangyuanpei/.../methylVI/02_build_methylvi_input.py
/home/lijia/jiangyuanpei/.../methylVI/03_train_methylvi_donor_batch.py
```

## 9. ALLCools上游比较

当前`03_cluster_allcools.py`与`yuanpei/cluster_5kbin.py`可以逐行比较。除了默认线程数和图片输出位置外，算法参数一致：

| 环节 | 共同参数 |
|---|---|
| MCDS score | `CGN`、`hypo-score` |
| binarize cutoff | 0.95 |
| LSI | `arpack`、seed 0 |
| significant PC | `p_cutoff=0.1` |
| neighbors | 25 |
| 初始Leiden | resolution 1.0 |
| t-SNE | Euclidean、perplexity 30、exaggeration -1 |
| Consensus repeats | 500 |
| Consensus resolution | 0.5 |
| min cluster size | 10 |
| consensus rate | 0.5 |
| train fraction | 0.5，最多500个细胞 |
| max iterations | 20 |

03阶段的实质区别只有：

- `yuanpei`默认50线程，当前默认32线程；
- 图片文件名和保存位置不同。

所以ALLCools聚类核心可以认定为同一种实现。

## 10. 两套流程的主要工程差异

| 比较项 | `yuanpei` | 当前流程 |
|---|---|---|
| 项目设计 | donor + disease | 5 IR + 5 NR |
| 调度器 | Slurm `sbatch` | `dsub` |
| 计算节点 | 固定`cu03` | 由dsub调度，本次为node-11 |
| CPU/内存 | 50 CPU、250 GB | 正式训练32 CPU、64 GiB |
| BLAS线程 | 各变量设为50 | 多进程阶段内部BLAS固定为1 |
| 原始细胞QC | 复现目录未记录MethSCAn provenance硬检查 | 核验300k QC和10份provenance |
| 样本数量硬检查 | donor可赋值即可 | 强制10样本、5 IR、5 NR |
| cell ID | donor前缀，规范化`-/_` | `sample_id__barcode` |
| 选中细胞 | 10,488 | 6,199 |
| bins | 272,521 | 231,648 |
| 注释匹配 | 9,390完整，1,098仅推断donor | 20260810新注释已验证：5,765匹配，434 Unannotated，17个Exclude |
| batch字段 | donor | sample_id |
| 生物学字段 | disease | condition（IR/NR） |
| supervised UMAP | 有，可选，weights 0.2–1.0 | 有，可选，本次weights 0.2、0.5、0.7、0.9 |
| 定量batch评估 | README未记录 | 当前暂不运行 |
| H5MU | README记录约14 GB | 约1.03 GiB |
| 版本记录 | MuData 0.3.2、Torch 2.5.1.post7 | MuData 0.3.9、Torch 2.12.1 |
| 模型入口 | 缺失源码，锁文件含`methyl-vi` | 明确为`scvi.external.METHYLVI` |

## 11. 注释是否参与模型训练

两套流程都说明核心模型不使用cell type或疾病/condition标签：

- `yuanpei`用donor作batch；cell type和disease仅用于绘图。
- 当前流程用sample_id作batch；cell type和IR/NR不传入MethylVI模型。

因此核心latent不是监督式cell type embedding。

`yuanpei`另外提供可选的supervised UMAP，并测试权重：

```text
0.2, 0.5, 0.7, 0.9, 1.0
```

这一步使用标签引导可视化，但README明确说它不是核心donor correction分析。当前流程将它作为独立的`supervised`阶段，基于已训练的`X_methylVI`使用cell type标签和0.2、0.5、0.7、0.9四个权重重算UMAP。这些图必须标注为“监督式UMAP”，不应与`06_train_methylvi.py`产生的无监督MethylVI UMAP混淆。

`yuanpei`的`run_supervised_umap.py`源码不在当前副本中，因此只能确认它接收0.2、0.5、0.7、0.9和1.0这些weights，不能证明其他UMAP参数与当前`08_plot_supervised_umap.py`逐项相同。当前脚本已显式固定`n_neighbors=15`、`min_dist=0.5`、`metric=euclidean`、`target_metric=categorical`和`seed=0`。

## 12. 两边结果能否直接比较

不能直接把两个项目的UMAP形状进行一一比较，因为：

- 细胞数不同；
- 保留bins不同；
- donor/sample构成不同；
- disease与IR/NR不是同一标签；
- 环境版本和线程资源不同；
- `yuanpei`训练源码缺失，部分内部参数未确认。

可以比较的是方法学框架和相同定义的质量指标，例如：

- batch silhouette；
- donor/sample邻域混合度；
- cell type silhouette和邻域纯度；
- 校正前后生物学结构保留；
- early stopping和loss曲线。

当前项目内部的有效对照应始终使用相同6,199个细胞：

```text
校正前：ALLCools H5AD的X_pca
校正后：MethylVI H5AD的X_methylVI
```

## 13. 当前项目特有的统计限制

当前`sample_id`与IR/NR condition完全绑定，因此batch与condition不是独立变量。使用`sample_id`校正可能同时削弱真实IR/NR差异。

`yuanpei`是否存在同样程度的donor/disease混杂，当前审计文件没有提供donor与disease交叉表，因此不能判断。

当前结果必须同时检查：

1. sample mixing是否改善；
2. cell type结构是否保留；
3. 各cell type内部的sample mixing；
4. 各cell type内部的IR/NR差异；
5. 不能把IR/NR完全混合作为唯一成功标准。

## 14. 最终判断

### 方法层面

两者本质相同。它们使用相同的核心思想：

```text
ALLCools保留5-kb bins
→ 从ALLC重建mc/cov
→ donor/sample作为batch covariate
→ 20维MethylVI latent
→ 邻居图、UMAP和Leiden
```

### 代码层面

- ALLCools聚类脚本可以确认核心代码和参数几乎完全相同。
- MethylVI核心脚本在`yuanpei`副本中缺失，因此只能确认接口设计和参数，不能确认逐行实现一致。

### 数据和工程层面

两者不完全相同。当前流程针对5 IR + 5 NR数据重写了样本识别、300k QC、输入审计、并行方式、断点和版本保护，并将supervised UMAP实现为独立的可选阶段。

最准确的一句话是：

> 当前流程是`yuanpei` MethylVI核心方法在当前IR/NR数据上的复现和增强版；统计主线相同，但不是对原脚本、数据和环境的原样复制。
