#!/usr/bin/env bash

# ============================================================================
# MethylVI 项目统一配置文件
#
# 入口脚本 mvi_04_run_pipeline.sh 会自动加载本文件。
# 如需更换数据集，可在加载本文件前覆盖对应的 MVI_* 环境变量。
# ============================================================================

# ----------------------------------------------------------------------------
# 1. 输入数据路径
# ----------------------------------------------------------------------------

# MethSCAn 上游数据根目录。
export MVI_DATA_ROOT="${MVI_DATA_ROOT:-/share/LCZX_Data/data/allcools}"

# 当前脚本目录；用于定位配置、日志和公共模块。
export MVI_REPRO="${MVI_REPRO:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# 当前 MethylVI 项目目录（即 scripts 的上一级目录）。
export MVI_PROJECT_ROOT="${MVI_PROJECT_ROOT:-$(cd "${MVI_REPRO}/.." && pwd)}"

# MethSCAn 上游脚本和流程说明所在目录，仅用于记录数据来源。
export MVI_METHSCAN_UPSTREAM="${MVI_METHSCAN_UPSTREAM:-/Users/luozhixiong/Library/Mobile Documents/com~apple~CloudDocs/Documents/PHD/脚本/Methscan/01_Upstream}"

# 通过 MethSCAn 300k 细胞 QC 后的 ALLCools 5-kb 输出目录。
# 与包含全部 58,534 个原始细胞的审计目录 methylvi_5kb 分开，避免混用。
export MVI_ALLCOOLS_OUTPUT="${MVI_ALLCOOLS_OUTPUT:-${MVI_DATA_ROOT}/methylvi_5kb_300k}"

# ALLCools 生成并筛选的 5-kb 聚类 H5AD 文件。
export MVI_H5AD="${MVI_H5AD:-${MVI_ALLCOOLS_OUTPUT}/mcg_5kb.clustered.h5ad}"

# hg38 canonical chromosome sizes；供 ALLCools generate-dataset 使用。
export MVI_CHROM_SIZES="${MVI_CHROM_SIZES:-${MVI_REPRO}/hg38.canonical.chrom.sizes}"

# 每个细胞一个 ALLC 软链接的平铺目录，文件名与 H5AD cell ID 匹配。
export MVI_ALLC_DIR="${MVI_ALLC_DIR:-${MVI_ALLCOOLS_OUTPUT}/input_allc}"

# SCANPY 导出的全细胞注释表。公共读取器会将其 sample、group 和
# cell_type_integrated 标准化为 sample_id、condition 和 cell_type。
export MVI_ANNOTATION="${MVI_ANNOTATION:-/share/home/rzli/SCANPY/20260714/result/annotation/02_cell_annotation_all_cells.csv}"

# 10 个样本的 sample_id/condition 元数据表。
export MVI_SAMPLE_METADATA="${MVI_SAMPLE_METADATA:-${MVI_REPRO}/mvi_01_sample_metadata.tsv}"

# ----------------------------------------------------------------------------
# 2. 输出路径
# ----------------------------------------------------------------------------

# 300k QC 细胞的 MethylVI 项目输出根目录。
export MVI_ROOT="${MVI_ROOT:-${MVI_DATA_ROOT}/methylVI_results_300k}"

# MethylVI 输入 H5MU，包含 mCG.layers['mc'] 和 mCG.layers['cov']。
export MVI_INPUT="${MVI_INPUT:-${MVI_ROOT}/methylvi_5kbin_input.h5mu}"

# 模型、latent、UMAP、Leiden 和训练记录的输出目录。
export MVI_RESULTS="${MVI_RESULTS:-${MVI_ROOT}/results_ir_nr}"

# 所有图像的统一输出目录；与 scripts 并列，不写入数据目录。
export MVI_FIGURES_DIR="${MVI_FIGURES_DIR:-${MVI_PROJECT_ROOT}/result}"

# 按分析阶段区分校正前和校正后图像。
export MVI_FIGURES_BEFORE_DIR="${MVI_FIGURES_BEFORE_DIR:-${MVI_FIGURES_DIR}/01_before_methylvi}"
export MVI_FIGURES_AFTER_DIR="${MVI_FIGURES_AFTER_DIR:-${MVI_FIGURES_DIR}/02_after_methylvi}"

# 输入审计 JSON 报告；编号与生成它的06脚本对应。
export MVI_AUDIT="${MVI_AUDIT:-${MVI_REPRO}/mvi_06_input_audit.json}"

# ----------------------------------------------------------------------------
# 3. 样本信息和元数据字段
# ----------------------------------------------------------------------------

# 从 cell ID 中提取 sample_id 的正则表达式。
# 支持：25110891_IR01_Met__barcode、IR01__barcode、IR01_cell。
# 捕获结果为 IR01/NR01。
export MVI_SAMPLE_ID_REGEX="${MVI_SAMPLE_ID_REGEX:-^(?:[^_]+_)?((?:IR|NR)[0-9][0-9])(?:_Met)?(?:__|_)}"

# H5MU.obs 中样本名称和分组名称对应的列名。
export MVI_SAMPLE_KEY="${MVI_SAMPLE_KEY:-sample_id}"
export MVI_CONDITION_KEY="${MVI_CONDITION_KEY:-condition}"

# 输入审计要求的样本总数、两组样本数量和 300k QC 后细胞数。
export MVI_EXPECTED_SAMPLES="${MVI_EXPECTED_SAMPLES:-10}"
export MVI_EXPECTED_IR="${MVI_EXPECTED_IR:-5}"
export MVI_EXPECTED_NR="${MVI_EXPECTED_NR:-5}"
export MVI_EXPECTED_CELLS="${MVI_EXPECTED_CELLS:-6199}"

# MethSCAn 细胞 QC 白名单设置。300k 表示每个细胞至少覆盖 300,000 个
# CpG 位点；同时要求最多 10,000,000 个位点和 min_meth=55。
export MVI_USE_FILTERED_CELLS="${MVI_USE_FILTERED_CELLS:-1}"
export MVI_QC_TAG="${MVI_QC_TAG:-minmeth55_maxmethnone_maxsites10000000_covdedupprob}"
export MVI_FILTER_THRESHOLD="${MVI_FILTER_THRESHOLD:-300k}"
export MVI_FILTER_MIN_SITES="${MVI_FILTER_MIN_SITES:-300000}"
export MVI_FILTER_MAX_SITES="${MVI_FILTER_MAX_SITES:-10000000}"
export MVI_FILTER_MIN_METH="${MVI_FILTER_MIN_METH:-55}"
export MVI_FILTER_MAX_METH="${MVI_FILTER_MAX_METH:-none}"

# ----------------------------------------------------------------------------
# 4. 批次校正设置
# ----------------------------------------------------------------------------

# 默认按 sample_id 去除样本级 batch 效应。
# 如果存在独立的技术批次列，例如 technical_batch，可改为：
# export MVI_BATCH_KEY=technical_batch
export MVI_BATCH_KEY="${MVI_BATCH_KEY:-sample_id}"

# ----------------------------------------------------------------------------
# 5. 软件环境和计算资源
# ----------------------------------------------------------------------------

# MethylVI Conda 环境名称和 Conda 初始化脚本路径。
export MVI_CONDA_ENV="${MVI_CONDA_ENV:-methylvi}"
export MVI_CONDA_INIT="${MVI_CONDA_INIT:-/share/home/rzli/miniconda3/etc/profile.d/conda.sh}"

# 已验收的 ALLCools 独立环境路径。
export MVI_ALLCOOLS_ENV="${MVI_ALLCOOLS_ENV:-/share/home/rzli/miniconda3/envs/allcools}"

# 设为 1 时只整理和核验 ALLC 输入，不生成 MCDS；默认正常运行。
export MVI_STAGE_ONLY="${MVI_STAGE_ONLY:-0}"

# CPU 线程数、内存记录值和 PyTorch 加速方式（auto/cpu/gpu）。
export MVI_THREADS="${MVI_THREADS:-32}"
export MVI_MEMORY_GB="${MVI_MEMORY_GB:-190}"
export MVI_ACCELERATOR="${MVI_ACCELERATOR:-auto}"

# ----------------------------------------------------------------------------
# 6. MethylVI 模型和下游分析参数
# ----------------------------------------------------------------------------

# 训练批大小和最大训练轮数；训练启用 early stopping。
export MVI_BATCH_SIZE="${MVI_BATCH_SIZE:-32}"
export MVI_MAX_EPOCHS="${MVI_MAX_EPOCHS:-500}"

# 随机种子，保证可复现。
export MVI_SEED="${MVI_SEED:-0}"

# 甲基化计数设置：5-kb bin 和 mCG/CGN context。
export MVI_BIN_SIZE="${MVI_BIN_SIZE:-5000}"
export MVI_MC_CONTEXT="${MVI_MC_CONTEXT:-CGN}"

# MethylVI 网络结构：latent 维度、hidden 层维度和 hidden 层数。
export MVI_N_LATENT="${MVI_N_LATENT:-20}"
export MVI_N_HIDDEN="${MVI_N_HIDDEN:-128}"
export MVI_N_LAYERS="${MVI_N_LAYERS:-1}"

# latent 空间下游分析：邻居数和 Leiden 聚类分辨率。
export MVI_NEIGHBORS="${MVI_NEIGHBORS:-15}"
export MVI_LEIDEN_RESOLUTION="${MVI_LEIDEN_RESOLUTION:-1.0}"
