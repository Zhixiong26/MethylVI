#!/usr/bin/env bash

# ============================================================================
# MethylVI 项目统一配置文件
#
# 入口脚本 mvi_05_run_pipeline.sh 会自动加载本文件。
# 如需更换数据集，可在加载本文件前覆盖对应的 MVI_* 环境变量。
# ============================================================================

# ----------------------------------------------------------------------------
# 1. 输入数据路径
# ----------------------------------------------------------------------------

# MethSCAn 上游数据根目录；如果实际数据在集群上，请改成对应路径。
export MVI_DATA_ROOT="${MVI_DATA_ROOT:-/path/to/your/methscan_project}"

# 当前脚本目录；用于定位配置、日志和公共模块。
export MVI_REPRO="${MVI_REPRO:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# MethSCAn 上游脚本和流程说明所在目录，仅用于记录数据来源。
export MVI_METHSCAN_UPSTREAM="${MVI_METHSCAN_UPSTREAM:-/Users/luozhixiong/Library/Mobile Documents/com~apple~CloudDocs/Documents/PHD/脚本/Methscan/01_Upstream}"

# ALLCools 生成并筛选的 5-kb 聚类 H5AD 文件。
export MVI_H5AD="${MVI_H5AD:-${MVI_DATA_ROOT}/mcg_5kb.clustered.h5ad}"

# 每个细胞一个 ALLC 文件的目录，文件名需要与 H5AD 的 cell ID 匹配。
export MVI_ALLC_DIR="${MVI_ALLC_DIR:-${MVI_DATA_ROOT}/input_allc}"

# 可选的逐细胞注释文件；为空时使用 cell ID 和样本元数据推断信息。
export MVI_ANNOTATION="${MVI_ANNOTATION:-}"

# 10 个样本的 sample_id/condition 元数据表。
export MVI_SAMPLE_METADATA="${MVI_SAMPLE_METADATA:-${MVI_REPRO}/mvi_01_sample_metadata.tsv}"

# ----------------------------------------------------------------------------
# 2. 输出路径
# ----------------------------------------------------------------------------

# MethylVI 项目输出根目录。
export MVI_ROOT="${MVI_ROOT:-${MVI_DATA_ROOT}/methylVI_results}"

# MethylVI 输入 H5MU，包含 mCG.layers['mc'] 和 mCG.layers['cov']。
export MVI_INPUT="${MVI_INPUT:-${MVI_ROOT}/methylvi_5kbin_input.h5mu}"

# 模型、latent、UMAP、Leiden 和训练记录的输出目录。
export MVI_RESULTS="${MVI_RESULTS:-${MVI_ROOT}/results_ir_nr}"

# 输入审计 JSON 报告。
export MVI_AUDIT="${MVI_AUDIT:-${MVI_REPRO}/mvi_input_audit.json}"

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

# 输入审计要求的样本总数和两组样本数量。
export MVI_EXPECTED_SAMPLES="${MVI_EXPECTED_SAMPLES:-10}"
export MVI_EXPECTED_IR="${MVI_EXPECTED_IR:-5}"
export MVI_EXPECTED_NR="${MVI_EXPECTED_NR:-5}"

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

# Conda 环境名称和初始化脚本路径。
export MVI_CONDA_ENV="${MVI_CONDA_ENV:-methVI}"
export MVI_CONDA_INIT="${MVI_CONDA_INIT:-}"

# CPU 线程数、内存记录值和 PyTorch 加速方式（auto/cpu/gpu）。
export MVI_THREADS="${MVI_THREADS:-50}"
export MVI_MEMORY_GB="${MVI_MEMORY_GB:-250}"
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
