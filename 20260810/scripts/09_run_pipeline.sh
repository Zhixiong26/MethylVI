#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/00_config.sh"

# 多进程阶段由 MVI_THREADS 控制任务数；每个进程内部的数学库固定为单线程。
# 训练阶段会在 Python 内单独设置 PyTorch 线程数。
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
export MPLBACKEND=Agg
mkdir -p "$HERE/logs" "$MVI_ROOT" "$MVI_RESULTS"

usage() {
    cat <<'EOF'
用法：bash 09_run_pipeline.sh {prepare|verify|build|train|plots|supervised|test|all}

  prepare      整理 ALLC、生成 MCDS 并运行 ALLCools 5-kb 聚类
  verify       审计 H5AD、ALLC、样本信息和当前 SCANPY 注释
  build        构建包含整数 mc/cov 层的 MethylVI H5MU
  train        训练 MethylVI 并生成 latent、UMAP 和 Leiden 结果
  plots        同时重画校正前和校正后的普通嵌入图
  supervised   生成 target_weight=0.2、0.5、0.7、0.9 的监督式 UMAP
  test         运行公共函数单元测试和两轮 CPU smoke test
  all          依次运行 verify、build、train、plots、supervised（不含 prepare/test）
EOF
}

activate_methylvi() {
    if [[ "${MVI_SKIP_CONDA:-0}" == "1" ]]; then
        return
    fi
    if [[ -n "${MVI_CONDA_INIT:-}" && -f "$MVI_CONDA_INIT" ]]; then
        # shellcheck disable=SC1090
        source "$MVI_CONDA_INIT"
    elif command -v conda >/dev/null 2>&1; then
        local conda_base
        conda_base=$(conda info --base)
        # shellcheck disable=SC1090
        source "$conda_base/etc/profile.d/conda.sh"
    else
        echo "ERROR：未找到 conda；仅在已激活环境中才能设置 MVI_SKIP_CONDA=1" >&2
        exit 1
    fi
    conda activate "$MVI_CONDA_ENV"
}

stage=${1:-}
case "$stage" in
  prepare)
    # 02 脚本直接调用 MVI_ALLCOOLS_ENV，不激活 MethylVI 环境。
    bash "$HERE/02_prepare_allcools.sh" \
      "$MVI_DATA_ROOT" "$MVI_ALLCOOLS_OUTPUT" "$MVI_CHROM_SIZES" \
      2>&1 | tee "$HERE/logs/02_prepare_allcools.log"
    ;;
  verify)
    activate_methylvi
    python "$HERE/04_verify_inputs.py" \
      2>&1 | tee "$HERE/logs/04_verify_inputs.log"
    ;;
  build)
    activate_methylvi
    python "$HERE/05_build_methylvi_input.py" --threads "$MVI_THREADS" \
      2>&1 | tee "$HERE/logs/05_build_methylvi_input.log"
    ;;
  train)
    activate_methylvi
    NUMBA_NUM_THREADS="$MVI_THREADS" python "$HERE/06_train_methylvi.py" \
      --threads "$MVI_THREADS" --epochs "$MVI_MAX_EPOCHS" \
      --batch-size "$MVI_BATCH_SIZE" --accelerator "$MVI_ACCELERATOR" \
      2>&1 | tee "$HERE/logs/06_train_methylvi.log"
    ;;
  plots)
    activate_methylvi
    python "$HERE/07_plot_embeddings.py" --stage all \
      2>&1 | tee "$HERE/logs/07_plot_embeddings.log"
    ;;
  supervised)
    activate_methylvi
    NUMBA_NUM_THREADS="$MVI_THREADS" python "$HERE/08_plot_supervised_umap.py" \
      --threads "$MVI_THREADS" \
      2>&1 | tee "$HERE/logs/08_plot_supervised_umap.log"
    ;;
  test)
    activate_methylvi
    PYTHONPATH="$HERE${PYTHONPATH:+:$PYTHONPATH}" \
      python "$HERE/tests/test_mvi_utils.py" \
      2>&1 | tee "$HERE/logs/test_mvi_utils.log"
    python "$HERE/tests/test_methylvi_smoke.py" \
      2>&1 | tee "$HERE/logs/test_methylvi_smoke.log"
    ;;
  all)
    bash "$0" verify
    bash "$0" build
    bash "$0" train
    bash "$0" plots
    bash "$0" supervised
    ;;
  *)
    usage
    exit 2
    ;;
esac
