#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/mvi_00_config.sh"

if [[ "${MVI_SKIP_CONDA:-0}" != "1" ]]; then
    if [[ -n "${MVI_CONDA_INIT:-}" && -f "$MVI_CONDA_INIT" ]]; then
        source "$MVI_CONDA_INIT"
    elif command -v conda >/dev/null 2>&1; then
        CONDA_BASE=$(conda info --base)
        source "$CONDA_BASE/etc/profile.d/conda.sh"
    else
        echo "ERROR: conda was not found; set MVI_SKIP_CONDA=1 only inside an activated environment" >&2
        exit 1
    fi
    conda activate "$MVI_CONDA_ENV"
fi

export OMP_NUM_THREADS="$MVI_THREADS"
export MKL_NUM_THREADS="$MVI_THREADS"
export OPENBLAS_NUM_THREADS="$MVI_THREADS"
export NUMEXPR_NUM_THREADS="$MVI_THREADS"
export NUMBA_NUM_THREADS="$MVI_THREADS"
export MPLBACKEND=Agg
mkdir -p "$HERE/logs" "$MVI_ROOT" "$MVI_RESULTS"

usage() {
    echo "Usage: bash mvi_05_run_pipeline.sh {verify|smoke|original-sample|build|train|plots|all}"
}

stage=${1:-}
case "$stage" in
  verify)
    python "$HERE/mvi_07_verify_inputs.py" 2>&1 | tee "$HERE/logs/mvi_07_verify_inputs.log"
    ;;
  smoke)
    python "$HERE/mvi_06_smoke_test.py" 2>&1 | tee "$HERE/logs/mvi_06_smoke_test.log"
    ;;
  original-sample)
    python "$HERE/mvi_08_plot_original_embedding.py" \
      2>&1 | tee "$HERE/logs/mvi_08_original_sample.log"
    ;;
  build)
    python "$HERE/mvi_09_build_input.py" --threads "$MVI_THREADS" \
      2>&1 | tee "$HERE/logs/mvi_09_build_input.log"
    ;;
  train)
    python "$HERE/mvi_10_train_model.py" \
      --threads "$MVI_THREADS" --epochs "$MVI_MAX_EPOCHS" \
      --batch-size "$MVI_BATCH_SIZE" --accelerator "$MVI_ACCELERATOR" \
      2>&1 | tee "$HERE/logs/mvi_10_train_model.log"
    ;;
  plots)
    python "$HERE/mvi_11_plot_celltype_sample.py" \
      2>&1 | tee "$HERE/logs/mvi_11_plot_celltype_sample.log"
    python "$HERE/mvi_12_plot_condition.py" \
      2>&1 | tee "$HERE/logs/mvi_12_plot_condition.log"
    ;;
  all)
    bash "$0" verify
    bash "$0" smoke
    bash "$0" original-sample
    bash "$0" build
    bash "$0" train
    bash "$0" plots
    ;;
  *) usage; exit 2 ;;
esac
