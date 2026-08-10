#!/usr/bin/env bash
#SBATCH --job-name=mvi_5kb_ir_nr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=250G
#SBATCH --time=7-00:00:00
#SBATCH --output=mvi_5kb_ir_nr.%j.out
#SBATCH --error=mvi_5kb_ir_nr.%j.err

set -euo pipefail
PIPE=$(cd "$(dirname "$0")" && pwd)
bash "$PIPE/mvi_05_run_pipeline.sh" all
