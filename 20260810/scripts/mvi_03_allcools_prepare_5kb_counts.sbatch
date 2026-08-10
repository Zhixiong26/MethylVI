#!/usr/bin/env bash
#SBATCH --job-name=allcools_5kb_counts
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=250G
#SBATCH --time=7-00:00:00
#SBATCH --output=allcools_5kb_counts.%j.out
#SBATCH --error=allcools_5kb_counts.%j.err

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
INPUT=${1:-${MVI_COV_DIR:-}}
OUTPUT=${2:-${MVI_ALLCOOLS_OUTPUT:-}}
CHROM_SIZES=${3:-${MVI_CHROM_SIZES:-}}
PYTHON_BIN=${PYTHON_BIN:-python}
ALLCOOLS_EXE=${ALLCOOLS_EXE:-allcools}
THREADS=${SLURM_CPUS_PER_TASK:-50}
INPUT_ALLC=${OUTPUT}/input_allc
ALLC_TABLE=${OUTPUT}/selected_cells.allc.tsv
MCDS=${OUTPUT}/mcg_5kb.mcds

export OMP_NUM_THREADS=${THREADS}
export MKL_NUM_THREADS=${THREADS}
export OPENBLAS_NUM_THREADS=${THREADS}
export NUMEXPR_NUM_THREADS=${THREADS}
export MPLBACKEND=Agg

[[ -n "$INPUT" && -d "$INPUT" ]] || { echo "ERROR: pass a cov directory as arg 1 or set MVI_COV_DIR" >&2; exit 1; }
[[ -n "$OUTPUT" ]] || { echo "ERROR: pass an output directory as arg 2 or set MVI_ALLCOOLS_OUTPUT" >&2; exit 1; }
[[ -s "$CHROM_SIZES" ]] || { echo "ERROR: chrom sizes not found: $CHROM_SIZES" >&2; exit 1; }
mkdir -p "$OUTPUT" "$INPUT_ALLC"

echo "[$(date)] host=$(hostname) input=$INPUT output=$OUTPUT threads=$THREADS"
{
    date
    "$PYTHON_BIN" --version
    "$ALLCOOLS_EXE" --version || true
    "$PYTHON_BIN" -c 'import ALLCools,anndata,numpy,pandas,scanpy,scipy,sklearn; print("ALLCools",ALLCools.__version__); print("anndata",anndata.__version__); print("numpy",numpy.__version__); print("pandas",pandas.__version__); print("scanpy",scanpy.__version__); print("scipy",scipy.__version__); print("sklearn",sklearn.__version__)'
} > "$OUTPUT/software_versions.txt" 2>&1

# Expected cov columns (no header): chrom, start, end, ratio_percent,
# methylated_count, unmethylated_count. Each CpG coordinate must occur once.
cov_to_allc() {
    local cov=$1 cell_id out tmp
    cell_id=$(basename "$cov")
    cell_id=${cell_id%.cov.gz}
    cell_id=${cell_id%.cov}
    out="$INPUT_ALLC/${cell_id}.allc.tsv.gz"
    if [[ -s "$out" && -s "$out.tbi" ]]; then
        return
    fi
    tmp="${out}.tmp.${SLURM_JOB_ID:-$$}"
    { if [[ "$cov" == *.gz ]]; then gzip -cd -- "$cov"; else cat -- "$cov"; fi; } | awk -v OFS='\t' '
        BEGIN { previous="" }
        {
            key=$1 OFS $2
            if (key == previous) {
                print "ERROR: duplicate CpG coordinate in " FILENAME ": " key > "/dev/stderr"
                exit 2
            }
            if (NF < 6 || $5 !~ /^[0-9]+$/ || $6 !~ /^[0-9]+$/) {
                print "ERROR: invalid cov row in " FILENAME " at line " NR > "/dev/stderr"
                exit 2
            }
            print $1, $2, "+", "CGN", $5, $5+$6, 1
            previous=key
        }
    ' | bgzip -@ 1 -c > "$tmp"
    mv "$tmp" "$out"
    tabix -f -s 1 -b 2 -e 2 "$out"
}
export INPUT_ALLC
export -f cov_to_allc
find "$INPUT" -maxdepth 1 \( -name '*.cov' -o -name '*.cov.gz' \) -print0 \
    | sort -z \
    | xargs -0 -r -n 1 -P "$THREADS" bash -euo pipefail -c 'cov_to_allc "$1"' _

find "$INPUT_ALLC" -maxdepth 1 -type f -name '*.allc.tsv.gz' -print0 \
    | sort -z \
    | awk -v RS='\0' -v OFS='\t' '{p=$0; sub(/^.*\//, "", $0); sub(/\.allc\.tsv\.gz$/, "", $0); print $0, p}' \
    > "$ALLC_TABLE"

n_cov=$(find "$INPUT" -maxdepth 1 \( -name '*.cov' -o -name '*.cov.gz' \) | wc -l)
n_allc=$(wc -l < "$ALLC_TABLE")
n_unique=$(cut -f1 "$ALLC_TABLE" | sort -u | wc -l)
if (( n_cov == 0 || n_cov != n_allc || n_allc != n_unique )); then
    echo "ERROR: input verification failed: cov=$n_cov allc=$n_allc unique=$n_unique" >&2
    exit 1
fi
echo "[$(date)] verified cells=$n_allc"

if [[ ! -e "$OUTPUT/mcds.COMPLETE" ]]; then
    "$ALLCOOLS_EXE" generate-dataset \
        --allc_table "$ALLC_TABLE" \
        --output_path "$MCDS" \
        --chrom_size_path "$CHROM_SIZES" \
        --obs_dim cell \
        --cpu "$THREADS" \
        --chunk_size 10 \
        --regions chrom5k 5000 \
        --quantifiers chrom5k count CGN \
        --quantifiers chrom5k hypo-score CGN cutoff=0.9
    touch "$OUTPUT/mcds.COMPLETE"
else
    echo "[$(date)] existing completed MCDS detected; skipping generation"
fi

"$PYTHON_BIN" "$SCRIPT_DIR/mvi_04_allcools_cluster_5kb.py" \
    --mcds "$MCDS" \
    --output "$OUTPUT" \
    --threads "$THREADS"
echo "[$(date)] completed"
