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
STAGED_COV=${OUTPUT}/staged_cov
COV_MANIFEST=${OUTPUT}/staged_cov_manifest.tsv

export OMP_NUM_THREADS=${THREADS}
export MKL_NUM_THREADS=${THREADS}
export OPENBLAS_NUM_THREADS=${THREADS}
export NUMEXPR_NUM_THREADS=${THREADS}
export MPLBACKEND=Agg

[[ -n "$INPUT" && -d "$INPUT" ]] || { echo "ERROR: pass a cov directory as arg 1 or set MVI_COV_DIR" >&2; exit 1; }
[[ -n "$OUTPUT" ]] || { echo "ERROR: pass an output directory as arg 2 or set MVI_ALLCOOLS_OUTPUT" >&2; exit 1; }
[[ -s "$CHROM_SIZES" ]] || { echo "ERROR: chrom sizes not found: $CHROM_SIZES" >&2; exit 1; }
mkdir -p "$OUTPUT" "$INPUT_ALLC"
mkdir -p "$STAGED_COV"

echo "[$(date)] host=$(hostname) input=$INPUT output=$OUTPUT threads=$THREADS"
{
    date
    "$PYTHON_BIN" --version
    "$ALLCOOLS_EXE" --version || true
    "$PYTHON_BIN" -c 'import ALLCools,anndata,numpy,pandas,scanpy,scipy,sklearn; print("ALLCools",ALLCools.__version__); print("anndata",anndata.__version__); print("numpy",numpy.__version__); print("pandas",pandas.__version__); print("scanpy",scanpy.__version__); print("scipy",scipy.__version__); print("sklearn",sklearn.__version__)'
} > "$OUTPUT/software_versions.txt" 2>&1

# 将输入整理为平铺的 cov 目录。
# 输入既可以是一个平铺目录，也可以是包含多个 *_Met/cov_dedup_probability
# 子目录的数据根目录。嵌套输入会把样本前缀加入 cell ID，避免不同样本的
# 相同 barcode 互相覆盖，例如 IR01__AAAC...cov.gz。
stage_cov_inputs() {
    local manifest_tmp="$COV_MANIFEST.tmp.${SLURM_JOB_ID:-$$}"
    local cov sample_dir sample_id cell_id staged_name staged_path existing
    local nested=0
    : > "$manifest_tmp"

    local flat_cov
    flat_cov=$(find "$INPUT" -maxdepth 1 -type f \( -name '*.cov' -o -name '*.cov.gz' \) -print -quit)
    if [[ -n "$flat_cov" ]]; then
        while IFS= read -r -d '' cov; do
            cell_id=$(basename "$cov")
            staged_name="$cell_id"
            staged_path="$STAGED_COV/$staged_name"
            if [[ -L "$staged_path" ]]; then
                existing=$(readlink "$staged_path")
                [[ "$existing" == "$cov" ]] || {
                    echo "ERROR: staged cov link points to a different file: $staged_path" >&2
                    exit 1
                }
            elif [[ -e "$staged_path" ]]; then
                echo "ERROR: staged cov target already exists: $staged_path" >&2
                exit 1
            else
                ln -s "$cov" "$staged_path"
            fi
            printf '%s\t%s\n' "$staged_name" "$cov" >> "$manifest_tmp"
        done < <(find "$INPUT" -maxdepth 1 -type f \( -name '*.cov' -o -name '*.cov.gz' \) -print0 | sort -z)
    else
        nested=1
        while IFS= read -r -d '' cov; do
            sample_dir=$(basename "$(dirname "$(dirname "$cov")")")
            if [[ "$sample_dir" =~ ^[0-9]+_((IR|NR)[0-9]{2})_Met$ ]]; then
                sample_id="${BASH_REMATCH[1]}"
            else
                echo "ERROR: cannot infer IR/NR sample ID from cov path: $cov" >&2
                exit 1
            fi
            cell_id=$(basename "$cov")
            cell_id=${cell_id%.cov.gz}
            cell_id=${cell_id%.cov}
            staged_name="${sample_id}__${cell_id}.cov.gz"
            staged_path="$STAGED_COV/$staged_name"
            if [[ -L "$staged_path" ]]; then
                existing=$(readlink "$staged_path")
                [[ "$existing" == "$cov" ]] || {
                    echo "ERROR: staged cov link points to a different file: $staged_path" >&2
                    exit 1
                }
            elif [[ -e "$staged_path" ]]; then
                echo "ERROR: staged cov target already exists: $staged_path" >&2
                exit 1
            else
                ln -s "$cov" "$staged_path"
            fi
            printf '%s\t%s\n' "$staged_name" "$cov" >> "$manifest_tmp"
        done < <(find "$INPUT" -type f -path '*/cov_dedup_probability/*.cov.gz' -print0 | sort -z)
    fi

    if [[ ! -s "$manifest_tmp" ]]; then
        echo "ERROR: no cov files found in flat input or nested cov_dedup_probability directories: $INPUT" >&2
        exit 1
    fi
    if [[ "$(cut -f1 "$manifest_tmp" | sort -u | wc -l)" -ne "$(wc -l < "$manifest_tmp")" ]]; then
        echo "ERROR: duplicate staged cell IDs detected" >&2
        exit 1
    fi
    mv "$manifest_tmp" "$COV_MANIFEST"
    if (( nested == 1 )); then
        echo "[$(date)] staged nested sample cov inputs: $(wc -l < "$COV_MANIFEST") cells"
    else
        echo "[$(date)] staged flat cov inputs: $(wc -l < "$COV_MANIFEST") cells"
    fi
}

stage_cov_inputs

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
find "$STAGED_COV" -maxdepth 1 -type l \( -name '*.cov' -o -name '*.cov.gz' \) -print0 \
    | sort -z \
    | xargs -0 -r -n 1 -P "$THREADS" bash -euo pipefail -c 'cov_to_allc "$1"' _

n_cov=$(wc -l < "$COV_MANIFEST")
: > "$ALLC_TABLE.tmp.${SLURM_JOB_ID:-$$}"
while IFS=$'\t' read -r staged_name _; do
    cell_id=${staged_name%.cov.gz}
    cell_id=${cell_id%.cov}
    allc_path="$INPUT_ALLC/${cell_id}.allc.tsv.gz"
    [[ -s "$allc_path" && -s "$allc_path.tbi" ]] || {
        echo "ERROR: missing converted ALLC for staged cell: $staged_name" >&2
        exit 1
    }
    printf '%s\t%s\n' "$cell_id" "$allc_path" >> "$ALLC_TABLE.tmp.${SLURM_JOB_ID:-$$}"
done < "$COV_MANIFEST"
mv "$ALLC_TABLE.tmp.${SLURM_JOB_ID:-$$}" "$ALLC_TABLE"
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
