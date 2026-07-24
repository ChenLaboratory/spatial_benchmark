
# Usage: sbatch 02_run_baysor_to_seu.sh <xenium|merscope>
# Converts Baysor's segmentation.csv to a Seurat object (see generate_seu.R).

PLATFORM=${1:?Usage: sbatch 02_run_baysor_to_seu.sh <xenium|merscope>}

msg () {
  datestr=$(date "+%Y-%m-%dT%H:%M:%S")
  echo "#[$datestr] $*"
}

SCRIPT_DIR=/path/to/analysis/cell_segmentation/baysor/${PLATFORM}
MTX_DIR=${SCRIPT_DIR}/results/mtx
BAYSOR_DIR=${SCRIPT_DIR}

declare -A SAMPLES_MAP
SAMPLES_MAP[xenium]="ER0114 ER0360 MH0007 MH0026 TN0177 TN0554"
SAMPLES_MAP[merscope]="ER0114 ER0360 MH0007 MH0026 TN0177"

read -ra SAMPLES <<< "${SAMPLES_MAP[${PLATFORM}]}"

MAP_SCRIPT=map_transcripts_baysor_xenium.py
if [ "${PLATFORM}" == "merscope" ]; then
  MAP_SCRIPT=map_transcripts_baysor_merscope.py
fi

CONF_CUTOFF=0        # 0 = keep all assigned transcripts; 0.9 = high-confidence only


# ── Step 1: Python — segmentation.csv → MTX matrix ──────────────────────────

for SAMPLE in "${SAMPLES[@]}"; do
  
  BAYSOR_CSV=${BAYSOR_DIR}/output_${PLATFORM}_${SAMPLE}/segmentation.csv
  OUT=${MTX_DIR}/${SAMPLE}_conf${CONF_CUTOFF}

  python3 ${SCRIPT_DIR}/${MAP_SCRIPT} \
    -baysor      ${BAYSOR_CSV} \
    -out         ${OUT} \
    -conf_cutoff ${CONF_CUTOFF}

done

# ── Step 2: R — MTX + cell stats → Seurat object (.rds) ─────────────────────

Rscript /path/to/analysis/cell_segmentation/baysor/generate_seu.R ${PLATFORM} ${CONF_CUTOFF}

