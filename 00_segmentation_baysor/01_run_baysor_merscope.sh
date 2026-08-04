# Generic per-sample template: Baysor segmentation on raw MERSCOPE transcripts.
# Sample-specific values (region directory, SLURM resources) differ -- see
# README.md for the full per-sample table -- fill in below before submitting.


SAMPLE="sample_id"
REGION_DIR="/path/to/merscope/region"   # dir containing detected_transcripts.csv
JULIA_NUM_THREADS=40

msg () { echo "#[$(date "+%Y-%m-%dT%H:%M:%S")] $*"; }

cd /path/to/analysis/cell_segmentation/baysor

msg "======== ${SAMPLE} ========"
mkdir -p output_merscope_${SAMPLE}

JULIA_NUM_THREADS=${JULIA_NUM_THREADS} baysor run \
  ${REGION_DIR}/detected_transcripts.csv \
  :cell_id \
  -c ./merscope.toml \
  -o ./output_merscope_${SAMPLE} \
  -p \
  --polygon-format=GeometryCollectionLegacy

msg "Done: ${SAMPLE}"
