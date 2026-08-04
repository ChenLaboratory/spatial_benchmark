# Generic per-sample template: Baysor segmentation on raw Xenium transcripts, then
# import the result back into a XeniumRanger bundle. 

SAMPLE="sample_id"
XENIUM_BUNDLE="/path/to/xenium/bundle"          # dir containing transcripts.parquet
JULIA_NUM_THREADS=50


msg () { echo "#[$(date "+%Y-%m-%dT%H:%M:%S")] $*"; }

# ── Step 1: Baysor segmentation ──────────────────────────────────────────────

msg "Running Baysor for ${SAMPLE}"
mkdir -p output_xenium_${SAMPLE}

JULIA_NUM_THREADS=${JULIA_NUM_THREADS} baysor run \
  ${XENIUM_BUNDLE}/transcripts.parquet \
  :cell_id \
  -c ./xenium.toml \
  -o ./output_xenium_${SAMPLE}/ \
  -p \
  --polygon-format=GeometryCollectionLegacy

