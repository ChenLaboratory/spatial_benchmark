# Usage: sbatch 04_run_cell_area_2d.sh <merscope|xenium>
# Computes 2D cell footprint area (xy_area_um2) from proseg boundary polygons
# and writes <sample>_proseg_xy_area.csv into the results directory. 

python3 zarr_cell_area_2d.py \
    -zarr proseg-output.zarr \
    -out  proseg_xy_area.csv
 
