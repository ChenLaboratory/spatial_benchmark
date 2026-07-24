
#!/bin/bash

module load micromamba/ 
eval "$(micromamba shell hook --shell bash)"
micromamba activate /vast/projects/CDFC/rlyu/Softwares/condaenvs/ficture
DIR="$(cd "$(dirname "$0")" && pwd)"
path=$1
s_n=$2
mu_scale=1 # If your data's coordinates are already in micrometer
key=Count
MJ=X # If your data is sorted by the Y-axis
env=/vast/projects/CDFC/rlyu/Softwares/condaenvs/ficture
gitpath=../../../../../Softwares/ficture/ # location of ficture
#SLURM_ACCOUNT= # For submitting jobs to slurm

input=${path}/${s_n}_roi.tsv.gz
output=${path}/batched.matrix.tsv.gz
rec=$(sbatch --job-name=${s_n}_vz1  --partition=regular --time=24:00:00 --cpus-per-task=1 ${DIR}/generic_I.sh input=${input} output=${output} MJ=${MJ} gitpath=${gitpath})
IFS=' ' read -ra ADDR <<< "$rec"
jobid1=${ADDR[3]}

# env=${env} 

nFactor=20 # Number of factors
sliding_step=2
train_nEpoch=5
train_width=12 # \sqrt{3} x the side length of the hexagon (um)
model_id=nF${nFactor}.d_${train_width} # An identifier kept in output file names
min_ct_per_feature=20 # Ignore genes with total count \< 20
R=10 # We use R random initializations and pick one to fit the full model
thread=16 # Number of threads to use
feature=${path}/feature.clean.tsv.gz



fit_width=12 # Often equal or smaller than train_width (um)
anchor_res=4 # Distance between adjacent anchor points (um)
radius=$(($anchor_res+1))
anchor_info=prj_${fit_width}.r_${anchor_res} # An identifier
coor=${path}/coordinate_minmax.tsv



#
# Prepare training minibatches, only need to run once if you plan to fit multiple models (say with different number of factors)
input=${path}/${s_n}_roi.tsv.gz
hexagon=${path}/${s_n}_roi_hexagon.d_${train_width}.tsv.gz
rec=$(sbatch --job-name=${s_n}_vz2  --partition=regular --dependency=afterok:${jobid1} --cpus-per-task=1 --time=24:00:00 ${DIR}/generic_II.sh env=${env} gitpath=${gitpath} key=${key} mu_scale=${mu_scale} major_axis=${MJ} path=${path} input=${input} output=${hexagon} width=${train_width} sliding_step=${sliding_step})
IFS=' ' read -ra ADDR <<< "$rec"
jobid2=${ADDR[3]}

# Model training
rec=$(sbatch --job-name=${s_n}_vz3  --partition=regular --dependency=afterok:${jobid2} --cpus-per-task=${thread} --time=24:00:00 ${DIR}/generic_III.sh env=${env} gitpath=${gitpath} key=${key} mu_scale=${mu_scale} major_axis=${MJ} path=${path} pixel=${input} hexagon=${hexagon} feature=${feature} model_id=${model_id} train_width=${train_width} nFactor=${nFactor} R=${R} train_nEpoch=${train_nEpoch} fit_width=${fit_width} anchor_res=${anchor_res} min_ct_per_feature=${min_ct_per_feature} thread=${thread})
IFS=' ' read -ra ADDR <<< "$rec"
jobid3=${ADDR[3]}

# --dependency=afterok:${jobid2}
# Pixel level decoding & visualization

rec=$(sbatch --job-name=${s_n}_vz4  --partition=regular --dependency=afterok:${jobid3} --cpus-per-task=${thread} --time=24:00:00  ${DIR}/generic_V.sh env=${env} gitpath=${gitpath} key=${key} mu_scale=${mu_scale} path=${path} model_id=${model_id} anchor_info=${anchor_info} radius=${radius} coor=${coor} pixel=${output} thread=${thread})
IFS=' ' read -ra ADDR <<< "$rec"
jobid4=${ADDR[3]}
#,${jobid1}