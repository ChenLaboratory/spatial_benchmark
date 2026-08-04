#!/bin/bash
#SBATCH --job-name=ct_run_ficture_ct
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=24:00:00


#--conda-frontend conda --use-conda
module load micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate /vast/projects/CDFC/rlyu/Softwares/condaenvs/snakemake7

snakemake -s 01_run_ficture.smk --rerun-incomplete -j 22 --rerun-trigger mtime --keep-going --retries 2 --cluster "sbatch --mem {resources.mem} --time=18:00:00 --cpus-per-task={resources.cpus} --output {params.slurm_out} --error {params.slurm_err} --partition=regular "
