#!/bin/bash
#source /mnt/software/Modules/current/init/bash
umask 002

. "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate shortreadvc
source variables.env

distributed=${1:-false}
threads=${2:-6}
createenv=${3:-false}

if [ $createenv = true ] ; then
    snakemake --conda-create-envs-only --conda-frontend mamba --use-conda -j ${threads} --snakefile process_flowcells.smk
else
    if [ $distributed = false ]; then
        snakemake --use-conda -j ${threads} -pr --snakefile process_flowcells.smk
    else
        snakemake --profile profile/slurm --snakefile process_flowcells.smk
    fi 
fi
