#!/bin/bash

. "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate shortreadvc
source variables.env


distributed=${1:-false}
threads=${2:-6}

if [ $distributed = false ]; then
    snakemake --use-conda -j ${threads} -pr --snakefile process_sample.smk
else
    snakemake --verbose --profile profile/slurm --snakefile process_sample.smk
fi 

