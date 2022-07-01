#!/bin/bash
source /mnt/software/Modules/current/init/bash

module load snakemake/5.26.1

flowcell=$1
sample=$2
threads=$3
createenv=${4-false}

if [ $createenv = true ] ; then
    snakemake --conda-create-envs-only --use-conda -j ${threads} --snakefile process_flowcell.smk
else
    snakemake --config "flowcell=$flowcell" --use-conda -j ${threads} -pr --snakefile process_flowcell.smk
fi
