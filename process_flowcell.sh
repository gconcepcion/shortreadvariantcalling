#!/bin/bash
source /mnt/software/Modules/current/init/bash

module load snakemake/5.26.1

flowcell=$1
sample=$2
threads=$3
targetfile="mapping/${sample}/${flowcell}_markdup.bam"

snakemake --use-conda -j ${threads} -pr --snakefile process_flowcell.smk $targetfile
