#!/bin/bash
source /mnt/software/Modules/current/init/bash

module load snakemake/5.26.1

sample=$1
threads=$2

snakemake --verbose --use-conda -j $threads -pr --snakefile process_sample.smk
