This is a snakemake workflow to run short read mapping and subsequent variant calling with gatk HaplotypeCaller.

Choice of dragen or bwa for mapping algorithm.

## Dependencies
 - snakemake
 - conda
 - dragen
 - bwa
 - samtools
 - 

Create a directory for the analysis and clone this repo into it:
```git clone ...
cd shortreadvariantcalling
```

### Link your reference files
```
ln -s /path/to/reference ./
```

### Edit your config.yaml & reference.yaml according to your layout
```
vi reference.yaml
vi config.yaml
```

### Create directory for your input fastqs and copy or link the data
```
mkdir data/chm13
ln -s /path/to/data/*.fastq.gz ./
```
### Initialize conda environments so there are no collisions when processing flowcells in parallel
```
snakemake -j 4 --conda-create-envs-only --use-conda --snakefile process_flowcell.smk mapping/chm13/XR0000000-BCC_L01_R1_Sample_Library_markdup.bam
```

### Process your input data by providing flowcell, sample name and number of requested cores for mapping
```
bash -ex process_flowcell.sh XR0000000-BCC_L01_R1_Sample_Library chm13 24
```

### Next process samples by providing sample name and requested threads:
``` 
bash -ex process_samples.sh chm13 24
```

