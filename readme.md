This is a snakemake workflow to run short read mapping and subsequent variant calling with gatk HaplotypeCaller.

Choice of dragen or bwa for mapping algorithm.

## Dependencies
 - snakemake
 - conda
 - dragen
 - bwa
 - samtools
 - snakemake

First clone the repo.
```
git clone https://github.com/gconcepcion/shortreadvariantcalling.git ./workflow
```

Next, prep the workspace by linking a reference and your input data. The reference directory must consist of a reference.fasta, it's *.fasta.fai index and reference.dict hash tables

If you don't have one already, first create a sequence dictionary
```
gatk CreateSequenceDictionary -R /path/to/reference/ref.fasta
```

Next, prep the workflow directory
```
cd workflow
ln -s /path/to/references ./
mkdir -p data/$SAMPLE && 
```

### Edit your config.yaml & reference.yaml accordingly
Edit your config.yaml for the sample being run and make any necessary adjustments.
Edit your reference.yaml and ensure necessary indices are in place


### Process your input sample data by specifying local and number of threads or selecting distributed.
```
bash -ex process_flowcells.sh 24 
```

### Next call all variants in a sample by providing sample name and requested threads for variant calling.
``` 
bash -ex process_samples.sh 24
```
