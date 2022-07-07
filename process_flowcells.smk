import os
import re

configfile: "reference.yaml"         # reference configuration
configfile: "config.yaml"         # reference configuration

mapper = config['mapper']
sample = config['sample']
adapter = config['adapter']
datadir = config['datadir']

ref = config['ref']['shortname']
refdict = config['ref']['dict']
reffasta = config['ref']['fasta']

print(f"Processing flowcells for Sample: {sample}")

targets = []

fastqs = list(Path(f"{datadir}/{sample}").glob("*.fastq.gz"))
print(fastqs)

### example flowcell name format: FB0026100-BCC_L02_R1
flowcell = re.compile(r'([A-Z][A-Z]\d+.*L\d+_R\d.*).fastq.gz')
flowcells = [flowcell.search(str(x)).group(1) for x in fastqs]
print(flowcells)

# align reads with dragmap and markduplicates with gatk
include: 'rules/alignment.smk'
# trimmed raw BAMS
if 'cutadapt' in config['flowcell_targets']:
    targets.extend(
        [f"sample/{sample}/cutadapt/{flowcell}_trimmed.{suffix}" for suffix in ['fastq.gz']] for flowcell in flowcells)

if 'seqtk' in config['flowcell_targets']:
    targets.extend([f"sample/{sample}/cutadapt/seqtk/{flowcell}_seqtkfqchk.txt"] for flowcell in flowcells)
    targets.extend([f"sample/{sample}/cutadapt/seqtk/{flowcell}_seqtkfqchk.txt.png"] for flowcell in flowcells)

if 'alignment' in config['flowcell_targets']:
    targets.extend([f"sample/{sample}/mapping/{flowcell}.bam"] for flowcell in flowcells)

if 'markdups' in config['flowcell_targets']:
    targets.extend(
        [f"sample/{sample}/mapping/{flowcell}_markdup.{suffix}" for suffix in ['bam', 'bam.bai'] for flowcell in flowcells])

if 'mosdepth' in config['flowcell_targets']:
    targets.extend([f"sample/{sample}/mapping/mosdepth/{flowcell}_markdup.{suffix}" for suffix in ['mosdepth.global.dist.txt',
                                                                                    'mosdepth.region.dist.txt',
                                                                                    'mosdepth.summary.txt',
                                                                                    'regions.bed.gz']
                                                                                    for flowcell in flowcells])

rule all:
    input: targets
