#import re
from pathlib import Path

configfile: "reference.yaml"         # reference configuration
configfile: "config.yaml"            # config options

reffasta = config['ref']['fasta']
refbase = Path(reffasta).stem
#intervals = config['ref']['index']
sample = config['sample']
qual = config['bwafilter']

print(f"Reference: {reffasta}")
print(f"Sample: {sample}")

targets = []
abams = list(Path(f"sample/{sample}/mapping").glob("*_markdup.bam"))

abam_paths = [i.as_posix() for i in abams]
mosdepthout = [os.path.splitext(i)[0] + ".mosdepth.summary.txt" for i in abam_paths]

flowcells = [os.path.splitext(os.path.basename(i))[0] for i in abam_paths]

print(f"Using the following movies: {flowcells}")

#sample_ids, flowcell_ids = glob_wildcards("sbb_data/{sample}/{flowcell}")

include: 'rules/variants.smk'
# Call variants from bams

if 'variants' in config['sample_targets']:
    targets.extend([f"sample/{sample}/gatk/{sample}.raw.vcf.gz"])

if 'coverage' in config['sample_targets']:
    targets.extend([f"sample/{sample}/gatk/{sample}.depthofcov.txt"])

if 'filter' in config['sample_targets']:
    targets.extend([f"sample/{sample}/gatk/{sample}.SNPs.filtered.vcf.gz"])
    targets.extend([f"sample/{sample}/gatk/{sample}.INDELs.filtered.vcf.gz"])

if 'updateref' in config['sample_targets']:
    targets.extend([f"sample/{sample}/gatk/{sample}.SNPs.vcf.gz"])
    targets.extend([f"sample/{sample}/gatk/{sample}.INDELs.vcf.gz"])
    targets.extend([f"sample/{sample}/gatk/{refbase}_updated.fasta"])

print(targets)

rule all:
    input: targets
