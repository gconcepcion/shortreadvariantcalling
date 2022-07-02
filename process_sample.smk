import re
from pathlib import Path

configfile: "reference.yaml"         # reference configuration
configfile: "config.yaml"            # config options

reffasta = config['ref']['fasta']
refbase = Path(reffasta).stem
intervals = config['ref']['index']
sample = config['sample']
qual = config['bwafilter']

print(f"Reference: {reffasta}")
print(f"Sample: {sample}")

targets = []
abams = list(Path(f"mapping/{sample}").glob("*_markdup.bam"))

abam_paths = [i.as_posix() for i in abams]
mosdepthout = [os.path.splitext(i)[0] + ".mosdepth.summary.txt" for i in abam_paths]

flowcells = [os.path.splitext(os.path.basename(i))[0] for i in abam_paths]

print(flowcells)
#print(expand(f"mapping/{sample}/mosdepth/{{flowcell}}.mosdepth.summary.txt", flowcell=flowcells))

print(f"Using the following movies: {flowcells}")

sample_ids, flowcell_ids = glob_wildcards("sbb_data/{sample}/{flowcell}")

include: 'rules/variants.smk'
# Call variants from bams

if 'variants' in config['sample_targets']:
    targets.extend([f"gatk/{sample}/{sample}.raw.vcf.gz"])

if 'coverage' in config['sample_targets']:
    targets.extend([f"gatk/{sample}/{sample}.depthofcov.txt"])

if 'filter' in config['sample_targets']:
    targets.extend([f"gatk/{sample}/{sample}.SNPs.filtered.vcf.gz"])
    targets.extend([f"gatk/{sample}/{sample}.INDELs.filtered.vcf.gz"])

if 'updateref' in config['sample_targets']:
    targets.extend([f"gatk/{sample}/{sample}.SNPs.vcf.gz"])
    targets.extend([f"gatk/{sample}/{sample}.INDELs.vcf.gz"])
    targets.extend([f"gatk/{sample}/{refbase}_updated.fasta"])
#    targets.extend([f"gatk/{sample}/{sample}_updated.fasta"])

print(targets)
#for i in targets:
#    print(i)

rule all:
    input: targets
