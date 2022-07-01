import re
from pathlib import Path

configfile: "reference.yaml"         # reference configuration
configfile: "config.yaml"            # config options

reffasta = config['ref']['fasta']
sample = config['sample']
qv = config['variantfilter']

print(f"Reference: {reffasta}")
print(f"Sample: {sample}")

pattern = re.compile(r'mapping/(?P<sample>[A-Za-z0-9_-]+)/(?P<flowcell>.+).bam')
abam_dict = {}

for infile in Path(f"mapping/{sample}/").glob('*markdup*.bam'):
    match = pattern.search(str(infile))
    if match:
        abam_dict[match.group('flowcell')] = str(infile)

movies = list(abam_dict.keys())
abams = list(abam_dict.values())

print("movies available {m}".format(m=movies))
print("paths: {a}".format(a=abams))
    
include: 'rules/variants.smk'

#rule all:
#  input:
#    expand(f"sbb_data/{sample}/{flowcell}.fastq.gz",sample=sample)
