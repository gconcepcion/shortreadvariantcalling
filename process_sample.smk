import re
from pathlib import Path

configfile: "reference.yaml"         # reference configuration
configfile: "config.yaml"            # config options

reffasta = config['ref']['fasta']
sample = config['sample']
print("Reference: {r}".format(r=reffasta))
print("Sample: {s}".format(s=sample))


#ipattern = re.compile(r'samples/(?P<sample>[A-Za-z0-9_-]+)/aligned/(?P<movie>m\d{5}U?e?_\d{6}_\d{6})\.(?P<reference>.*).bam')
pattern = re.compile(r"dragen/{s}/(?P<run>.*).bam".format(s=sample))
abam_dict = {}

for infile in Path("dragen/{s}/".format(s=sample)).glob('*markdup*.bam'):
    match = pattern.search(str(infile))
    if match:
        abam_dict[match.group('run')] = str(infile)

movies = list(abam_dict.keys())
abams = list(abam_dict.values())

print("movies available {m}".format(m=movies))
print("paths: {a}".format(a=abams))
    
include: 'rules/variants.smk'

#rule all:
#  input:
#    expand(f"sbb_data/{sample}/{flowcell}.fastq.gz",sample=sample)
