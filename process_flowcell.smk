
configfile: "reference.yaml"         # reference configuration
configfile: "config.yaml"         # reference configuration

temp = '/scratch'
mapper = config['mapper']
sample = config['sample']
adapter = config['adapter']
datadir = config['datadir']

flowcell = config['flowcell']
id = 1
ref = config['ref']['shortname']
refdict = config['ref']['dict']
reffasta = config['ref']['fasta']


# align reads with dragmap and markduplicates with gatk
include: 'rules/alignment.smk'

print(f"flowcell: {flowcell}")
print(f"For Sample: {sample}")

targets = []
fastqinput = os.path.join(f"{datadir}/{sample}/{flowcell}.fastq.gz")

targets.append(fastqinput)

# align reads with dragmap or bwa
include: 'rules/alignment.smk'
# trimmed raw BAMS
if 'cutadapt' in config['flowcell_targets']:
    targets.extend(
        [f"cutadapt/{sample}/{flowcell}_trimmed.{suffix}" for suffix in ['fastq.gz']])

if 'seqtk' in config['flowcell_targets']:
    targets.extend([f"cutadapt/{sample}/seqtk/{flowcell}_seqtkfqchk.txt"])
    targets.extend([f"cutadapt/{sample}/seqtk/{flowcell}_seqtkfqchk.txt.png"])

if 'alignment' in config['flowcell_targets']:
    targets.extend([f"mapping/{sample}/{flowcell}.bam"])

if 'markdups' in config['flowcell_targets']:
    targets.extend(
        [f"mapping/{sample}/{flowcell}_markdup.{suffix}" for suffix in ['bam', 'bam.bai']])

if 'mosdepth' in config['flowcell_targets']:
    targets.extend([f"mapping/{sample}/mosdepth/{flowcell}_markdup.{suffix}" for suffix in ['mosdepth.global.dist.txt',
                                                                                    'mosdepth.region.dist.txt',
                                                                                    'mosdepth.summary.txt',
                                                                                    'regions.bed.gz']])

rule all:
    input: targets
