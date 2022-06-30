
configfile: "reference.yaml"         # reference configuration
configfile: "config.yaml"         # reference configuration

temp = '/scratch'
mapper = config['mapper']
sample = config['sample']
id=1
#readgroup="@RG\\tID:{id}\\tSM:{sample}\\tLB:{id}_{sample}\\tPL:sbb".format(id=id, sample=sample)
ref = config['ref']['shortname']
refdict = config['ref']['dict']
reffasta = config['ref']['fasta']


#align reads with dragmap and markduplicates with gatk
include: 'rules/alignment.smk'

