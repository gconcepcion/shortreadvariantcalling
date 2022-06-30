if mapper == "bwa":
    print("Selected BWA")
    ruleorder: bwa > dragmap
else:
    print("Selected dragmap")
    ruleorder: dragmap > bwa


rule cutadapt:
  input: "sbb_data/{sample}/{flowcell}.fastq.gz"
  output: "cutadapt/{sample}/{flowcell}_trimmed.fastq.gz"
  params:
    adapter = "ATCGATTCGTGCTCGATGAACCGGGCGCTTA"
  benchmark: "benchmarks/{sample}/cutadapt_{flowcell}.tsv"
  log: "logs/{sample}/cutadapt_{flowcell}.log"
  conda: "envs/cutadapt.yaml"
  threads: 24
  message: "Running cutadapt with {threads} threads on {input}."
  shell:
    """
    cutadapt --cores {threads} -a {params.adapter} --overlap 8 -o {output} {input}
    """

rule dragmap:
  input: "cutadapt/{sample}/{flowcell}_trimmed.fastq.gz"
  output: "mapping/{sample}/{flowcell}.bam"
  params:
    refdict = "reference/"
  benchmark: "benchmarks/{sample}/dragen_{flowcell}.tsv"
  log: "logs/{sample}/dragen_{flowcell}.log"
  conda: "envs/dragen.yaml"
  threads: 24
  message: "Running dragmap with {threads} threads on {input}."
  shell:
    """
    dragen-os --RGSM {wildcards.sample} --num-threads {threads} -r {params.refdict} -1 {input} --output-file-prefix {output} 2> {log} | 
    samtools view -bS | 
    samtools sort -o {output}
    """

rule bwa:
   input: "sbb_data/{sample}/{flowcell}_trimmed.fastq.gz"
   output: "mapping/{sample}/{flowcell}.bam"
   params:
     readgroup="@RG\\tID:{id}\\tSM:{sample}\\tLB:{id}_{sample}\\tPL:sbb".format(id=id, sample=sample)   
   benchmark: "benchmarks/{sample}/dragen_{flowcell}.tsv"
   log: "logs/{sample}/bwa_{flowcell}.log"
   conda: "envs/bwa.yaml"
   threads: 24
   message: "Running bwa-mem with {threads} threads on {input}."
   shell:
     """
     bwa mem -R "{params.readgroup}" -t {threads} {reffasta} {input} | samtools sort -@{threads} -o {output} -   
     """

rule markduplicates:
   input: "mapping/{sample}/{flowcell}.bam"
   output: "mapping/{sample}/{flowcell}_markdup.bam"
   params:
     metrics_out = "mapping/{sample}/marked_dup_metrics_{flowcell}.txt"
   benchmark: "benchmarks/{sample}/mark_dups_{flowcell}.tsv"
   log: "logs/{sample}/markdups_{flowcell}.log"
   conda: "envs/gatk.yaml"
   threads: 24
   message: "Running gatk MarkDuplicatesSpark with {threads} threads on {input}."
   shell:
     """
     gatk MarkDuplicatesSpark --conf "spark.executor.cores={threads}" --tmp-dir {temp} --input {input} --output {output} -M {params.metrics_out}
     """
   
