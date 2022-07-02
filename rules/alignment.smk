if mapper == "bwa":
    print("Selected BWA")
    ruleorder: bwa > dragmap
else:
    print("Selected dragmap")
    ruleorder: dragmap > bwa


rule cutadapt:
    input: "sbb_data/{sample}/{flowcell}.fastq.gz"
    output: "cutadapt/{sample}/{flowcell}_trimmed.fastq.gz"
    benchmark: "benchmarks/{sample}/cutadapt_{flowcell}.tsv"
    log: "logs/{sample}/cutadapt_{flowcell}.log"
    conda: "envs/cutadapt.yaml"
    threads: 24
    message: "Running cutadapt with {threads} threads on {input}."
    shell:
        """
        cutadapt --cores {threads} -a {adapter} --overlap 8 -o {output} {input} > {log} 2>&1
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
        dragen-os --RGSM {wildcards.sample} \
                  --num-threads {threads} \
                  -r {params.refdict} \
                  -1 {input} \
                  --output-file-flowcell \
                  {output} \
                  2> {log} | \
                  samtools view -bS | \
                  samtools sort -o {output}
        """

rule bwa:
    input: "cutadapt/{sample}/{flowcell}_trimmed.fastq.gz"
    output: "mapping/{sample}/{flowcell}.bam"
    params:
        readgroup = "@RG\\tID:{flowcell}\\tSM:{sample}\\tLB:{flowcell}_{sample}\\tPL:sbb"
    benchmark: "benchmarks/{sample}/dragen_{flowcell}.tsv"
    log: "logs/{sample}/bwa_{flowcell}.log"
    conda: "envs/bwa.yaml"
    threads: 24
    message: "Running bwa-mem with {threads} threads on {input}."
    shell:
        """
        bwa mem -R "{params.readgroup}" \
                -t {threads} \
                {reffasta} \
                {input} | \
                samtools sort -@{threads} \
                              -o {output} -   
        """

rule markduplicates:
    input: "mapping/{sample}/{flowcell}.bam"
    output: bam = "mapping/{sample}/{flowcell}_markdup.bam", bai = "mapping/{sample}/{flowcell}_markdup.bam.bai"
    params:
        metrics_out = "mapping/{sample}/marked_dup_metrics_{flowcell}.txt"
    benchmark: "benchmarks/{sample}/mark_dups_{flowcell}.tsv"
    log: "logs/{sample}/markdups_{flowcell}.log"
    conda: "envs/gatk.yaml"
    threads: 24
    message: "Running gatk MarkDuplicatesSpark with {threads} threads on {input}."
    shell:
        """
        gatk MarkDuplicatesSpark \
                    --conf "spark.executor.cores={threads}" \
                    --tmp-dir {temp} \
                    --input {input} \
                    --output {output.bam} \
                    -M {params.metrics_out}
        """

rule mosdepth:
    input: "mapping/{sample}/{flowcell}_markdup.bam"
    output:
        "mapping/{sample}/mosdepth/{flowcell}_markdup.mosdepth.global.dist.txt",
        "mapping/{sample}/mosdepth/{flowcell}_markdup.mosdepth.region.dist.txt",
        "mapping/{sample}/mosdepth/{flowcell}_markdup.mosdepth.summary.txt",
        "mapping/{sample}/mosdepth/{flowcell}_markdup.regions.bed.gz"
    log: "logs/{sample}/mosdepth/{flowcell}.log"
    benchmark: "benchmarks/{sample}/mosdepth/{flowcell}.tsv"
    params:
        by = "500",
        flowcell = "mapping/{sample}/mosdepth/{flowcell}_markdup",
        extra = "--no-per-base --use-median"
    threads: 4
    conda: "envs/mosdepth.yaml"
    message: "Executing {rule}: Calculating coverage of {input} using mosdepth."
    shell:
        """
        (mosdepth --threads {threads} --by {params.by} \
            {params.extra} {params.flowcell} {input}) > {log} 2>&1
        """

rule seqtkfqchk:
    input: "cutadapt/{sample}/{flowcell}_trimmed.fastq.gz"
    output: "cutadapt/{sample}/seqtk/{flowcell}_seqtkfqchk.txt"
    log: "logs/{sample}/seqtk/{flowcell}.log"
    benchmark: "benchmarks/{sample}/seqtk/{flowcell}.tsv"
    params:
        qfilter = 20  # 20 is program default
    threads: 4
    conda: "envs/seqtk.yaml"
    message: "Executing {rule}: Calculating basic read stats on {input} using seqtk."
    shell:
        """
        seqtk fqchk -q {params.qfilter} {input} >{output} 2> {log}
        """

rule seqtkfqchkplot:
    input: "cutadapt/{sample}/seqtk/{flowcell}_seqtkfqchk.txt"
    output: "cutadapt/{sample}/seqtk/{flowcell}_seqtkfqchk.txt.png"
    log: "logs/{sample}/seqtk/{flowcell}_plot.log"
    params:
        qfilter = 20  # 20 is program default
    threads: 4
    conda: "envs/pyenv.yaml"
    message: "Executing {rule}: Calculating basic read stats on {input} using seqtk."
    shell:
        """
        python scripts/plots.py {input} 2> {log}      
        """
