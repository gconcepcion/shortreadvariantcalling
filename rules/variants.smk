rule haplotypecaller:
    input:
        bams = abams
    output: "gatk/{sample}/{sample}.raw.vcf.gz"
    params:
        ref = reffasta,
        inbams = ' '.join([f"-I {i}" for i in abams])
    benchmark: "benchmarks/{sample}/gatk_{sample}.tsv"
    log: "logs/{sample}/gatk_{sample}.log"
    conda: "envs/gatk.yaml"
    threads: 24
    message: "Running HaplotypeCaller with {threads} threads on {input}."
    shell:
        """
        gatk --java-options "-Xmx4g" HaplotypeCaller \
          --native-pair-hmm-threads {threads} \
          -R {params.ref} \
          {params.inbams} \
          -O {output}
        """

rule aggregatecoverage:
    ### count total bases from all aligned bam files (mosdepth output)
    input: expand(f"mapping/{sample}/mosdepth/{{flowcell}}.mosdepth.summary.txt", flowcell=flowcells)
    output: "gatk/{sample}/{sample}.depthofcov.txt"
    benchmark: "benchmarks/{sample}/{sample}_gatk_depthofcoverage.tsv"
    log: "logs/{sample}/{sample}_gatk_depthofcoverage.log"
    threads: 4
    message: "Running DepthOfCoverage with {threads} threads on {input}."
    run:
        totals = []
        with open(f"{output}", "w+") as o:
            for bam in input:
                with open(bam, 'r') as i:
                    contents = i.read().splitlines()
                flowcell = Path(bam).stem 
                for y in contents:
                    if y.startswith("total_"):
                        o.write(flowcell+"\t"+y+"\n")
                        totals.append(int(y.split()[2]))
            totalbp = sum(totals)
            o.write(f"total_bp\t{totalbp}")


rule selectsnps:
    input: vcf = "gatk/{sample}/{sample}.raw.vcf.gz",
           ref = f"{reffasta}"
    output: "gatk/{sample}/{sample}.SNPs.vcf.gz"
    params:
        filter = f"QUAL < {qual}"
    conda: "envs/gatk.yaml"
    log: "logs/{sample}/gatk_extractsnps_{sample}.log"
    message: "Extracting SNPs from {input}"
    shell:
        """
        gatk SelectVariants \
            -R {input.ref}  \
            -V {input.vcf} \
            -select-type SNP \
            -O {output} \
            2> {log}
        """

rule selectindels:
    input: vcf = "gatk/{sample}/{sample}.raw.vcf.gz",
           ref = f"{reffasta}"
    output: "gatk/{sample}/{sample}.INDELs.vcf.gz"
    params:
        filter = f"QUAL < {qual}"
    conda: "envs/gatk.yaml"
    log: "logs/{sample}/gatk_extractsnps_{sample}.log"
    message: "Extracting INDELs from {input}"
    shell:
        """
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            -select-type INDEL \
            -O {output} \
            2> {log}
        """

rule variantfilterINDELs:
    input: "gatk/{sample}/{sample}.INDELs.vcf.gz"
    output: "gatk/{sample}/{sample}.INDELs.filtered.vcf.gz"
    params:
        qual = f"QUAL < {qual}"
    log: "logs/{sample}/gatk_INDELs_{sample}.log"
    conda: "envs/gatk.yaml"
    message: "Running gatk VariantFiltration with {threads} threads on {input}."
    shell:
        """
        gatk VariantFiltration \
            -V {input} \
            --filter-expression "{params.qual}" \
            --filter-name "QUALFilter" \
            -O {output} \
            2> {log}
        """

rule variantfilterSNPs:
    input: "gatk/{sample}/{sample}.SNPs.vcf.gz"
    output: "gatk/{sample}/{sample}.SNPs.filtered.vcf.gz"
    params:
        qual = f"QUAL < {qual}"
    log: "logs/{sample}/gatk_SNPs_{sample}.log"
    conda: "envs/gatk.yaml"
    message: "Running gatk VariantFiltration with {threads} threads on {input}."
    shell:
        """
        gatk VariantFiltration \
            -V {input} \
            --filter-expression "{params.qual}" \
            --filter-name "QUALFilter" \
            -O {output} \
            2> {log}
        """

rule updatereference:
    input: vcf = "gatk/{sample}/{sample}.INDELs.vcf.gz",
           ref = f"{reffasta}"
    output: "gatk/{sample}/{refbase}_updated.fasta"
    conda: "envs/gatk.yaml"
    log: "logs/{sample}/gatk_{refbase}_updatedref_{sample}.log"
    message: "Updating reference with new variants from {input}."
    shell:
        """
        gatk FastaAlternateReferenceMaker \
           -R {input.ref} \
           -O {output} \
           -V {input.vcf} \
           2> {log}
        """
