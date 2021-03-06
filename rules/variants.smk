rule haplotypecaller:
    input:
        bams = abams
    output: "sample/{sample}/gatk/{sample}.raw.vcf.gz"
    params:
        ref = reffasta,
        inbams = ' '.join([f"-I {i}" for i in abams]),
        minmapq = {minmapq},
        ploidy = config['ploidy'],
        intervals = intervals
    benchmark: "benchmarks/{sample}/gatk_{sample}.tsv"
    log: "logs/{sample}/gatk_{sample}.log"
    conda: "envs/gatk.yaml"
    threads: 24
    message: "Running HaplotypeCaller with {threads} threads on {input}."
    shell:
        """
        gatk --java-options "-Xmx4g" HaplotypeCaller \
          --native-pair-hmm-threads {threads} \
          --minimum-mapping-quality {params.minmapq} \
          --ploidy {params.ploidy} \
          --intervals {params.intervals} \
          -R {params.ref} \
          {params.inbams} \
          -O {output} \
          2> {log}
        """

rule aggregatecoverage:
    ### count total bases from all aligned bam files (mosdepth output)
    input: expand(f"sample/{sample}/mapping/mosdepth/{{flowcell}}.mosdepth.summary.txt", flowcell=flowcells)
    output: "sample/{sample}/gatk/{sample}.depthofcov.txt"
    benchmark: "benchmarks/{sample}/{sample}_gatk_depthofcoverage.tsv"
    log: "logs/{sample}/{sample}_gatk_depthofcoverage.log"
    threads: 4
    message: "Aggregating flowcell coverage"
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
    input: vcf = "sample/{sample}/gatk/{sample}.raw.vcf.gz",
           ref = f"{reffasta}"
    output: "sample/{sample}/gatk/{sample}.SNPs.vcf.gz"
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
    input: vcf = "sample/{sample}/gatk/{sample}.raw.vcf.gz",
           ref = f"{reffasta}"
    output: "sample/{sample}/gatk/{sample}.INDELs.vcf.gz"
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

rule updatereference:
    input: vcf = "sample/{sample}/gatk/{sample}.INDELs.vcf.gz",
           ref = f"{reffasta}"
    output: "sample/{sample}/gatk/{refbase}_updated.fasta"
    conda: "envs/gatk.yaml"
    log: "logs/{sample}/gatk_{refbase}_updatedref_{sample}.log"
    message: "Updating reference with INDEL variants from {input}."
    shell:
        """
        gatk FastaAlternateReferenceMaker \
           -R {input.ref} \
           -O {output} \
           -V {input.vcf} \
           2> {log}
        """

rule updatereferenceSNPs:
    input: vcf = "sample/{sample}/gatk/{sample}.SNPs.vcf.gz",
           ref = f"{reffasta}"
    output: "sample/{sample}/gatk/{refbase}_SNPsupdated.fasta"
    conda: "envs/gatk.yaml"
    log: "logs/{sample}/gatk_{refbase}_SNPsupdatedref_{sample}.log"
    message: "Updating reference with new SNP variants from {input}."
    shell:
        """
        gatk FastaAlternateReferenceMaker \
           -R {input.ref} \
           -O {output} \
           -V {input.vcf} \
           2> {log}
        """

rule updatereferenceALL:
    input: vcf = "sample/{sample}/gatk/{sample}.raw.vcf.gz",
           ref = f"{reffasta}"
    output: "sample/{sample}/gatk/{refbase}_ALLupdated.fasta"
    conda: "envs/gatk.yaml"
    log: "logs/{sample}/gatk_{refbase}_ALLupdatedref_{sample}.log"
    message: "Updating reference with all variants from {input}."
    shell:
        """
        gatk FastaAlternateReferenceMaker \
           -R {input.ref} \
           -O {output} \
           -V {input.vcf} \
           2> {log}
        """

rule selectHQindels:
    input: vcf = "sample/{sample}/gatk/{sample}.INDELs.vcf.gz"
    output: "sample/{sample}/gatk/{sample}.HQ.INDELs.vcf.gz"
    params:
        filter = f"QD<{qdfilter}"
    conda: "envs/gatk.yaml"
    log: "logs/{sample}/gatk_extractsnps_{sample}.log"
    message: "Extracting INDELs from {input}"
    shell:
        """
        gatk VariantFiltration \
            -V {input.vcf} \
            --filterExpression {params.filter} \
            -O {output} \
            2> {log}
        """

rule sampleaccuracy:
    input: 
        rawvcf = 'sample/{sample}/gatk/{sample}.raw.vcf.gz',
        snpvcf = 'sample/{sample}/gatk/{sample}.SNPs.vcf.gz',
        indelvcf = 'sample/{sample}/gatk/{sample}.INDELs.vcf.gz'
    output: "sample/{sample}/gatk/{sample}.acc"
    params: 
        reflen = reflen
    log: "logs/{sample}/{sample}_accuracy.log"
    message: "Calculate overall accuracy (QV) from {input}"
    run:
        import gzip
        import math

        with open(f"{output}", 'w+') as o:
            for vcf in [input.snpvcf, input.indelvcf, input.rawvcf]:

                with gzip.open(vcf, 'rb') as i:
                    contents = i.read().splitlines()
                variantdata = []
                for line in contents:
                    if not line.startswith(b'#'):
                        variantdata.append(line)   
                errors = len(variantdata)
                pcalled = errors / int(params.reflen)
                qv = -10 * math.log10(pcalled)
                o.write(f"{vcf}\t{qv}\t{pcalled}\n")   
