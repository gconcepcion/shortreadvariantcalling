rule haplotypecaller:
  input:
    bams = abams
  output: "gatk/{sample}/{flowcell}.raw.vcf.gz"
  params:
    ref = reffasta,
    inbams = ' '.join(["-I {i}".format(i=i) for i in abams])
  benchmark: "benchmarks/{sample}/gatk_{flowcell}.tsv"
  log: "logs/{sample}/gatk_{flowcell}.log"
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

rule variantfilter:
  input: "gatk/{sample}/{flowcell}.raw.vcf.gz"
  output: "gatk/{sample}/{flowcell}.filtered.vcf.gz"
  params:
    qual = "QUAL < 10.4139"
  log: "logs/{sample}/gatk_{flowcell}.log"
  conda: "envs/gatk.yaml"
  message: "Running gatk VariantFiltration with {threads} threads on {input}."
  shell:
    """
    gatk VariantFiltration \
        -V {input} \
        --filter-expression "{params.qual}" \
        --filter-name "DRAGENHardQUAL" \
        -O {output}
    """
