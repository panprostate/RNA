# rule VC_create_dictIndex:
#     input:
#         reference=config["reference"]
#     output:
#         refDict=config["dict"]
#         fai="resources/genome.fai"
#     threads: 2
#     resources:
#         mem_mb=4000,
#         runtime_min="00:30:00"
#     conda:
#         "../envs/variantCalling.yaml"
#     shell:
#         """
#         gatk CreateSequenceDictionary -R {input.reference}
#         samtools faidx {input.reference}
#         """

rule VC_create_intervalList:
    input:
        gtf=config["gtf"],
        refDict=config["dict"]
    output:
        tmp=temp("results/variantCalling/intervalList/exome.bed"),
        intervals="results/variantCalling/intervalList/exons.interval_list"
    threads: 2
    resources:
        mem_mb=4000,
        runtime_min="00:30:00"
    script:
        "scripts/createInterval.R"



rule VC_markDuplicates:
    input:
        bam="results/sortedBams/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        dedup=temp("results/variantCalling/markDuplicates/{sample}.bam"),
        index=temp("results/variantCalling/markDuplicates/{sample}.bam.bai"),
        metrics="results/variantCalling/metrics/{sample}_dedup.metrics"
    params:
        compression=config["compression_level"]
    threads: 2
    resources:
        mem_mb=config["mem_markDuplicates"],
        runtime_min=config["rt_markDuplicates"]
    benchmark:
        "benchmark/markDuplicates/{sample}.tsv"
    log:
        "logs/markDuplicates/{sample}.log"
    conda:
        "../envs/variantCalling.yaml"
    shell:
        """
        gatk \
        MarkDuplicates \
        --INPUT {input.bam} \
        --OUTPUT {input.dedup} \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT \
        --METRICS_FILE output.metrics \
        --COMPRESSION_LEVEL {params.compression} 2>> {log}
        """

rule VC_splitNCigars:
    input:
        reference=config["reference"],
        refDict=config["dict"],
        refIndex=config["faidx"],
        bam="results/variantCalling/markDuplicates/{sample}.bam"
    output:
        sacr=temp("results/variantCalling/splitNCigars/{sample}.bam"),
        index=temp("results/variantCalling/splitNCigars/{sample}.bam.bai")
    params:
        compression=config["compression_level"]
    threads: 2
    resources:
        mem_mb=config["mem_splitNCigars"],
        runtime_min=config["rt_splitNCigars"]
    benchmark:
        "benchmark/splitNCigars/{sample}.tsv"
    log:
        "logs/splitNCigars/{sample}.log"
    conda:
        "../envs/variantCalling.yaml"
    shell:
        """
        gatk \
        SplitNCigarReads \
        -R {input.reference} \
        -I {input.bam} \
        -O {output.sacr} \
        --COMPRESSION_LEVEL {params.compression} 2>> {log}
        """

rule VC_baseRecalibrator:
    input:
        reference=config["reference"],
        refDict=config["dict"],
        refIndex=config["faidx"],
        bam="results/variantCalling/splitNCigars/{sample}.bam",
        bai="results/variantCalling/splitNCigars/{sample}.bam.bai",
        dbSNP=config["dbSNP"],
        knowIndels=config["knowIndels"]
    output:
        table="results/variantCalling/recalibration/{sample}.tbl"
    threads: 2
    resources:
        mem_mb=config["mem_baseRecalibrator"],
        runtime_min=config["rt_baseRecalibrator"]
    benchmark:
        "benchmark/baseRecalibrator/{sample}.tsv"
    log:
        "logs/baseRecalibrator/{sample}.log"
    conda:
        "../envs/variantCalling.yaml"
    shell:
        """
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
        -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
        -Xloggc:gc_log.log -Xms4000m" \
        BaseRecalibrator \
        -R {input.reference} \
        -I {input.bam} \
        --use-original-qualities \
        -O {output.table} \
        -known-sites {input.dbSNP} \
        -known-sites {input.knowIndels} 2>> {log}
        """

rule VC_applyBQSR:
    input:
        reference=config["reference"],
        refDict=config["dict"],
        refIndex=config["faidx"],
        bam="results/variantCalling/splitNCigars/{sample}.bam",
        bai="results/variantCalling/splitNCigars/{sample}.bam.bai",
        table="results/variantCalling/recalibration/{sample}.tbl"
    output:
        recalbam="results/variantCalling/recalibration/{sample}.bam",
        bai="results/variantCalling/recalibration/{sample}.bam.bai",
        report="results/variantCalling/metrics/{sample}_recalibration.metrics"
    params:
        compression=config["compression_level"]
    threads: 2
    resources:
        mem_mb=config["mem_applyBQSR"],
        runtime_min=config["rt_applyBQSR"]
    benchmark:
        "benchmark/applyBQSR/{sample}.tsv"
    log:
        "logs/applyBQSR/{sample}.log"
    conda:
        "../envs/variantCalling.yaml"
    shell:
        """
        gatk \
        --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
        -XX:+PrintGCDetails -Xloggc:gc_log.log \
        -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m" \
        ApplyBQSR \
        --add-output-sam-program-record \
        -R {input.reference} \
        -I {input.recalbam} \
        --use-original-qualities \
        -O {output.bam} \
        --bqsr-recal-file {output.report} \
        --COMPRESSION_LEVEL {params.compression} 2>> {log}
        """

rule VC_haplotypeCaller:
    input:
        reference=config["reference"],
        refDict=config["dict"],
        refIndex=config["faidx"],
        bam="results/variantCalling/recalibration/{sample}.bam",
        bai="results/variantCalling/recalibration/{sample}.bam.bai",
        dbSNP=config["dbSNP"],
        intervals="results/variantCalling/intervalList/exons.interval_list"
    output:
        vcf=temp("results/variantCalling/vcf/{sample}.vcf.gz")
    threads: 2
    resources:
        mem_mb=config["mem_haplotypeCaller"],
        runtime_min=config["rt_haplotypeCaller"]
    benchmark:
        "benchmark/haplotypeCaller/{sample}.tsv"
    log:
        "logs/haplotypeCaller/{sample}.log"
    conda:
        "../envs/variantCalling.yaml"
    shell:
        """
        gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        HaplotypeCaller \
        -R {input.reference} \
        -I {input.bam} \
        -L {input.intervals} \
        -O {output.vcf} \
        -dont-use-soft-clipped-bases \
        --standard-min-confidence-threshold-for-calling 20 \
        --dbsnp {dbSNP} 2>> {log}
        """

rule VC_filterVCF:
    input:
        reference=config["reference"],
        refDict=config["dict"],
        refIndex=config["faidx"],
        vcf="results/variantCalling/vcf/{sample}.vcf.gz"
    output:
        vcf_f="results/variantCalling/vcf_filtered/{sample}.vcf.gz",
        vcfIdx="results/variantCalling/vcf_filtered/{sample}.vcf.gz.tbi"
    threads: 2
    resources:
        mem_mb=config["mem_haplotypeCaller"],
        runtime_min=config["rt_haplotypeCaller"]
    benchmark:
        "benchmark/filterVCF/{sample}.tsv"
    log:
        "logs/filterVCF/{sample}.log"
    conda:
        "../envs/variantCalling.yaml"
    shell:
        """
        gatk \
        VariantFiltration \
        --R {input.reference} \
        --V {input.vcf} \
        --window 35 \
        --cluster 3 \
        --filter-name "FS" \
        --filter "FS > 30.0" \
        --filter-name "QD" \
        --filter "QD < 2.0" \
        -O {output.vcf_f}
        """    
