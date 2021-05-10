rule trim:
    input:
        get_fastq
    output:
        fastq1=temp("results/trim/{sample}/{sample}_{unit}_R1_trimmed.fastq.gz"),
        fastq2=temp("results/trim/{sample}/{sample}_{unit}_R2_trimmed.fastq.gz"),
        qc="results/qc/{sample}/{sample}_{unit}_trim.log"
    params:
        adapters=config["adapters"],
        compression=config["compression_level"]
    threads: config["ncpus_trim"]
    resources:
        mem_mb=config["mem_trim"],
        runtime_min=config["rt_trim"]
    benchmark:
        "benchmark/trimqc/{sample}_{unit}.tsv"
    log:
        "logs/trimqc/{sample}_{unit}.log"
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "cutadapt"
        " {params.adapters}"
        " -m 31"
        " --compression-level {params.compression}"
        " -o {output.fastq1}"
        " -p {output.fastq2}"
        " -j {threads}"
        " {input} > {output.qc} 2>>{log}"

rule star_index:
    input:
        reference=config["reference"],
        gtf=config["gtf"]
    output:
        index="results/index/hg19/SA"
    params:
        idx="results/index/hg19/"
    threads: config["ncpus_STAR_index"]
    resources:
        mem_mb=config["mem_STAR_index"],
        runtime_min=config["rt_STAR_index"]
    benchmark:
        "benchmark/STARindex/STARindex.tsv"
    log:
        "logs/STARindex/STARindex.log"
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "STAR"
        " --runMode genomeGenerate"
        " --genomeDir {params.idx}"
        " --genomeFastaFiles {input.reference}"
        " --sjdbOverhang 250"
        " --sjdbGTFfile {input.gtf}"
        " --runThreadN {threads} 2>&1 | tee -a {log}"

rule star_1pass:
    input:
        index="results/index/hg19/SA",
        f1="results/trim/{sample}/{sample}_{unit}_R1_trimmed.fastq.gz",
        f2="results/trim/{sample}/{sample}_{unit}_R2_trimmed.fastq.gz"
    output:
        SJ="results/STAR_1p/{sample}_{unit}.SJ.out.tab"
    params:
        idx="results/index/hg19"
    threads: config["ncpus_STAR_align"]
    resources:
        mem_mb=config["mem_STAR_align"],
        runtime_min=config["rt_STAR_align"]
    benchmark:
        "benchmark/STAR_1p/{sample}_{unit}.tsv"
    log:
        "logs/STAR_1p/{sample}_{unit}.log"
    conda:
        "../envs/RNAseq.yaml"
    shell:
        """
        STAR \
        --genomeDir {params.idx} \
        --readFilesIn {input.f1} {input.f2} \
        --runThreadN 8 \
        --readFilesCommand zcat \
        --outFilterMultimapScoreRange 1 \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 10 \
        --alignIntronMax 500000 \
        --alignMatesGapMax 1000000 \
        --sjdbScore 2 \
        --alignSJDBoverhangMin 5 \
        --genomeLoad NoSharedMemory \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterScoreMinOverLread 0.33 \
        --sjdbOverhang 250 \
        --outSAMstrandField intronMotif \
        --peOverlapNbasesMin 10 \
        --alignSplicedMateMapLminOverLmate 0.5 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimSegmentMin 10 \
        --chimOutType Junctions WithinBAM HardClip \
        --chimJunctionOverhangMin 10 \
        --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 \
        --chimScoreSeparation 1 \
        --chimSegmentReadGapMax 3 \
        --chimMultimapNmax 20 \
        --outSAMtype None \
        --outSAMmode None \
        --outFileNamePrefix STAR_1p/{wildcards.sample} 2>&1 | tee -a {log}
        """

rule star_2pass:
    input:
        f1="results/trim/{sample}/{sample}_{unit}_R1_trimmed.fastq.gz",
        f2="results/trim/{sample}/{sample}_{unit}_R2_trimmed.fastq.gz",
        index="results/index/hg19/SA",
        SJ=gather_SJ
    output:
        bam=temp("results/STAR_2p/{sample}_{unit}.Aligned.out.bam"),
        star_logs=multiext("STAR_2p/{sample}_{unit}", ".Log.final.out", ".Log.out", ".Chimeric.out.junction")
    params:
        compression=config["compression_level"],
        idx="results/index/hg19"
    threads: config["ncpus_STAR_align"]
    resources:
        mem_mb=config["mem_STAR_align"],
        runtime_min=config["rt_STAR_align"]
    benchmark:
        "benchmark/STAR_2p/{sample}_{unit}.tsv"
    log:
        "logs/STAR_2p/{sample}_{unit}.log"
    conda:
        "../envs/RNAseq.yaml"
    shell:
        """
        STAR \
        --genomeDir {params.idx} \
        --readFilesIn {input.f1} {input.f2} \
        --runThreadN 8 \
        --readFilesCommand zcat \
        --outBAMCompression {params.compression} \
        --outFilterMultimapScoreRange 1 \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 10 \
        --alignIntronMax 500000 \
        --alignMatesGapMax 1000000 \
        --sjdbScore 2 \
        --alignSJDBoverhangMin 5 \
        --genomeLoad NoSharedMemory \
        --sjdbFileChrStartEnd {input.SJ} \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterScoreMinOverLread 0.33 \
        --sjdbOverhang 250 \
        --outSAMstrandField intronMotif \
        --peOverlapNbasesMin 10 \
        --alignSplicedMateMapLminOverLmate 0.5 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimSegmentMin 10 \
        --chimOutType Junctions WithinBAM HardClip \
        --chimJunctionOverhangMin 10 \
        --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 \
        --chimScoreSeparation 1 \
        --chimSegmentReadGapMax 3 \
        --chimMultimapNmax 20 \
        --outSAMtype BAM Unsorted \
        --outSAMmode None \
        --outWigType wiggle \
        --outWigStrand Unstranded \
        --outWigNorm None \
        --outFileNamePrefix STAR_2p/{wildcards.sample} 2>&1 | tee -a {log}
        """


rule sortBam:
    input:
        bam="results/STAR_2p/{sample}_{unit}.Aligned.out.bam"
    output:
        sbam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
        bai="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam.bai"
    threads: config["ncpus_sortBam"]
    resources:
        mem_mb=config["mem_sortBam"],
        runtime_min=config["rt_sortBam"]
    benchmark:
        "benchmark/sortBam/{sample}_{unit}.tsv"
    log:
        "logs/sortBam/{sample}_{unit}.log"
    conda:
        "../envs/RNAseq.yaml"
    shell:
        """
        samtools sort -@ {threads} -m 2G {input.bam} -o {output.sbam} 2>>{log}
        samtools index {output.sbam}
        """

rule mergeBam:
    input:
        bams=gather_bams
    output:
        mbam="results/mergedBam/{sample}.Aligned.sortedByCoord.out.bam",
        bai="results/mergedBam/{sample}.Aligned.sortedByCoord.out.bam.bai"
    params:
        bams=gather_bams
    threads: config["ncpus_mergeBam"]
    resources:
        mem_mb=config["mem_mergeBam"],
        runtime_min=config["rt_mergeBam"]
    benchmark:
        "benchmark/mergedBam/{sample}.tsv"
    log:
        "logs/mergedBam/{sample}.log"
    conda:
        "../envs/RNAseq.yaml"
    shell:
        """
        INPUT=({input.bams})
        if ((${{#INPUT[@]}} == 1)); then
            cp {input.bams} {output.mbam}
            cp {input.bams}.bai {output.mbam}.bai
        else
            samtools merge -f -1 {output.mbam} {input.bams} 2>>{log}
            samtools index {output.mbam}
        fi
        """


rule salmon_index:
    input:
        transcripts=config["transcripts"],
        reference=config["reference"]
    output:
        "results/index/salmon_hg19/ctable.bin"
    params:
        idx="results/index/salmon_hg19/"
    threads: config["ncpus_salmonIndex"]
    resources:
        mem_mb=config["mem_salmonIndex"],
        runtime_min=config["rt_salmonIndex"]
    benchmark:
        "benchmark/salmonIndex/salmonIndex.tsv"
    log:
        "logs/salmonIndex/salmonIndex.log"
    conda:
        "../envs/RNAseq.yaml"
    shell:
        """
        cat {input.transcripts} {input.reference} > data/gentrome.fa.gz
        grep "^>" <(gunzip -c {input.reference}) | cut -d " " -f 1 > data/decoys.txt
        salmon index -t gentrome.fa.gz -i {params.idx} --decoys data/decoys.txt -k 31 2>&1 | tee -a {log}
        rm -f data/gentrome.fa.gz data/decoys.txt
        """

rule salmon_quant:
    input:
        idx="results/index/salmon_hg19/ctable.bin",
        f1=gather_salmon_input1,
        f2=gather_salmon_input2
    output:
        tc="results/salmon/{sample}/quant.sf",
    params:
        idx="results/index/salmon_hg19/",
        dir="results/salmon/{sample}/"
    threads: config["ncpus_salmonQuant"]
    resources:
        mem_mb=config["mem_salmonQuant"],
        runtime_min=config["rt_salmonQuant"]
    benchmark:
        "benchmark/salmonQuant/{sample}.tsv"
    log:
        "logs/salmonQuant/{sample}.log"
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "salmon quant -p {threads} -i {params.idx} -l A -1 {input.f1} -2 {input.f2} --validateMappings -o {params.dir} 2>&1 | tee -a {log}"

rule featureCounts:
    input:
        bam="results/mergedBam/{sample}.Aligned.sortedByCoord.out.bam",
        gtf=config["gtf"]
    output:
        gene_counts="results/counts/featureCounts/{sample}_geneCounts.tab",
        exon_counts="results/counts/featureCounts/{sample}_exonCounts.tab"
    params:
        strand=config["fc_strand"]
    threads: config["ncpus_fc"]
    resources:
        mem_mb=config["mem_fc"],
        runtime_min=config["rt_fc"]
    benchmark:
        "benchmark/featureCounts/{sample}.tsv"
    log:
        "logs/featureCounts/{sample}.log"
    script:
        "scripts/featurecounts.R"

rule TxImport:
    input:
        quant="results/salmon/{sample}/quant.sf",
        gtf=config["gtf"]
    output:
        gene_raw="results/counts/salmon/{sample}_geneCounts.tab",
        gene_scaled="results/counts/salmon/{sample}_geneCounts_scaled.tab",
        transcript_raw="results/counts/salmon/{sample}_transcriptCounts.tab",
        transcript_TPM="results/counts/salmon/{sample}_transcriptTPM.tab"
    threads: 2
    resources:
        mem_mb=config["mem_txi"],
        runtime_min=config["rt_txi"]
    benchmark:
        "benchmark/txImport/{sample}.tsv"
    log:
        "logs/txImport/{sample}.log"
    script:
        "scripts/txImport.R"

rule arriba:
    input:
        bam="results/mergedBam/{sample}.Aligned.sortedByCoord.out.bam",
        reference=config["reference"],
        gtf=config["gtf"],
        bl="resources/blacklist_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz",
        kf="resources/known_fusions_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz",
        pd="resources/protein_domains_hg19_hs37d5_GRCh37_v2.1.0.gff3"
    output:
        fusion="results/fusion/arriba/{sample}.tsv",
        discarded="results/fusion/arriba/{sample}_discarded.tsv"
    threads: 2
    resources:
        mem_mb=config["mem_arriba"],
        runtime_min=config["rt_arriba"]
    benchmark:
        "benchmark/arriba/{sample}.tsv"
    log:
        "logs/arriba/{sample}.log"
    conda:
        "../envs/RNAseq.yaml"
    shell:
        """
        arriba \
        -x {input.bam} \
        -a {input.reference} \
        -g {input.gtf} \
        -b {input.bl} \
        -k {input.kf} \
        -t {input.kf} \
        -p {input.pd} \
        -o {output.fusion} -O {output.discarded} 2>&1 | tee -a {log}
        """

# rule STARfusion:
#     input:
#         bam="results/sortedBams/{sample}.Aligned.sortedByCoord.out.bam",
#         cj="results/STAR_2p/{sample}.Chimeric.out.junction",
#         ctat_lib="resources/ctat+lib"
#     output:
#         fusion="results/fusion/STAR_fusion/{sample}_star-fusion.fusion_predictions.tsv"
#     threads: config["ncpus_STARfusion"]
#     resources:
#         mem_mb=config["mem_STARfusion"],
#         runtime_min=config["rt_STARfusion"]
#     benchmark:
#         "benchmark/STARfusion/{sample}.tsv"
#     log:
#         "logs/STARfusion/{sample}.log"
#     conda:
#         "../envs/RNAseq.yaml"
#     shell:
#         """
#         STAR-Fusion --genome_lib_dir {input.ctat_lib} \
#         -J {input.cj} \
#         --output_dir fusion/STAR_fusion/{wildcards.sample}_
#         """

rule megadepth:
    input:
        bed=config["paSites"],
        bam="results/mergedBam/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        paQuant="results/paQuant/{sample}_paQuant.tab"
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "megadepth {input.bam} --annotation {input.bed} --op sum > {output.paQuant}"