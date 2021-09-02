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
        "../../envs/RNAseq.yaml"
    priority: 1
    shell:
        "cutadapt"
        " {params.adapters}"
        " -m 31"
        " --compression-level {params.compression}"
        " -o {output.fastq1}"
        " -p {output.fastq2}"
        " -j {threads}"
        " {input} > {output.qc} 2>>{log}"

rule star_index_build:
    input:
        reference="resources/Homo_sapiens_assembly19_1000genomes_decoy.fasta",
        gtf="resources/gencode.v38lift37.annotation.gtf"
    output:
        index="results/index/STARindex_hg19/SA"
    params:
        idx="results/index/STARindex_hg19/"
    threads: config["ncpus_STAR_index"]
    resources:
        mem_mb=config["mem_STAR_index"],
        runtime_min=config["rt_STAR_index"]
    benchmark:
        "benchmark/STARindex/STARindex.tsv"
    log:
        "logs/STARindex/STARindex.log"
    conda:
        "../../envs/RNAseq.yaml"
    priority: 2
    shell:
        "STAR"
        " --runMode genomeGenerate"
        " --genomeDir {params.idx}"
        " --genomeFastaFiles {input.reference}"
        " --sjdbOverhang 250"
        " --sjdbGTFfile {input.gtf}"
        " --runThreadN {threads} 2>&1 | tee -a {log}"

rule salmon_index_build:
    input:
        transcripts="resources/gencode.v38lift37.transcripts.fa",
        reference="resources/Homo_sapiens_assembly19_1000genomes_decoy.fasta"
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
        "../../envs/RNAseq.yaml"
    priority: 2
    shell:
        """
        set -e
        cat {input.transcripts} {input.reference} > resources/gentrome.fa
        grep "^>" {input.reference} | cut -d " " -f 1 > resources/decoys.txt
        sed -i -e 's/>//g' resources/decoys.txt
        salmon index -t resources/gentrome.fa -i {params.idx} -p {threads} --decoys resources/decoys.txt --gencode -k 31 2>&1 | tee -a {log}
        rm -f resources/gentrome.fa resources/decoys.txt
        """

rule star_1pass:
    input:
        index=star_index,
        f1="results/trim/{sample}/{sample}_{unit}_R1_trimmed.fastq.gz",
        f2="results/trim/{sample}/{sample}_{unit}_R2_trimmed.fastq.gz"
    output:
        SJ="results/STAR_1p/{sample}_{unit}SJ.out.tab",
        outputs=temp(multiext("results/STAR_1p/{sample}_{unit}", "Log.progress.out", "Log.final.out", "Log.out"))
    params:
        idx=lambda wildcards, input: input.index[:-2]
    threads: config["ncpus_STAR_align"]
    resources:
        mem_mb=config["mem_STAR_align"],
        runtime_min=config["rt_STAR_align"]
    benchmark:
        "benchmark/STAR_1p/{sample}_{unit}.tsv"
    log:
        "logs/STAR_1p/{sample}_{unit}.log"
    conda:
        "../../envs/RNAseq.yaml"
    priority: 1
    shell:
        """
        STAR \
        --genomeDir {params.idx} \
        --readFilesIn {input.f1} {input.f2} \
        --runThreadN {threads} \
        --readFilesCommand zcat \
        --outFilterMultimapScoreRange 1 \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 10 \
        --outFilterType BySJout \
        --alignIntronMax 500000 \
        --alignIntronMin 20 \
        --alignMatesGapMax 1000000 \
        --sjdbScore 2 \
        --alignSJDBoverhangMin 5 \
        --genomeLoad NoSharedMemory \
        --outFilterMatchNminOverLread 0.1 \
        --outFilterScoreMinOverLread 0.1 \
        --sjdbOverhang 250 \
        --outSAMstrandField intronMotif \
        --peOverlapNbasesMin 10 \
        --alignSplicedMateMapLminOverLmate 0.5 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --outSAMtype None \
        --outSAMmode None \
        --outFileNamePrefix results/STAR_1p/{wildcards.sample}_{wildcards.unit} 2>&1 | tee -a {log}
        """

rule star_2pass:
    input:
        f1="results/trim/{sample}/{sample}_{unit}_R1_trimmed.fastq.gz",
        f2="results/trim/{sample}/{sample}_{unit}_R2_trimmed.fastq.gz",
        index=star_index,
        SJ=gather_SJ
    output:
        bam=temp("results/STAR_2p/{sample}_{unit}Aligned.out.bam"),
        star_logs=multiext("results/STAR_2p/{sample}_{unit}", "SJ.out.tab", "Log.final.out", "Log.out", "Chimeric.out.junction")
    params:
        compression=config["compression_level"],
        idx=lambda wildcards, input: input.index[:-2]
    threads: config["ncpus_STAR_align"]
    resources:
        mem_mb=config["mem_STAR_align"],
        runtime_min=config["rt_STAR_align"]
    benchmark:
        "benchmark/STAR_2p/{sample}_{unit}.tsv"
    log:
        "logs/STAR_2p/{sample}_{unit}.log"
    conda:
        "../../envs/RNAseq.yaml"
    priority: 1
    shell:
        """
        STAR \
        --genomeDir {params.idx} \
        --readFilesIn {input.f1} {input.f2} \
        --runThreadN {threads} \
        --readFilesCommand zcat \
        --outBAMcompression {params.compression} \
        --outFilterType BySJout \
        --outFilterMultimapScoreRange 1 \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 10 \
        --alignIntronMax 500000 \
        --alignIntronMin 20 \
        --alignMatesGapMax 1000000 \
        --sjdbScore 2 \
        --outSAMtype BAM Unsorted \
        --alignSJDBoverhangMin 5 \
        --genomeLoad NoSharedMemory \
        --sjdbFileChrStartEnd {input.SJ} \
        --outFilterMatchNminOverLread 0.1 \
        --outFilterScoreMinOverLread 0.1 \
        --sjdbOverhang 250 \
        --outSAMstrandField intronMotif \
        --peOverlapNbasesMin 10 \
        --alignSplicedMateMapLminOverLmate 0.5 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimSegmentMin 10 \
        --chimOutJunctionFormat 1 \
        --chimOutType Junctions WithinBAM HardClip \
        --chimJunctionOverhangMin 10 \
        --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 \
        --chimScoreSeparation 1 \
        --chimSegmentReadGapMax 3 \
        --chimMultimapNmax 20 \
        --outFileNamePrefix results/STAR_2p/{wildcards.sample}_{wildcards.unit} 2>&1 | tee -a {log}
        """

rule sortBam:
    input:
        bam="results/STAR_2p/{sample}_{unit}Aligned.out.bam"
    output:
        sbam=temp("results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam"),
        bai=temp("results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam.bai")
    threads: config["ncpus_sortBam"]
    resources:
        mem_mb=config["mem_sortBam"],
        runtime_min=config["rt_sortBam"]
    benchmark:
        "benchmark/sortBam/{sample}_{unit}.tsv"
    log:
        "logs/sortBam/{sample}_{unit}.log"
    conda:
        "../../envs/RNAseq.yaml"
    priority: 1
    shell:
        """
        samtools sort -@ {threads} {input.bam} -o {output.sbam} 2>>{log}
        samtools index {output.sbam}
        """

rule mergeBam:
    input:
        bams=gather_bams,
        chims=gather_chims
    output:
        mbam="results/mergedBam/{sample}.Aligned.sortedByCoord.out.bam",
        bai="results/mergedBam/{sample}.Aligned.sortedByCoord.out.bam.bai",
        chim="results/mergedBam/{sample}Chimeric.out.junction"
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
        "../../envs/RNAseq.yaml"
    priority: 1
    shell:
        """
        set -e
        INPUT=({input.bams})
        if ((${{#INPUT[@]}} == 1)); then
            cp {input.bams} {output.mbam}
            cp {input.bams}.bai {output.mbam}.bai
            cp {input.chims} {output.chim}
        else
            samtools merge -f -1 {output.mbam} {input.bams} 2>>{log}
            samtools index {output.mbam}
            # Also merge the chimeric files for star-fusion
            cat {input.chims[0]} | head -n -2 > {output.chim}
            files=({input.chims})
            unset files[0]
            files=("${{files[@]}}")
            rm -f results/mergedBam/{wildcards.sample}chims results/mergedBam/{wildcards.sample}counts
            touch results/mergedBam/{wildcards.sample}chims
            cat {input.chims[0]} | tail -n 1 > results/mergedBam/{wildcards.sample}counts
            for file in $files; do
                cat $file | tail -n +2 | head -n -2 >> results/mergedBam/{wildcards.sample}chims
                cat $file | tail -n 1 >> results/mergedBam/{wildcards.sample}counts
            done
            cat results/mergedBam/{wildcards.sample}chims >> {output.chim}
            cat {input.chims[0]} | tail -n 2 |head -n 1 >> {output.chim}
            awk -F' ' '{{Nreads+=$3; NreadsUnique+=$5; NreadsMulti+=$7}}END{{print "# Nreads " Nreads "\t" "NreadsUnique " NreadsUnique "\t" "NreadsMulti " NreadsMulti}}' results/mergedBam/{wildcards.sample}counts >> {output.chim}
            rm -f results/mergedBam/{wildcards.sample}chims results/mergedBam/{wildcards.sample}counts
        fi
        """

rule salmon_quant:
    input:
        index=salmon_index,
        f1=gather_salmon_input1,
        f2=gather_salmon_input2
    output:
        tc="results/salmon/{sample}/quant.sf",
        flen="results/salmon/{sample}/libParams/flenDist.txt",
        meta="results/salmon/{sample}/aux_info/meta_info.json"
    params:
        idx=lambda wildcards, input: input.index[:-10],
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
        "../../envs/RNAseq.yaml"
    priority: 1
    shell:
        "salmon quant -p {threads} -i {params.idx} -l A -1 {input.f1} -2 {input.f2} --validateMappings --rangeFactorizationBins 4 --gcBias -o {params.dir} 2>&1 | tee -a {log}"

rule featureCounts:
    input:
        bam="results/mergedBam/{sample}.Aligned.sortedByCoord.out.bam",
        gtf="resources/gencode.v38lift37.annotation.gtf",
        fc="resources/FANTOM_CAT.lv3_robust.gtf"
    output:
        gene_counts="results/counts/featureCounts/{sample}_geneCounts_gencode.tsv",
        exon_counts="results/counts/featureCounts/{sample}_exonCounts_gencode.tsv",
        stats="results/counts/featureCounts/{sample}_geneCounts_gencode.tsv.summary",
        gene_counts_fc = "results/counts/featureCounts/{sample}_geneCounts_fc.tsv",
        exon_counts_fc = "results/counts/featureCounts/{sample}_exonCounts_fc.tsv",
        stats_fc = "results/counts/featureCounts/{sample}_geneCounts_fc.tsv.summary"
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
    conda:
        "../../envs/RNAseq.yaml"
    priority: 1
    shell:
        """
        featureCounts -p -a {input.gtf} -T {threads} -s {params.strand} -t exon -g gene_id -o {output.gene_counts} {input.bam}
        featureCounts -p -a {input.fc} -T {threads} -s {params.strand} -t exon -g gene_id -o {output.gene_counts_fc} {input.bam}
        featureCounts -p -f -O -a {input.gtf} -T {threads} -s {params.strand} -t exon -g gene_id -o {output.exon_counts} {input.bam}
        featureCounts -p -f -O -a {input.fc} -T {threads} -s {params.strand} -t exon -g gene_id -o {output.exon_counts_fc} {input.bam}
        """

rule TxImport:
    input:
        quant="results/salmon/{sample}/quant.sf",
        tx2gene="resources/tx2gene.tsv.gz",
        gtf="resources/gencode.v38lift37.annotation.gtf"
    output:
        gene_raw="results/counts/salmon/{sample}_geneCounts.tsv",
        gene_scaled="results/counts/salmon/{sample}_geneCounts_scaled.tsv",
        transcript_raw="results/counts/salmon/{sample}_transcriptCounts.tsv",
        transcript_TPM="results/counts/salmon/{sample}_transcriptTPM.tsv"
    threads: 2
    resources:
        mem_mb=config["mem_txi"],
        runtime_min=config["rt_txi"]
    benchmark:
        "benchmark/txImport/{sample}.tsv"
    log:
        "logs/txImport/{sample}.log"
    conda:
        "../../envs/RNAseq.yaml"
    priority: 1
    script:
        "../../scripts/txImport.R"

rule arriba:
    input:
        bam="results/mergedBam/{sample}.Aligned.sortedByCoord.out.bam",
        reference="resources/Homo_sapiens_assembly19_1000genomes_decoy.fasta",
        gtf="resources/gencode.v38lift37.annotation.gtf",
        bl="resources/blacklist_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz",
        kf="resources/known_fusions_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz",
        pd="resources/protein_domains_hg19_hs37d5_GRCh37_v2.1.0.gff3"
    output:
        fusion="results/fusion/arriba/{sample}.tsv",
        discarded="results/fusion/arriba/{sample}_discarded.tsv.gz"
    threads: 1
    resources:
        mem_mb=config["mem_arriba"],
        runtime_min=config["rt_arriba"]
    benchmark:
        "benchmark/arriba/{sample}.tsv"
    log:
        "logs/arriba/{sample}.log"
    conda:
        "../../envs/RNAseq.yaml"
    priority: 1
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
        -o {output.fusion} -O results/fusion/arriba/{wildcards.sample}_discarded.tsv 2>&1 | tee -a {log}
        gzip -9 -c results/fusion/arriba/{wildcards.sample}_discarded.tsv > {output.discarded}
        rm -f results/fusion/arriba/{wildcards.sample}_discarded.tsv
        """

rule STARfusion:
    input:
        cj="results/mergedBam/{sample}Chimeric.out.junction",
        ctat_lib = "resources/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa"
    output:
        fusion="results/fusion/STAR_fusion/{sample}/star-fusion.fusion_predictions.tsv"
    params:
        ctat_path = "resources/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir"
    threads: config["ncpus_STARfusion"]
    resources:
        mem_mb=config["mem_STARfusion"],
        runtime_min=config["rt_STARfusion"]
    benchmark:
        "benchmark/STARfusion/{sample}.tsv"
    log:
        "logs/STARfusion/{sample}.log"
    conda:
        "../../envs/STARfusion.yaml"
    priority: 1
    shell:
        """
        workflow/scripts/fixChim.awk < {input.cj} > {input.cj}.fix
        mv -f {input.cj}.fix {input.cj}
        STAR-Fusion --genome_lib_dir {params.ctat_path} \
        -J {input.cj} \
        --output_dir results/fusion/STAR_fusion/{wildcards.sample}
        """

rule megadepth:
    input:
        bed="resources/pasites_hg19.bed",
        bam="results/mergedBam/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        paQuant="results/paQuant/{sample}_paQuant.tsv.gz"
    threads: config["ncpus_megadepth"]
    resources:
        mem_mb=config["mem_megadepth"],
        runtime_min=config["rt_megadepth"]
    conda:
        "../../envs/RNAseq.yaml"
    priority: 1
    shell:
        """
        megadepth {input.bam} --annotation {input.bed} --op sum > results/paQuant/{wildcards.sample}_paQuant.tsv
        gzip results/paQuant/{wildcards.sample}_paQuant.tsv
        """
