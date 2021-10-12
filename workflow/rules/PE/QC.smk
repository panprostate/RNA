rule rseqc_gtf2bed:
    input:
        gtf="resources/gencode.v38lift37.annotation.gtf",
    output:
        bed="results/qc/rseqc/annotation.bed",
        db=temp("results/qc/rseqc/annotation.db")
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/rseqc/rseqc_gtf2bed.log",
    conda:
        "../../envs/RNAseq.yaml"
    script:
        "../../scripts/gtf2bed.py"

rule fastqc:
    input:
        bam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam"
    output:
        "results/qc/fastqc/{sample}/{sample}_{unit}_fastqc.zip"
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/fastqc/{sample}_{unit}.log",
    params:
        dir="results/qc/fastqc/{sample}/",
    conda:
        "../../envs/RNAseq.yaml"
    shell:
        """
        fastqc -o {params.dir} {input.bam}
        mv results/qc/fastqc/{wildcards.sample}/{wildcards.sample}_{wildcards.unit}.Aligned.sortedByCoord.out_fastqc.zip {output} 
        """

rule rseqc_junction_annotation:
    input:
        bam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.junctionanno.junction.bed"
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/rseqc/rseqc_junction_annotation/{sample}_{unit}.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix="results/qc/rseqc/{sample}_{unit}.junctionanno",
    conda:
        "../../envs/RNAseq.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"


rule rseqc_junction_saturation:
    input:
        bam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.junctionsat.junctionSaturation_plot.pdf",
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/rseqc/rseqc_junction_saturation/{sample}_{unit}.log",
    params:
        extra=r"-q 255",
        prefix="results/qc/rseqc/{sample}_{unit}.junctionsat",
    conda:
        "../../envs/RNAseq.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule rseqc_stat:
    input:
        "results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
    output:
        "results/qc/rseqc/{sample}_{unit}.stats.txt",
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/rseqc/rseqc_stat/{sample}_{unit}.log",
    conda:
        "../../envs/RNAseq.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"


rule rseqc_infer:
    input:
        bam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.infer_experiment.txt",
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/rseqc/rseqc_infer/{sample}_{unit}.log",
    conda:
        "../../envs/RNAseq.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_innerdis:
    input:
        bam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.inner_distance_freq.inner_distance.txt",
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/rseqc/rseqc_innerdis/{sample}_{unit}.log",
    params:
        prefix="results/qc/rseqc/{sample}_{unit}.inner_distance_freq",
    conda:
        "../../envs/RNAseq.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.readdistribution.txt",
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/rseqc/rseqc_readdis/{sample}_{unit}.log",
    conda:
        "../../envs/RNAseq.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
    input:
        "results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
    output:
        "results/qc/rseqc/{sample}_{unit}.readdup.DupRate_plot.pdf",
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/rseqc/rseqc_readdup/{sample}_{unit}.log",
    params:
        prefix="results/qc/rseqc/{sample}_{unit}.readdup",
    conda:
        "../../envs/RNAseq.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
    input:
        "results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
    output:
        "results/qc/rseqc/{sample}_{unit}.readgc.GC_plot.pdf",
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    log:
        "logs/rseqc/rseqc_readgc/{sample}_{unit}.log",
    params:
        prefix="results/qc/rseqc/{sample}_{unit}.readgc",
    conda:
        "../../envs/RNAseq.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"


rule bamqc:
    input:
        expand("results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",zip, sample=samples_trim, unit=unit_trim),
        expand("results/qc/rseqc/{sample}_{unit}.junctionanno.junction.bed",zip, sample=samples_trim, unit=unit_trim),
        expand("results/qc/rseqc/{sample}_{unit}.junctionsat.junctionSaturation_plot.pdf",zip, sample=samples_trim, unit=unit_trim),
        expand("results/qc/rseqc/{sample}_{unit}.infer_experiment.txt",zip, sample=samples_trim, unit=unit_trim),
        expand("results/qc/rseqc/{sample}_{unit}.stats.txt",zip, sample=samples_trim, unit=unit_trim),
        expand("results/qc/rseqc/{sample}_{unit}.inner_distance_freq.inner_distance.txt",zip, sample=samples_trim, unit=unit_trim),
        expand("results/qc/rseqc/{sample}_{unit}.readdistribution.txt",zip, sample=samples_trim, unit=unit_trim),
        expand("results/qc/rseqc/{sample}_{unit}.readdup.DupRate_plot.pdf",zip, sample=samples_trim, unit=unit_trim),
        expand("results/qc/rseqc/{sample}_{unit}.readgc.GC_plot.pdf",zip, sample=samples_trim, unit=unit_trim),
        expand("results/STAR_2p/{sample}_{unit}Log.final.out",zip, sample=samples_trim, unit=unit_trim)
    output:
        out="results/qc/bam_qc.html",
        bam_stat="results/qc/bam_qc_data/multiqc_rseqc_bam_stat.txt",
        bam_infer_exp="results/qc/bam_qc_data/multiqc_rseqc_infer_experiment.txt",
        bam_read_distribution="results/qc/bam_qc_data/multiqc_rseqc_read_distribution.txt",
        star_stat="results/qc/bam_qc_data/multiqc_star.txt",
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    params:
    	outDir="results/qc",
        name="bam_qc.html"
    log:
        "logs/bam_qc.log",
    conda:
        "../../envs/RNAseq.yaml"
    shell:
        "multiqc"
        " --force"
        " -o {params.outDir}"
        " -n {params.name}"
        " {input}"

rule sample_qc:
    input:
        gcg=expand("results/counts/featureCounts/{sample}_geneCounts_gencode.tsv.summary", sample=SAMPLES),
        gcf=expand("results/counts/featureCounts/{sample}_geneCounts_fc.tsv.summary", sample=SAMPLES),
        ecg=expand("results/counts/featureCounts/{sample}_exonCounts_gencode.tsv.summary", sample=SAMPLES),
        ecf=expand("results/counts/featureCounts/{sample}_exonCounts_fc.tsv.summary", sample=SAMPLES),
        salmon=expand("results/salmon/{sample}/aux_info/meta_info.json",sample=SAMPLES),
        salmon_length=expand("results/salmon/{sample}/libParams/flenDist.txt",sample=SAMPLES) 
    output:
        geneCount_gc="results/qc/featureCounts_gene_gencode.html",
        geneCount_fc="results/qc/featureCounts_gene_fc.html",
        exonCount_gc="results/qc/featureCounts_exon_gencode.html",
        exonCount_fc="results/qc/featureCounts_exon_fc.html",
        salmon="results/qc/salmon.html"
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    params:
    	outDir="results/qc",
        geneCount_gc="featureCounts_gene_gencode.html",
        geneCount_fc="featureCounts_gene_fc.html",
        exonCount_gc="featureCounts_exon_gencode.html",
        exonCount_fc="featureCounts_exon_fc.html",
        salmon="salmon.html"
    log:
        "logs/samples_qc.log",
    conda:
        "../../envs/RNAseq.yaml"
    shell:
        """
        set -e
        multiqc --force -o {params.outDir} -n {params.geneCount_gc} {input.gcg}
        multiqc --force -o {params.outDir} -n {params.geneCount_fc} {input.gcf}
        multiqc --force -o {params.outDir} -n {params.exonCount_gc} {input.ecg}
        multiqc --force -o {params.outDir} -n {params.exonCount_fc} {input.ecf}
        multiqc --force -o {params.outDir} -n {params.salmon} {input.salmon} {input.salmon_length}
        """

rule reads_qc:
    input:
        expand("results/qc/fastqc/{sample}/{sample}_{unit}_fastqc.zip",zip, sample=samples_trim, unit=unit_trim)
    output:
        out="results/qc/reads_qc.html",
        reads_dup="results/qc/reads_qc_data/mqc_fastqc_sequence_counts_plot_1.txt",
        reads_dup_distribution="results/qc/reads_qc_data/mqc_fastqc_sequence_duplication_levels_plot_1.txt",
        reads_length_distribution="results/qc/reads_qc_data/mqc_fastqc_sequence_length_distribution_plot_1.txt",
        reads_quality_scores="results/qc/reads_qc_data/mqc_fastqc_per_sequence_quality_scores_plot_1.txt",
        reads_qc_distribution="results/qc/reads_qc_data/mqc_fastqc_per_sequence_gc_content_plot_Percentages.txt",
    resources:
        mem_mb=config["mem_qc"],
        runtime_min=config["rt_qc"]
    params:
    	outDir="results/qc",
        name="reads_qc.html"
    log:
        "logs/read_qc.log",
    conda:
        "../../envs/RNAseq.yaml"
    shell:
        "multiqc"
        " --force"
        " -o {params.outDir}"
        " -n {params.name}"
        " {input}"