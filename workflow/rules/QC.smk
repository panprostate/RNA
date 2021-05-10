rule rseqc_gtf2bed:
    input:
        gtf=config["gtf"],
    output:
        bed="results/qc/rseqc/annotation.bed",
        db=temp("results/qc/rseqc/annotation.db"),
    log:
        "logs/rseqc_gtf2bed.log",
    conda:
        "../envs/RNAseq.yaml"
    script:
        "../scripts/gtf2bed.py"


rule rseqc_junction_annotation:
    input:
        bam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.junctionanno.junction.bed",
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_annotation/{sample}_{unit}.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix="results/qc/rseqc/{sample}_{unit}.junctionanno",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"


rule rseqc_junction_saturation:
    input:
        bam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.junctionsat.junctionSaturation_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_saturation/{sample}_{unit}.log",
    params:
        extra=r"-q 255",
        prefix="results/qc/rseqc/{sample}_{unit}.junctionsat",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule rseqc_stat:
    input:
        "results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
    output:
        "results/qc/rseqc/{sample}_{unit}.stats.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_stat/{sample}_{unit}.log",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"


rule rseqc_infer:
    input:
        bam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.infer_experiment.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_infer/{sample}_{unit}.log",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_innerdis:
    input:
        bam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.inner_distance_freq.inner_distance.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_innerdis/{sample}_{unit}.log",
    params:
        prefix="results/qc/rseqc/{sample}_{unit}.inner_distance_freq",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam="results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.readdistribution.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_readdis/{sample}_{unit}.log",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
    input:
        "results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
    output:
        "results/qc/rseqc/{sample}_{unit}.readdup.DupRate_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_readdup/{sample}_{unit}.log",
    params:
        prefix="results/qc/rseqc/{sample}_{unit}.readdup",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
    input:
        "results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam",
    output:
        "results/qc/rseqc/{sample}_{unit}.readgc.GC_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_readgc/{sample}_{unit}.log",
    params:
        prefix="results/qc/rseqc/{sample}_{unit}.readgc",
    conda:
        "../envs/RNAseq.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"


rule multiqc:
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
            expand("logs/rseqc/rseqc_junction_annotation/{sample}_{unit}.log",zip, sample=samples_trim, unit=unit_trim)
    output:
        "results/qc/multiqc_report.html"
    log:
        "logs/multiqc.log",
    wrapper:
        "0.31.1/bio/multiqc"
