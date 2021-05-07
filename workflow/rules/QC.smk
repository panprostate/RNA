rule rseqc_gtf2bed:
    input:
        gtf=config["gtf"],
    output:
        bed="results/qc/rseqc/annotation.bed",
        db=temp("results/qc/rseqc/annotation.db"),
    log:
        "logs/rseqc_gtf2bed.log",
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"


rule rseqc_junction_annotation:
    input:
        bam="results/STAR_2p/{sample}.Aligned.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}.junctionanno.junction.bed",
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_annotation/{sample}.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix=lambda w, output: strip_suffix(output[0], ".junction.bed"),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"


rule rseqc_junction_saturation:
    input:
        bam="results/STAR_2p/{sample}.Aligned.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}.junctionsat.junctionSaturation_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_saturation/{sample}.log",
    params:
        extra=r"-q 255",
        prefix=lambda w, output: strip_suffix(
            output[0], ".junctionSaturation_plot.pdf"
        ),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule rseqc_stat:
    input:
        "results/STAR_2p/{sample}.Aligned.out.bam",
    output:
        "results/qc/rseqc/{sample}.stats.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_stat/{sample}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"


rule rseqc_infer:
    input:
        bam="results/STAR_2p/{sample}.Aligned.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}.infer_experiment.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_infer/{sample}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_innerdis:
    input:
        bam="results/STAR_2p/{sample}.Aligned.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}.inner_distance_freq.inner_distance.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_innerdis/{sample}.log",
    params:
        prefix=lambda w, output: strip_suffix(output[0], ".inner_distance.txt"),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam="results/STAR_2p/{sample}.Aligned.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}.readdistribution.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_readdis/{sample}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
    input:
        "results/STAR_2p/{sample}.Aligned.out.bam",
    output:
        "results/qc/rseqc/{sample}.readdup.DupRate_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_readdup/{sample}.log",
    params:
        prefix=lambda w, output: strip_suffix(output[0], ".DupRate_plot.pdf"),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
    input:
        "results/STAR_2p/{sample}.Aligned.out.bam",
    output:
        "results/qc/rseqc/{sample}.readgc.GC_plot.pdf",
    priority: 1
    log:
        "logs/rseqc/rseqc_readgc/{sample}.log",
    params:
        prefix=lambda w, output: strip_suffix(output[0], ".GC_plot.pdf"),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"


rule multiqc:
    input:
        expand(
            "results/STAR_2p/{sample}.Aligned.out.bam",
            sample=SAMPLES,
        ),
        expand(
            "results/qc/rseqc/{sample}.junctionanno.junction.bed",
            sample=SAMPLES,
        ),
        expand(
            "results/qc/rseqc/{sample}.junctionsat.junctionSaturation_plot.pdf",
            sample=SAMPLES,
        ),
        expand(
            "results/qc/rseqc/{sample}.infer_experiment.txt",
            sample=SAMPLES,
        ),
        expand(
            "results/qc/rseqc/{sample}.stats.txt",
            sample=SAMPLES,
        ),
        expand(
            "results/qc/rseqc/{sample}.inner_distance_freq.inner_distance.txt",
            sample=SAMPLES,
        ),
        expand(
            "results/qc/rseqc/{sample}.readdistribution.txt",
            sample=SAMPLES,
        ),
        expand(
            "results/qc/rseqc/{sample}.readdup.DupRate_plot.pdf",
            sample=SAMPLES,
        ),
        expand(
            "results/qc/rseqc/{sample}.readgc.GC_plot.pdf",
            sample=SAMPLES,
        ),
        expand(
            "logs/rseqc/rseqc_junction_annotation/{sample}.log",
            sample=SAMPLES,
        ),
    output:
        "results/qc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    wrapper:
        "0.31.1/bio/multiqc"
