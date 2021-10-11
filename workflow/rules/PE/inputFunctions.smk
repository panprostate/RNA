def star_index(wildcards):
    if config["buildIndex"] == True:
        return "results/index/STARindex_hg19/SA"
    else:
        return "resources/STARindex_hg19/SA"

def salmon_index(wildcards):
    if config["buildIndex"] == True:
        return "results/index/salmon_hg19/ctable.bin"
    else:
        return "resources/salmon_hg19/ctable.bin"

def get_fastq(wildcards):
    return df.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]]

def gather_SJ_filter(wildcards):
    SAMPLES=df.loc[:,"sample_name"]
    UNITS=df.loc[:,"unit_name"]
    return expand("results/STAR_1p/{sample}_{unit}SJ.out.tab", zip, sample=SAMPLES, unit=UNITS)

def gather_SJ(wildcards):
    if config["SJ_filter"] == False:
        SAMPLES=df.loc[:,"sample_name"]
        UNITS=df.loc[:,"unit_name"]
        return expand("results/STAR_1p/{sample}_{unit}SJ.out.tab", zip, sample=SAMPLES, unit=UNITS)
    else:
        SAMPLES=df.loc[:,"sample_name"]
        UNITS=df.loc[:,"unit_name"]
        return expand("results/STAR_1p/{sample}_{unit}.filtered.SJ.out.tab", zip, sample=SAMPLES, unit=UNITS)

def gather_bams(wildcards):
    UNITS=df.loc[wildcards.sample, "unit_name"]
    bams=expand("results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam", unit=UNITS, sample=wildcards.sample)
    return bams

def gather_chims(wildcards):
    UNITS = df.loc[wildcards.sample, "unit_name"]
    chims = expand("results/STAR_2p/{sample}_{unit}Chimeric.out.junction", unit=UNITS, sample=wildcards.sample)
    return chims

def gather_salmon_input1(wildcards):
    UNITS=df.loc[wildcards.sample, "unit_name"]
    fastq1=expand("results/trim/{{sample}}/{{sample}}_{unit}_R1_trimmed.fastq.gz", unit=UNITS)
    return fastq1

def gather_salmon_input2(wildcards):
    UNITS=df.loc[wildcards.sample, "unit_name"]
    fastq2=expand("results/trim/{{sample}}/{{sample}}_{unit}_R2_trimmed.fastq.gz", unit=UNITS)
    return fastq2

def VC_gather_bams(wildcards):
    UNITS=df.loc[wildcards.sample, "unit_name"]
    bams=expand("results/variantCalling/mergedBam/{sample}_{unit}.bam", unit=UNITS, sample=wildcards.sample)
    return bams
