def get_fastq(wildcards):
    return df.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]]

def gather_SJ(wildcards):
    SAMPLES=df.loc[:,"sample_name"]
    UNITS=df.loc[:,"unit_name"]
    return expand("results/STAR_1p/{sample}_{unit}.SJ.out.tab", zip, sample=SAMPLES, unit=UNITS)

def gather_bams(wildcards):
    UNITS=df.loc[wildcards.sample, "unit_name"]
    bams=expand("results/sortedBams/{sample}_{unit}.Aligned.sortedByCoord.out.bam", unit=UNITS, sample=wildcards.sample)
    return bams

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