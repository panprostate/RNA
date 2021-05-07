def get_fastq(wildcards):
    return df.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]]

def gather_SJ(wildcards):
    SAMPLES=df.loc[:,"sample_name"]
    UNITS=df.loc[:,"unit_name"]
    return expand("results/STAR_1p/{sample}_{unit}.SJ.out.tab", zip, sample=SAMPLES, unit=UNITS)

def gather_salmon_input1(wildcards):
    UNITS=df.loc[wildcards.sample, ["unit_name"]]
    fastq1=expand("results/trim/{{sample}}/{{sample}}_{unit}_R1_trimmed.fastq.gz", unit=UNITS)
    return fastq1

def gather_salmon_input2(wildcards):
    UNITS=df.loc[wildcards.sample, ["unit_name"]]
    fastq2=expand("results/trim/{{sample}}/{{sample}}_{unit}_R2_trimmed.fastq.gz", unit=UNITS)
    return fastq2

def gather_STAR_input1(wildcards):
    UNITS=df.loc[wildcards.sample, ["unit_name"]]
    fastq1=expand("results/trim/{{sample}}/{{sample}}_{unit}_R1_trimmed.fastq.gz", unit=UNITS)
    fastq1=",".join(fastq1)
    return fastq1

def gather_STAR_input2(wildcards):
    UNITS=df.loc[wildcards.sample, ["unit_name"]]
    fastq2=expand("results/trim/{{sample}}/{{sample}}_{unit}_R2_trimmed.fastq.gz", unit=UNITS)
    fastq2=",".join(fastq2)
    return fastq2

def get_RG_groups(wildcards):
    UNITS=df.loc[wildcards.sample, ["unit_name"]]
    fastq1=expand("results/trim/{{sample}}/{{sample}}_{unit}_R1_trimmed.fastq.gz", unit=UNITS)
    RG=["ID:"+re.sub("_R1_trimmed.fastq.gz", "", i) for i in fastq1]
    RG=" , ".join(RG)
    return RG