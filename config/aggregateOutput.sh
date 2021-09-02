#!/usr/bin/sh
set -e

if [ -e results/COMPLETED ]
then
    mkdir -p shared/counts/gencode/gene/featureCounts
    mkdir -p shared/counts/gencode/gene/salmon/raw
    mkdir -p shared/counts/gencode/gene/salmon/scaled
    mkdir -p shared/counts/gencode/exon
    mkdir -p shared/counts/gencode/transcript/raw
    mkdir -p shared/counts/gencode/transcript/TPM
    mkdir -p shared/counts/fc/gene
    mkdir -p shared/counts/fc/exon
    mkdir -p shared/counts/pa
    mkdir -p shared/fusion/arriba
    mkdir -p shared/fusion/starfusion
    mkdir -p shared/VariantCalling
    mkdir -p shared/QC/featureCounts
    mv results/counts/featureCounts/*_geneCounts_gencode.tsv shared/counts/gencode/gene/featureCounts
    mv results/counts/featureCounts/*_geneCounts_fc.tsv shared/counts/fc/gene/
    mv results/counts/featureCounts/*_exonCounts_gencode.tsv shared/counts/gencode/exon/
    mv results/counts/featureCounts/*_exonCounts_fc.tsv shared/counts/fc/exon/
    mv results/counts/salmon/*_geneCounts_scaled.tsv shared/counts/gencode/gene/salmon/scaled
    mv results/counts/salmon/*_geneCounts.tsv shared/counts/gencode/gene/salmon/raw
    mv results/counts/salmon/*_transcriptCounts.tsv shared/counts/gencode/transcript/raw
    mv results/counts/salmon/*_transcriptTPM.tsv shared/counts/gencode/transcript/TPM
    mv results/paQuant/*_paQuant.tsv.gz shared/counts/pa
    mv results/fusion/arriba/*.tsv shared/fusion/arriba
    mv results/fusion/STAR_fusion/* shared/fusion/starfusion
    mv results/variantCalling/vcf_filtered/*.vcf.gz shared/VariantCalling
    mv results/qc/samples_qc.html results/qc/reads_qc.html results/qc/bam_qc.html shared/QC
    mv results/counts/featureCounts/*.summary shared/QC/featureCounts
    gzip -r shared/*
else
    echo "*** PIPELINE HAS NOT BEEN COMPLETED! MOVING FILES WILL BREAK SNAKEMAKE FILE TRACKING – SO MAKE SURE IT IS COMPLETE BEFORE MOVING ***"
fi
