if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("Rsubread")) {
    BiocManager::install("Rsubread")
}
require(Rsubread)

gene_quant <- featureCounts(files=snakemake@input["bam"], annot.ext = snakemake@input["gtf"], GTF.attrType = "gene_id",  isGTFAnnotationFile = T, strandSpecific = snakemake@params["strand"], isPairedEnd = T, countReadPairs = T,  nthreads = snakemake@threads, requireBothEndsMapped = F)
exon_quant <- featureCounts(files=snakemake@input["bam"], annot.ext = snakemake@input["gtf"], GTF.attrType = "exon_id",  isGTFAnnotationFile = T, strandSpecific = snakemake@params["strand"], isPairedEnd = T, countReadPairs = T,  nthreads = snakemake@threads, requireBothEndsMapped = F)

write.table(gene_quant$stats, sep ="\t", file=snakemake@output["stats"], col.names=T, row.names=F)
write.csv(gene_quant$counts, file=snakemake@output["gene_counts"],  col.names=T, row.names=T)
write.csv(exon_quant$counts, file=snakemake@output["exon_counts"],  col.names=T, row.names=T)