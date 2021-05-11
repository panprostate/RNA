options(scipen=999)
require(tximport)
require(GenomicFeatures)
transcripts <- read.table(snakemake@input[["quant"]], sep ="\t", header = T)
write.csv(transcripts[,c("Name", "NumReads")], file=snakemake@output[["transcript_raw"]], col.names = T, row.names = F)
write.csv(transcripts[,c("Name", "TPM")], file=snakemake@output[["transcript_TPM"]], col.names = T, row.names = F)

tx2gene <- read.table(snakemake@input[["tx2gene"]])

### Read counts and summarize to gene level
txi <- tximport(snakemake@input[["quant"]], type = "salmon", tx2gene = tx2gene)
txiScaled <- tximport(snakemake@input[["quant"]], type = "salmon", tx2gene = tx2gene,countsFromAbundance = "lengthScaledTPM")
### Rename columns
nms <- basename(gsub("\\/quant.sf", "", snakemake@input[["quant"]]))
mat <- txi$counts
matScaled <- txiScaled$counts
colnames(mat) <- nms
colnames(matScaled) <- nms

write.csv(mat, file=snakemake@output[["gene_raw"]], row.names = T, col.names = T)
write.csv(matScaled, file=snakemake@output[["gene_scaled"]], row.names = T, col.names = T)
