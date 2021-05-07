gtf <- read.table(snakemake@input["gtf"], sep="\t")
gtf <- subset(gtf, V3 == "exon")
gtf <- data.frame(chrom=gtf[,'V1'], start=gtf[,'V4']-1, end=gtf[,'V5'])
write.table(gtf, file=snakemake@output["tmp"], quote = F, sep="\t", col.names = F, row.names = F)

system(paste0("gatk BedToIntervalList -I ", snakemake@output["tmp"], " -O ", snakemake@output["intervals"], " -SD ", snakemake@input["refDict"]))