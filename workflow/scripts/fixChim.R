args <- commandArgs(trailingOnly=TRUE)
tbl <- read.csv(args[[1]], sep="\t")

tbl$chr_donorA <- sapply(tbl$chr_donorA, function(x) {
    if (x %in% c(1:22, "X", "Y")){
        paste0("chr", x)
    } else if (x == "MT") {
        "chrM"
    } else {
        x
    }
})

tbl$chr_acceptorB <- sapply(tbl$chr_acceptorB, function(x) {
    if (x %in% c(1:22, "X", "Y")){
        paste0("chr", x)
    } else if (x == "MT") {
        "chrM"
    } else {
        x
    }
})

write.table(tbl, sep = "\t", na="", file= args[[1]], quote=F, row.names = F)
