
library(Cairo)
library(qqman)

SILENT_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.silent.txt"

silent = read.table(SILENT_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)

silent$max_p = apply(silent[, c("ddd.biallelic_silent_p", "exac.biallelic_silent_p")], 1, max)

Cairo(file="recessive.synonymous_enrichment.qqplot.pdf", type="pdf", height=15, width=15, units="cm")
qq(silent$max_p)
dev.off()
