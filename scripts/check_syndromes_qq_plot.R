
library(qqman)
library(Cairo)

PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.syndrome_tests.permuted.txt"
OUTPUT = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.syndrome_tests.qq.pdf"

p_values = read.table(PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)

Cairo(OUTPUT, type="pdf", height=15, width=15, units="cm")
qq(p_values$syndrome, ylim=c(0, 4), xlim=c(0, 4))
dev.off()
