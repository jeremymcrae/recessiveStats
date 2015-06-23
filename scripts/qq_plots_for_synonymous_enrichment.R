### quick script to check the impact of updating the ExAC allele frequencies to
### take different DDD ethnicities into account.

library(Cairo)
library(qqman)

# define the paths to the tests results using rare synonymous variants
INITIAL_SYNONYMOUS_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.silent.txt"
MULTIPLE_ETHNICITIES_SYNONYMOUS_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.silent.exac.txt"
AUTOZYGOSITY_SYNONYMOUS_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.silent.autozygosity.txt"
MULTIPLE_ETHNICITIES_AUTOZYGOSITY_SYNONYMOUS_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.silent.exac_and_autozygosity.txt"

# load the results for the synonymous tests
initial = read.table(INITIAL_SYNONYMOUS_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)
multiple_ethnicities = read.table(MULTIPLE_ETHNICITIES_SYNONYMOUS_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)
autozygosity = read.table(AUTOZYGOSITY_SYNONYMOUS_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)
multiple_ethnicities_and_autozygosity = read.table(MULTIPLE_ETHNICITIES_AUTOZYGOSITY_SYNONYMOUS_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# make QQ plots of the results from testing the synonymous variants
Cairo(file="initial_vs_updated_tests.qq_plot.pdf", type="pdf", height=15, width=15, units="cm")
qq(initial$exac.biallelic_p, las=1, xlim=c(0, 5), ylim=c(0,5))
par(new=TRUE)
qq(multiple_ethnicities$exac.biallelic_p,
    xlim=c(0, 5), ylim=c(0,5), col="red", axes=FALSE, xlab=NA, ylab=NA)
par(new=TRUE)
qq(autozygosity$exac.biallelic_p,
    xlim=c(0, 5), ylim=c(0,5), col="green", axes=FALSE, xlab=NA, ylab=NA)
par(new=TRUE)
qq(multiple_ethnicities_and_autozygosity$exac.biallelic_p,
    xlim=c(0, 5), ylim=c(0,5), col="gray", axes=FALSE, xlab=NA, ylab=NA)
legend("topleft",
    legend=c("initial",
        "adjusted for multiple ethnicities",
        "adjusted for autozygosity rates",
        "adjusted for multiple ethnicities and autozygosity"),
    col=c("black", "red", "green", "gray"), pch=19, bty="n", cex=0.8)
dev.off()
