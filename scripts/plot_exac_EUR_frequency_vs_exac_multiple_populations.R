### quick script to check the impact of updating the ExAC allele frequencies to
### take different DDD ethnicities into account.

library(Cairo)
library(qqman)

# define the paths to the initial and updated test results
INITIAL_FREQUENCIES_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.adjusted_AC.txt"
MULTIPLE_ETHNICITIES_FREQUENCIES_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.exac_populations.txt"
AUTOZYGOSITY_FREQUENCIES_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.autozygosity.txt"
MULTIPLE_ETHNICITIES_AND_AUTOZYGOSITY_FREQUENCIES_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.exac_and_autozygosity.txt"

# define the paths to the tests results using rare synonymous variants
INITIAL_SYNONYMOUS_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.silent.txt"
MULTIPLE_ETHNICITIES_SYNONYMOUS_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.silent.exac.txt"
AUTOZYGOSITY_SYNONYMOUS_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.silent.autozygosity.txt"
MULTIPLE_ETHNICITIES_AUTOZYGOSITY_SYNONYMOUS_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.silent.exac_and_autozygosity.txt"

# load all the subsets from the modified tests
initial = read.table(INITIAL_FREQUENCIES_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)
multiple_ethnicities = read.table(MULTIPLE_ETHNICITIES_FREQUENCIES_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)
autozygosity = read.table(AUTOZYGOSITY_FREQUENCIES_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)
multiple_ethnicities_and_autozygosity = read.table(MULTIPLE_ETHNICITIES_AND_AUTOZYGOSITY_FREQUENCIES_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# load the results for the synonymous tests
initial_silent = read.table(INITIAL_SYNONYMOUS_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)
multiple_ethnicities_silent = read.table(MULTIPLE_ETHNICITIES_SYNONYMOUS_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)
autozygosity_silent = read.table(AUTOZYGOSITY_SYNONYMOUS_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)
multiple_ethnicities_and_autozygosity_silent = read.table(MULTIPLE_ETHNICITIES_AUTOZYGOSITY_SYNONYMOUS_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)

plot_frequency_comparison <-function(initial, updated, p_name, title) {
    # make sure we can distingush the updated p-values from the initial p-values
    names(initial) = paste("initial", names(initial), sep=".")
    names(updated) = paste("updated", names(updated), sep=".")
    
    frequencies = merge(initial, updated, by.x="initial.hgnc", by.y="updated.hgnc")
    
    # check the correlation between P-values from the initial and updated tests
    r_squared = cor(-log10(frequencies[[paste("initial.exac.", p_name, sep="")]]),
        -log10(frequencies[[paste("updated.exac.", p_name, sep="")]]),
        use="pairwise.complete.obs") ** 2
    
    max_val = max(c(-log10(frequencies[[paste("initial.exac.", p_name, sep="")]])),
        -log10(frequencies[[paste("updated.exac.", p_name, sep="")]]),
        na.rm=TRUE)
    
    # plot the comparison of P-values between the old and updated tests
    plot(-log10(frequencies[[paste("initial.exac.", p_name, sep="")]]),
        -log10(frequencies[[paste("updated.exac.", p_name, sep="")]]),
        main=title,
        xlab=paste("-log10(", p_name, "(old))"),
        ylab=paste("-log10(", p_name, "(new))"),
        xlim=c(0, max_val), ylim=c(0, max_val), las=1)
    abline(0, 1)
    legend("bottomright", legend=c("y=x"), col=c("black"), lty=1, bty="n")
    legend("topleft", legend=paste("R^2 = ", signif(r_squared, 3), sep=""), bty="n")
    # dev.off()
}

Cairo(file="updated_recessive_test_comparison.pdf", type="pdf", height=15, width=15, units="cm")
plot_frequency_comparison(initial, multiple_ethnicities, "biallelic_lof_p", "initial vs multiple ethnicities")
plot_frequency_comparison(initial, multiple_ethnicities, "lof_func_p", "initial vs multiple ethnicities")
plot_frequency_comparison(initial, autozygosity, "biallelic_lof_p", "initial vs autozygosity")
plot_frequency_comparison(initial, multiple_ethnicities_and_autozygosity, "biallelic_lof_p", "initial vs multiple ethnicities and autozygosity")
plot_frequency_comparison(initial_silent, multiple_ethnicities_silent, "biallelic_silent_p", "initial synonymous vs multiple ethnicities synonymous")
plot_frequency_comparison(initial_silent, autozygosity_silent, "biallelic_silent_p", "initial synonymous vs autozygosity synonymous")
plot_frequency_comparison(initial_silent, multiple_ethnicities_and_autozygosity_silent, "biallelic_silent_p", "initial synonymous vs multiple ethnicities and autozygosity synonymous")
dev.off()

# and make QQ plots of the results from testing the synonymous variants
Cairo(file="initial_vs_updated_tests.qq_plot.pdf", type="pdf", height=15, width=15, units="cm")
qq(initial_silent$exac.biallelic_silent_p, las=1, xlim=c(0, 5), ylim=c(0,5))
par(new=TRUE)
qq(multiple_ethnicities_silent$exac.biallelic_silent_p,
    xlim=c(0, 5), ylim=c(0,5), col="red", axes=FALSE, xlab=NA, ylab=NA)
par(new=TRUE)
qq(autozygosity_silent$exac.biallelic_silent_p,
    xlim=c(0, 5), ylim=c(0,5), col="green", axes=FALSE, xlab=NA, ylab=NA)
par(new=TRUE)
qq(multiple_ethnicities_and_autozygosity_silent$exac.biallelic_silent_p,
    xlim=c(0, 5), ylim=c(0,5), col="gray", axes=FALSE, xlab=NA, ylab=NA)
legend("topleft",
    legend=c("initial",
        "adjusted for multiple ethnicities",
        "adjusted for autozygosity rates",
        "adjusted for multiple ethnicities and autozygosity"),
    col=c("black", "red", "green", "gray", "brown", "gold"), pch=19, bty="n", cex=0.8)
dev.off()
