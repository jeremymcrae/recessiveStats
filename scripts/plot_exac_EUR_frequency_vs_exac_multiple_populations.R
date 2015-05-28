### quick script to check the impact of updating the ExAC allele frequencies to
### take different DDD ethnicities into account.

library(Cairo)

INITIAL_FREQUENCIES_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.adjusted_AC.txt"
UPDATED_FREQUENCIES_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.exac_populations.txt"


initial = read.table(INITIAL_FREQUENCIES_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)
updated = read.table(UPDATED_FREQUENCIES_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)

names(initial) = paste("initial", names(initial), sep=".")
names(updated) = paste("updated", names(updated), sep=".")


frequencies = merge(initial, updated, by.x="initial.hgnc", by.y="updated.hgnc", all.x=TRUE)

# check the correlation between P-values from the
cor(-log10(frequencies$initial.exac.biallelic_lof_p),
    -log10(frequencies$updated.exac.biallelic_lof_p)) ** 2

cor(-log10(frequencies$initial.exac.lof_func_p),
    -log10(frequencies$updated.exac.lof_func_p)) ** 2

# plot the comparison of P-values between the old and updated tests
Cairo(file="updated_recessive_test_comparison.pdf", type="pdf", height=15, width=30, units="cm")
par(mfrow = c(1,2))
plot(-log10(frequencies$initial.exac.biallelic_lof_p),
    -log10(frequencies$updated.exac.biallelic_lof_p),
    xlab="-log10(biallelic Lof P (old test))", ylab="-log10(biallelic Lof P (new test))", las=1)
abline(0, 1)

plot(-log10(frequencies$initial.exac.lof_func_p),
    -log10(frequencies$updated.exac.lof_func_p),
    xlab="-log10(LoF/func P (old test))", ylab="-log10(LoF/func P (new test))", las=1)
abline(0, 1)
legend("bottomright", legend=c("y=x"), col=c("black"), lty=1)
dev.off()
