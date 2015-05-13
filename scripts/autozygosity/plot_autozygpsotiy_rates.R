
library(Cairo)

ALL_PROBANDS_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/autozygosity.all_genes.all_probands.bjob_output.txt"
CONSANG_PROBANDS_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/autozygosity.all_genes.consanguinous_probands.bjob_output.txt"

get_autozygosity <- function(path) {
    # read the autozygosity output (which is a bjob output file)
    probands = scan(path, what="", sep="\n")
    probands = probands[grepl("^autozygosity_rate", probands)]

    # convert the data to a dataframe
    autozygous = read.table(text=probands, sep="\t", header=FALSE, fill=TRUE)
    names(autozygous) = c("temp", "hgnc", "chrom", "rate")
    autozygous$rate = gsub("warning messages from top-level task callback 1", "", autozygous$rate)
    autozygous$rate = gsub("Loading required package: proto", "", autozygous$rate)
    autozygous$rate = as.numeric(as.character(autozygous$rate))
    
    # select the necessary columns, and order by gene synmbol
    autozygous = autozygous[, c("hgnc", "chrom", "rate")]
    autozygous = autozygous[order(autozygous$hgnc), ]
    
    return(autozygous)
}

plot_autozygosity <- function(all_probands, consang) {
    # exclude chrX, since males always appear autozygous, and I hadn't removed
    # them for this analysis run
    all_probands = all_probands[!all_probands$chrom == "X", ]
    all_density = density(all_probands$rate, na.rm=TRUE)
    consang = consang[!consang$chrom == "X", ]
    consang_density = density(consang$rate, na.rm=TRUE)

    max_y = max(all_density$y, consang_density$y)
    max_x = max(all_density$x, consang_density$x)
    min_x = min(all_density$x, consang_density$x)

    # plot the autozygosity rates
    Cairo(file="autozygosity_rates.pdf", type="pdf", height=15, width=15, units="cm")
    plot(all_density,
        xlab="Autozygosity rate in DDD probands",
        main="autozygosity rates for autosomal genes",
        xlim=c(min_x, max_x), ylim=c(0, max_y), las=1,
        cex=1.3, cex.lab=1.3, cex.axis=1.3, tck=-0.03)
    lines(consang_density, col="red")
    legend("topright", legend=c("all probands", "consanguinous probands"),
        col=c("black", "red"), bty="n", lwd=1)
    dev.off()

}

all_autozygous = get_autozygosity(ALL_PROBANDS_PATH)
consang_autozygous = get_autozygosity(CONSANG_PROBANDS_PATH)
write.table(all_autozygous, file="autozygosity.all_genes.all_probands.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(consang_autozygous, file="autozygosity.all_genes.consanguinous_probands.tsv", sep="\t", row.names=FALSE, quote=FALSE)

plot_autozygosity(all_autozygous, consang_autozygous)
