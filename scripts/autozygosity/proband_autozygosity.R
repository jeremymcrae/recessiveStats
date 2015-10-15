

library(recessiveStats)
library(argparse)

main <- function() {
    parser = ArgumentParser()
    parser$add_argument("--proband", help="Proband to analyse.")
    parser$add_argument("--bcf", help="Path to bcf to analyse.")

    args = parser$parse_args()

    proband_roh = check_sample_autozygosity_genome_wide(args$vcf, args$proband)
    write.table(proband_roh,
        file=file.path("data-raw", "autozygosity", args$proband),
        sep="\t", row.names=FALSE, quote=FALSE)
}

main()
