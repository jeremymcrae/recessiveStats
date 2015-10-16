

library(recessiveStats)
library(argparse)

main <- function() {
    parser = ArgumentParser()
    parser$add_argument("--proband", help="Proband to analyse.")
    parser$add_argument("--bcf", help="Path to bcf to analyse.")
    parser$add_argument("--output", help="Path to write output to.")

    args = parser$parse_args()

    proband_roh = check_sample_autozygosity_genome_wide(args$bcf, args$proband)
    write.table(proband_roh, file=args$output, sep="\t", row.names=FALSE, quote=FALSE)
}

main()
