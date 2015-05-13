

library(recessiveStats)
library(argparse)

parser = ArgumentParser()
parser$add_argument("--hgnc", help="HGNC symbol for gene to analyse.")
parser$add_argument("--chrom", help="chromosome that the gene is on.")
parser$add_argument("--consang", default=FALSE, action="store_true", help="whether to restrict to the consanguinous probands.")

args = parser$parse_args()

rate = get_autozygous_rate(args$hgnc, args$chrom, args$consang)
cat(paste("autozygosity_rate", args$hgnc, args$chrom, rate, "RATE_END", sep="\t"))
