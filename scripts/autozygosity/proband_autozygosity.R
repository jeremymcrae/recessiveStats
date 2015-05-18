

library(recessiveStats)
library(argparse)

parser = ArgumentParser()
parser$add_argument("--proband", help="HGNC symbol for gene to analyse.")

args = parser$parse_args()

bcf = "/lustre/scratch113/projects/ddd/users/jm33/ddd_4k.bcf"
proband = check_sample_autozygosity_genome_wide(bcf, args$proband)
write.table(proband, file=file.path("data-raw", args$proband) sep="\t", row.names=FALSE, quote=FALSE)
