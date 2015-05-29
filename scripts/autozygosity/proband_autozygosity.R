

library(recessiveStats)
library(argparse)

parser = ArgumentParser()
parser$add_argument("--proband", help="Proband to analyse.")

args = parser$parse_args()

vcf = "/lustre/scratch113/projects/ddd/users/jm33/ddd_4k.bcftools.vcf.gz"
proband_roh = check_sample_autozygosity_genome_wide(vcf, args$proband)
write.table(proband_roh, file=file.path("data-raw", "autozygosity", args$proband), sep="\t", row.names=FALSE, quote=FALSE)
