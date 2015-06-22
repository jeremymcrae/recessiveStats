# the intention is to run this as a LSF job array:
# bsub \
#   -J alleles[1-2685]%100 \
#   -o recessive.bjob_output.%I.txt \
#   -R "select[mem>15000] rusage[mem=15000]" \
#   -M 15000 \
#   bash -c "/software/R-3.1.2/bin/Rscript scripts/analyse_recessive_gene.ddd.R -p \$LSB_JOBINDEX"

library(recessiveStats)
library(jsonlite)

library(getopt)
spec = matrix(c('position', 'p', 2, "integer"), byrow=TRUE, ncol=4);
opt = getopt(spec);

COHORT_N = list("AFR"=100, "EAS"=14, "NFE"=2688, "SAS"=270)
LAST_BASE_RULE = FALSE
RECESSIVE_DIR = "/nfs/users/nfs_j/jm33/apps/recessiveStats"
RECESSIVE_COUNTS_PATH = file.path(RECESSIVE_DIR,"data-raw/recessive_counts_by_gene.silent.txt")
PROBANDS_PATH = file.path(RECESSIVE_DIR,"data-raw/recessive_probands_by_gene.silent.json")

# define all the output paths
RESULTS_DIR = file.path(RECESSIVE_DIR, "results")
RESULTS_PREFIX = "recessive.allele_frequency_tests.silent"
UNCORRECTED_PATH = file.path(RESULTS_DIR, paste(RESULTS_PREFIX, opt$position, "txt", sep="."))
AUTOZYGOSITY_PATH = file.path(RESULTS_DIR, paste(RESULTS_PREFIX, "autozygosity", opt$position, "txt", sep="."))
POP_STRUCTURE_PATH = file.path(RESULTS_DIR, paste(RESULTS_PREFIX, "exac", opt$position, "txt", sep="."))
POP_STRUCTURE_AND_AUTOZYGOSITY_PATH = file.path(RESULTS_DIR, paste(RESULTS_PREFIX, "exac_and_autozygosity", opt$position, "txt", sep="."))

# load the datasets
recessive_genes = read.table(RECESSIVE_COUNTS_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)
probands = fromJSON(PROBANDS_PATH)

autozygous_rates = read.table(file.path(RECESSIVE_DIR, "results/autozygosity.all_genes.all_probands.tsv"),
    sep="\t", header=TRUE, stringsAsFactors=FALSE)

gene = sort(unique(recessive_genes$gene))[opt$position]
    
row = recessive_genes[recessive_genes$gene == gene, ]
proband_ids = probands[[gene]]
autozygosity = autozygous_rates$rate[autozygous_rates$hgnc == gene]
# allow for mising autozygosity rate estimates
if (length(autozygosity) == 0) {
    autozygosity = 0
}

chrom = row$chrom
biallelic_lof = 0
lof_func = 0
biallelic_silent = row$silent_sum
biallelic_func = 0

run_test <- function(cohort_n, autozygosity, output_path) {
# get results for when we are testing the effects of
    result = try(analyse_inherited_enrichment(gene, chrom, biallelic_lof,
        biallelic_func, biallelic_silent, lof_func,
        probands=proband_ids, cohort_n=cohort_n, check_last_base=LAST_BASE_RULE,
        autozygous_rate=autozygosity))
    
    if (class(result) == "try-error") { stop("encountered try-error") }
    
    # join the result for the gene to the larger set of gene results
    gene_data = data.frame(hgnc=gene, chrom=chrom, func_count=lof_func,
        lof_count=biallelic_lof, silent_count=biallelic_silent)
    result = cbind(gene_data, data.frame(result))
    
    write.table(result, file=output_path, sep="\t", row.names=FALSE, quote=FALSE)
}

# run the test for all of the setup permutations
run_test(sum(unlist(COHORT_N)), 0, UNCORRECTED_PATH)
run_test(sum(unlist(COHORT_N)), autozygosity, AUTOZYGOSITY_PATH)
run_test(COHORT_N, 0, POP_STRUCTURE_PATH)
run_test(COHORT_N, autozygosity, POP_STRUCTURE_AND_AUTOZYGOSITY_PATH)
