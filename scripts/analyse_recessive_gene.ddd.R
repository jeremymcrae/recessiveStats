# the intention is to run this as a LSF job array:
# bsub \
#   -J alleles[1-77]%100 \
#   -o recessive.bjob_output.%I.txt \
#   -R "select[mem>15000] rusage[mem=15000]" \
#   -M 15000 \
#   bash -c "/software/R-3.1.2/bin/Rscript scripts/analyse_recessive_gene.ddd.R -p \$LSB_JOBINDEX"

library(recessiveStats)
library(jsonlite)

library(getopt)
spec = matrix(c('position', 'p', 2, "integer"), byrow=TRUE, ncol=4);
opt = getopt(spec);

COHORT_N = list("AFR"=109, "EAS"=15, "NFE"=2799, "SAS"=297)
LAST_BASE_RULE = TRUE
RECESSIVE_DIR = "/nfs/users/nfs_j/jm33/apps/recessiveStats"
if (!LAST_BASE_RULE) {
    RECESSIVE_COUNTS_PATH = file.path(RECESSIVE_DIR,"data-raw/recessive_counts_by_gene.txt")
    PROBANDS_PATH = file.path(RECESSIVE_DIR,"data-raw/recessive_probands_by_gene.json")
    OUTPUT_PATH = file.path(RECESSIVE_DIR,paste("results/recessive.allele_frequency_tests.", opt$position, ".txt", sep=""))
} else {
    RECESSIVE_COUNTS_PATH = file.path(RECESSIVE_DIR,"data-raw/recessive_counts_by_gene.last_base_rule.txt")
    PROBANDS_PATH = file.path(RECESSIVE_DIR,"data-raw/recessive_probands_by_gene.last_base_rule.json")
    OUTPUT_PATH = file.path(RECESSIVE_DIR,paste("results/recessive.allele_frequency_tests.last_base_rule.", opt$position, ".txt", sep=""))
}

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
if (length(autozygosity) == 0) { autozygosity = 0 }

chrom = row$chrom
biallelic_lof = row$lof_sum
lof_func = row$lof_func
biallelic_func = 0

result = try(analyse_inherited_enrichment(gene, chrom, biallelic_lof, biallelic_func, lof_func,
    probands=proband_ids, cohort_n=COHORT_N, check_last_base=LAST_BASE_RULE,
    autozygous_rate=autozygosity))

if (class(result) == "try-error") { }

# join the result for the gene to the larger set of gene results
gene_data = data.frame(hgnc=gene, chrom=chrom, func_count=lof_func,
    lof_count=biallelic_lof)
result = cbind(gene_data, data.frame(result))

write.table(result, file=OUTPUT_PATH, sep="\t", row.names=FALSE, quote=FALSE)
