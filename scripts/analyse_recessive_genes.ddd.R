
library(recessiveStats)
library(jsonlite)

COHORT_N = list("AFR"=109, "EAS"=15, "NFE"=2799, "SAS"=297)
LAST_BASE_RULE = TRUE
RECESSIVE_DIR = "/nfs/users/nfs_j/jm33/apps/recessiveStats"
if (!LAST_BASE_RULE) {
    RECESSIVE_COUNTS_PATH = file.path(RECESSIVE_DIR,"data-raw/recessive_counts_by_gene.txt")
    PROBANDS_PATH = file.path(RECESSIVE_DIR,"data-raw/recessive_probands_by_gene.json")
    OUTPUT_PATH = file.path(RECESSIVE_DIR,"results/recessive.allele_frequency_tests.exac_populations.txt")
} else {
    RECESSIVE_COUNTS_PATH = file.path(RECESSIVE_DIR,"data-raw/recessive_counts_by_gene.last_base_rule.txt")
    PROBANDS_PATH = file.path(RECESSIVE_DIR,"data-raw/recessive_probands_by_gene.last_base_rule.json")
    OUTPUT_PATH = file.path(RECESSIVE_DIR,"results/recessive.allele_frequency_tests.last_base_rule.exac_populations.txt")
}

# load the datasets
recessive_genes = read.table(RECESSIVE_COUNTS_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)
probands = fromJSON(PROBANDS_PATH)

# start a blank dataframe for the results
results = data.frame(gene=character(0), chrom=character(0),
    func_count=character(0), lof_count=character(0),
    ddd.lof=character(0), ddd.functional=character(0),
    ddd.lof_p=character(0), ddd.func_p=character(0),
    exac.lof=character(0), exac.functional=character(0),
    exac.lof_p=character(0), exac.func_p=character(0))

for (gene in sort(unique(recessive_genes$gene))) {
    
    row = recessive_genes[recessive_genes$gene == gene, ]
    proband_ids = probands[[gene]]
    
    chrom = row$chrom
    biallelic_lof = row$lof_sum
    lof_func = row$lof_func
    biallelic_func = 0
    
    result = try(analyse_inherited_enrichment(gene, chrom, biallelic_lof, biallelic_func, lof_func,
        probands=proband_ids, cohort_n=COHORT_N, check_last_base=LAST_BASE_RULE))
    
    if (class(result) == "try-error") { next }
    
    # join the result for the gene to the larger set of gene results
    gene_data = data.frame(hgnc=gene, chrom=chrom, func_count=lof_func,
        lof_count=biallelic_lof)
    result = cbind(gene_data, data.frame(result))
    results = rbind(results, result)
}

write.table(results, file=OUTPUT_PATH, sep="\t", row.names=FALSE, quote=FALSE)
