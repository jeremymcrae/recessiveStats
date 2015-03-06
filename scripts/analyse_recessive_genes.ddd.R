
library(recessiveStats)

RECESSIVE_COUNTS_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/data-raw/recessive_counts_per_gene.last_base_rule.txt"
PROBANDS_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/data-raw/recessive_probands_per_gene.last_base_rule.txt"
OUTPUT_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive_frequency_testing.last_base_rule.txt"
COHORT_N = 3072

# load the datasets
recessive_genes = read.table(RECESSIVE_COUNTS_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)
probands = read.table(PROBANDS_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# start a blank dataframe for the results
results = data.frame(gene=character(0), chrom=character(0),
    func_count=character(0), lof_count=character(0), ddd_lof_p=character(0),
    exac_lof_p=character(0), ddd_func_p=character(0),  exac_func_p=character(0),
    ddd_lof=character(0), ddd_func=character(0), exac_lof=character(0),
    exac_func=character(0))

for (gene in sort(unique(recessive_genes$gene))) {
    
    row = recessive_genes[recessive_genes$gene == gene, ]
    proband_ids = probands$proband[probands$gene == gene]
    
    chrom = row$chrom
    lof_lof = row$lof_sum
    lof_func = row$lof_func
    
    result = try( analyse_inherited_enrichment(gene, chrom, lof_lof, lof_func, probands=proband_ids, cohort_n=COHORT_N, check_last_base=TRUE))
    
    if (class(result) == "try-error") { next }
    
    # join the result for the gene to the larger set of gene results
    gene_data = data.frame(hgnc=gene, chrom=chrom, func_count=lof_func,
        lof_count=lof_lof)
    result = cbind(gene_data, data.frame(result))
    results = rbind(results, result)
}

write.table(results, file=OUTPUT_PATH, sep="\t", row.names=FALSE, quote=FALSE)
