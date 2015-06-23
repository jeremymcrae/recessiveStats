# combines all the statistics for recessive testing

library(Cairo)

GENOTYPES_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.last_base_rule.txt"
PHENOTYPES_PATH = "/nfs/users/nfs_j/jm33/apps/hpo_similarity/results/recessive.hpo_similarity.last_base_rule.txt"
# GENOTYPES_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.txt"
# PHENOTYPES_PATH = "/nfs/users/nfs_j/jm33/apps/hpo_similarity/results/recessive.hpo_similarity.txt"
EXAC_VS_DDD_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/exac_vs_ddd.pdf"
OUTPUT_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.combined_tests.ver3.txt"

#' function to combine p values, using Fisher's method
#'
#' @param x vector of P values for a gene
#' @export
#'
#' @return single P value for gene
fishersMethod <- function(x) {
    x = x[!is.na(x)]
    if (length(x) == 0) { return(NA) }
    
    adjusted = pchisq(-2 * sum(log(x)), df=2 * length(x), lower.tail=FALSE)
    
    return(adjusted)
}

open_genotype_p_values <- function(path) {
    
    p_values = read.table(path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    # get the ExAC P values
    p_values$p_lof = p_values$exac.biallelic_lof_p
    p_values$p_func =  p_values$exac.lof_func_p
    
    # select the most powerful test for each gene
    p_values$p_genotype = apply(p_values, 1, function(x) min(as.numeric(x[["p_lof"]]), as.numeric(x[["p_func"]])))
    
    return(p_values)
}

open_phenotype_p_values <- function(path) {
    p_values = read.table(path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    p_values$p_phenotype = p_values$hpo_similarity_p_value
    p_values$hpo_similarity_p_value = NULL
    
    return(p_values)
}

main <- function() {
    genotypic = open_genotype_p_values(GENOTYPES_PATH)
    phenotypic = open_phenotype_p_values(PHENOTYPES_PATH)
    
    all = merge(genotypic, phenotypic, by="hgnc", all=TRUE)
    all$p_combined = apply(all, 1, function(x) fishersMethod(c(as.numeric(x[["p_genotype"]]), as.numeric(x[["p_phenotype"]]))))
    
    # exclude a few genes with known problems
    all = all[!(all$hgnc %in% c("CPSF1", "LINC01551", "WWOX")), ]
    
    # sort the genes by significant, round the values etc
    all = all[order(all$p_combined), ]
    all = all[, c("hgnc", "chrom", "func_count", "lof_count", "p_genotype", "p_phenotype", "p_combined")]
    # all$p_genotype = signif(all$p_genotype, 3)
    # all$p_phenotype = signif(all$p_phenotype, 3)
    # all$p_combined = signif(all$p_combined, 3)
    
    write.table(all, file=OUTPUT_PATH, sep="\t", row.names=FALSE, quote=FALSE)
    
}

main()
