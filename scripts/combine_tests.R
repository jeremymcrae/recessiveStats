# combines all the statistics for recessive testing

library(Cairo)

GENOTYPES_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.last_base_rule.txt"
PHENOTYPES_PATH = "/nfs/users/nfs_j/jm33/apps/hpo_similarity/results/recessive.hpo_similarity.last_base_rule.txt"
# GENOTYPES_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.allele_frequency_tests.txt"
# PHENOTYPES_PATH = "/nfs/users/nfs_j/jm33/apps/hpo_similarity/results/recessive.hpo_similarity.txt"
EXAC_VS_DDD_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/exac_vs_ddd.pdf"
OUTPUT_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.combined_tests.txt"

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

plot_ddd_vs_exac <- function(genotypic_p) {
    # plot the DDD P values versus the ExAC P values, so we can demonstrate the
    # two are very correlated
    max_p = max(c(-log10(genotypic_p$ddd.biallelic_lof_p), -log10(genotypic_p$exac.biallelic_lof_p)))
    Cairo(EXAC_VS_DDD_PATH, type="pdf", height=15, width=15, units="cm")
    plot(-log10(genotypic_p$ddd.biallelic_lof_p), -log10(genotypic_p$exac.biallelic_lof_p), las=1,
        xlab="-log10(P) for DDD LoF", ylab="-log10(P) for ExAC LoF",
        xlim=c(0, max_p), ylim=c(0, max_p))
        
    plot(-log10(genotypic_p$ddd.lof_func_p), -log10(genotypic_p$exac.lof_func_p), las=1,
        xlab="-log10(P) for DDD func", ylab="-log10(P) for ExAC func",
        xlim=c(0, max_p), ylim=c(0, max_p))
            
    dev.off()
}

open_genotype_p_values <- function(path) {
    
    p_values = read.table(path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    # get the maximum P value from the DDD P value and the ExaC P values
    p_values$p_lof = apply(p_values, 1, function(x) max(as.numeric(x[["ddd.biallelic_lof_p"]]), as.numeric(x[["exac.biallelic_lof_p"]])))
    p_values$p_func = apply(p_values, 1, function(x) max(as.numeric(x[["ddd.lof_func_p"]]), as.numeric(x[["exac.lof_func_p"]])))
    
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
    plot_ddd_vs_exac(genotypic)
    
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
