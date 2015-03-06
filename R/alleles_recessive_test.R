

#' test for enrichment of inherited variants
#'
#' @param hgnc HGNC symbol for a gene.
#' @param chrom chromosome that the gene is on.
#' @param lof_lof number of probands with inherited Lof/LoF variants in the gene.
#' @param lof_func number of probands with inherited Lof/Func variants in the gene.
#' @param probands vector of probands who have inherited LoF/LoF or LoF/Func
#'     variants in the gene, or NA.
#' @param cohort_n number of probands in population.
#' @param check_last_base whether to correct missense or synonymous G alleles at
#'     the last base of exons to a LoF consequence.
#' @export
#'
#' @return a vector of P values from tests using the DDD population, the ExAC
#'     population, under LoF and functional tests.
analyse_inherited_enrichment <- function(hgnc, chrom, lof_lof, lof_func, probands=NA, cohort_n=3072, check_last_base=TRUE) {
    cat("extracting ddd frequencies\n")
    ddd = get_ddd_variants_for_gene(hgnc, chrom, probands, check_last_base=check_last_base)
    cat("extracting ExAC frequencies\n")
    exac = get_exac_variants_for_gene(hgnc, chrom, check_last_base=check_last_base)
    
    ddd = get_cumulative_frequencies(ddd)
    exac = get_cumulative_frequencies(exac)
    
    ddd$lof_rate = ddd$lof * ddd$lof
    ddd$lof_func_rate = ddd$lof * (ddd$lof + ddd$functional)
    exac$lof_rate = exac$lof * ddd$lof
    exac$lof_func_rate = exac$lof * (exac$lof + exac$functional)
    
    # get the probability of getting more than or equal to the number of
    # observed inherited events
    ddd_lof_p = pbinom(lof_lof - 1, cohort_n, prob=ddd$lof_rate, lower.tail=FALSE)
    exac_lof_p = pbinom(lof_lof - 1, cohort_n, prob=exac$lof_rate, lower.tail=FALSE)
    ddd_func_p = pbinom(lof_func + lof_lof - 1, cohort_n, prob=ddd$lof_func_rate, lower.tail=FALSE)
    exac_func_p = pbinom(lof_func + lof_lof - 1, cohort_n, prob=exac$lof_func_rate, lower.tail=FALSE)
    
    p_values = list(ddd_lof_p=ddd_lof_p, exac_lof_p=exac_lof_p,
        ddd_func_p=ddd_func_p, exac_func_p=exac_func_p, ddd_lof=ddd$lof,
        ddd_func=ddd$func, exac_lof=exac$lof, exac_func=exac$func)
    
    return(p_values)
}
