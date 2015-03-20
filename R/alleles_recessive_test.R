

#' test for enrichment of inherited variants in the DDD and ExAC datasets
#'
#' @param hgnc HGNC symbol for a gene.
#' @param chrom chromosome that the gene is on.
#' @param lof_lof number of probands with inherited biallelic LoF variants in the gene.
#' @param lof_func number of probands with inherited Lof/Func variants in the gene.
#' @param probands vector of probands who have inherited recessive variants in
#'     the gene, or NULL.
#' @param cohort_n number of probands in population.
#' @param check_last_base whether to correct missense or synonymous G alleles at
#'     the last base of exons to a LoF consequence.
#' @export
#'
#' @return a list of P values from tests using the DDD population, the ExAC
#'     population, under LoF and functional tests.
analyse_inherited_enrichment <- function(hgnc, chrom, lof_lof, lof_func, probands=NULL, cohort_n=3072, check_last_base=TRUE) {
    
    cat("extracting ddd frequencies\n")
    ddd = try(get_ddd_variants_for_gene(hgnc, chrom, probands, check_last_base=check_last_base), silent=TRUE)
    if (class(ddd) != "try-error") {
        ddd = get_cumulative_frequencies(ddd)
        ddd = test_enrichment(ddd, lof_lof, lof_func, cohort_n)
    } else {
        ddd=list(lof=NA, func=NA, lof_p=NA, func_p=NA)
    }
    
    cat("extracting ExAC frequencies\n")
    exac = get_exac_variants_for_gene(hgnc, chrom, check_last_base=check_last_base)
    exac = get_cumulative_frequencies(exac)
    exac = test_enrichment(exac, lof_lof, lof_func, cohort_n)
    
    p_values = list(ddd=ddd, exac=exac)
    
    return(p_values)
}

#' test for enrichment of inherited variants
#'
#' @param freq list of cumulative frequencies of variation in populations for
#'     rare LoF variants, and rare functional variants.
#' @param biallelic_lof number of probands with inherited Lof/LoF variants in the gene.
#' @param lof_func number of probands with inherited Lof/Func variants in the gene.
#' @param cohort_n number of probands in population.
#' @export
#'
#' @return a list of P values from tests, under LoF and functional tests.
test_enrichment <- function(freq, biallelic_lof, lof_func, cohort_n) {
    lof_rate = freq$lof * freq$lof
    lof_func_rate = freq$lof * (freq$lof + freq$functional)
    
    # get the probability of getting more than or equal to the number of
    # observed inherited events
    freq$lof_p = pbinom(biallelic_lof - 1, cohort_n, prob=lof_rate, lower.tail=FALSE)
    freq$func_p = pbinom(lof_func + biallelic_lof - 1, cohort_n, prob=lof_func_rate, lower.tail=FALSE)
    
    return(freq)
}
