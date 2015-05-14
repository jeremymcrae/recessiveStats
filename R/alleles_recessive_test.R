

#' test for enrichment of inherited variants in the DDD and ExAC datasets
#'
#' @param hgnc HGNC symbol for a gene.
#' @param chrom chromosome that the gene is on.
#' @param biallelic_lof number of probands with inherited biallelic LoF variants in the gene.
#' @param biallelic_func number of probands with inherited biallelic functional variants in the gene.
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
analyse_inherited_enrichment <- function(hgnc, chrom, biallelic_lof, biallelic_func, lof_func, probands=NULL, cohort_n=3072, check_last_base=FALSE) {
    
    cat("extracting ddd frequencies\n")
    ddd = try(get_ddd_variants_for_gene(hgnc, chrom, probands, check_last_base=check_last_base), silent=TRUE)
    if (class(ddd) != "try-error") {
        ddd = get_cumulative_frequencies(ddd)
        ddd = test_enrichment(ddd, biallelic_lof, biallelic_func, lof_func, sum(unlist(cohort_n)))
    } else {
        ddd=list(lof=NA, func=NA, biallelic_lof_p=NA, lof_func_p=NA, biallelic_func_p=NA)
    }
    
    cat("extracting ExAC frequencies\n")
    exac = get_exac_variants_for_gene(hgnc, chrom, check_last_base=check_last_base)
    exac = get_cumulative_frequencies(exac)
    
    if (!is.list(cohort_n)) {
        exac = test_enrichment(exac, biallelic_lof, biallelic_func, lof_func, cohort_n)
    } else {
        
        max_n = biallelic_lof
        
        # get all the possible combinations of spreading the probands amongst
        # the ExAC different populations
        combinations = expand.grid("AFR"=seq(0, max_n), "EAS"=seq(0, max_n),
            "EUR"=seq(0, max_n), "SAS"=seq(0, max_n))
        # select only the combinations where the families are dispersed
        # correctly across the populations (each row should sum to the number of
        # families, no more, no less).
        combinations = combinations[rowSums(combinations) == max_n, ]
        
        for (pos in 1:nrow(combinations)) {
            for (pop in names(combinations)) {
                count = combinations[[pop]][pos]
                frequency = exac[[pop]]
                size = cohort_n[[pop]]
                
                p_values = test_enrichment(frequency, count, 0, 0, size)
            }
        }
    }
    
    p_values = list(ddd=ddd, exac=exac)
    
    return(p_values)
}

#' test for enrichment of inherited variants
#'
#' @param freq list of cumulative frequencies of variation in a population for
#'     rare LoF variants, and rare functional variants.
#' @param biallelic_lof number of probands with inherited Lof/LoF variants in the gene.
#' @param biallelic_func number of probands with inherited func/func variants in the gene.
#' @param lof_func number of probands with inherited Lof/Func variants in the gene.
#' @param cohort_n number of probands in population.
#' @export
#'
#' @return a list of P values from tests, under LoF and functional tests.
test_enrichment <- function(freq, biallelic_lof, biallelic_func, lof_func, cohort_n) {
    lof_rate = freq$lof ** 2
    lof_func_rate = freq$lof ** 2 + (2 * freq$lof * (1 - freq$lof) * freq$functional)
    func_rate = freq$functional ** 2
    
    # get the probability of getting more than or equal to the number of
    # observed inherited events
    freq$biallelic_lof_p = pbinom(biallelic_lof - 1, cohort_n, prob=lof_rate, lower.tail=FALSE)
    freq$lof_func_p = pbinom(lof_func + biallelic_lof - 1, cohort_n, prob=lof_func_rate, lower.tail=FALSE)
    freq$biallelic_func_p = pbinom(biallelic_func - 1, cohort_n, prob=func_rate, lower.tail=FALSE)
    
    return(freq)
}
