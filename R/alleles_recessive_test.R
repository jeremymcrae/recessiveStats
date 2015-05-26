

#' test for enrichment of inherited variants in the DDD and ExAC datasets
#'
#' @param hgnc HGNC symbol for a gene.
#' @param chrom chromosome that the gene is on.
#' @param biallelic_lof number of probands with inherited biallelic LoF variants
#'        in the gene.
#' @param biallelic_func number of probands with inherited biallelic functional
#'        variants in the gene.
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
        exac = test_enrichment_across_multiple_populations(exac, biallelic_lof, biallelic_func, lof_func, cohort_n)
    }
    
    p_values = list(ddd=ddd, exac=exac)
    
    return(p_values)
}

#' test for enrichment of inherited variants
#'
#' @param freq list of cumulative frequencies of variation in a population for
#'        rare LoF variants, and rare functional variants.
#' @param biallelic_lof number of probands with inherited Lof/LoF variants in
#'        the gene.
#' @param biallelic_func number of probands with inherited func/func variants in
#'        the gene.
#' @param lof_func number of probands with inherited Lof/Func variants in the gene.
#' @param cohort_n number of probands in population.
#' @export
#'
#' @return a list of P values from tests, under LoF and functional tests.
test_enrichment <- function(freq, biallelic_lof, biallelic_func, lof_func, cohort_n) {
    
    # get the probability of getting more than or equal to the number of
    # observed inherited events
    freq$biallelic_lof_p = biallelic_lof_enrichment(freq, biallelic_lof, cohort_n)
    freq$lof_func_p = lof_func_enrichment(freq, biallelic_lof + lof_func, cohort_n)
    freq$biallelic_func_p = biallelic_func_enrichment(freq, biallelic_func, cohort_n)
    
    return(freq)
}

#' test for enrichment of biallelic LoF inherited variants
#'
#' @param freq list of cumulative frequencies of variation in a population for
#'        rare LoF variants, and rare functional variants.
#' @param count number of probands with inherited variants in the gene.
#' @param cohort_n number of probands in population.
#' @export
#'
#' @return P-value from testing for biallelic LoF variants.
biallelic_lof_enrichment <- function(freq, count, cohort_n) {
    rate = freq$lof ** 2
    return(pbinom(count - 1, cohort_n, prob=rate, lower.tail=FALSE))
}

#' test enrichment of inherited biallelic LoF and compound heterozygous LoF/Func
#'
#' @param freq list of cumulative frequencies of variation in a population for
#'        rare LoF variants, and rare functional variants.
#' @param count number of probands with inherited variants in the gene.
#' @param cohort_n number of probands in population.
#' @export
#'
#' @return P-value from testing for biallelic LoF and Lof/Func variants.
lof_func_enrichment <- function(freq, count, cohort_n) {
    rate = freq$lof ** 2 + (2 * freq$lof * (1 - freq$lof) * freq$functional)
    return(pbinom(count - 1, cohort_n, prob=rate, lower.tail=FALSE))
}

#' test enrichment of inherited biallelic functional variants
#'
#' @param freq list of cumulative frequencies of variation in a population for
#'        rare LoF variants, and rare functional variants.
#' @param count number of probands with inherited variants in the gene.
#' @param cohort_n number of probands in population.
#' @export
#'
#' @return P-value from testing for biallelic functional variants.
biallelic_func_enrichment <- function(freq, count, cohort_n) {
    rate = freq$functional ** 2
    return(pbinom(count - 1, cohort_n, prob=rate, lower.tail=FALSE))
}

#' test for enrichment of inherited variants in multiple populations
#'
#' @param exac list of frequency estimates for each ExAC population.
#' @param biallelic_lof number of probands with inherited biallelic LoF variants in the gene.
#' @param biallelic_func number of probands with inherited biallelic functional variants in the gene.
#' @param lof_func number of probands with inherited Lof/Func variants in the gene.
#' @param cohort_n list of number of probands in each population.
#' @export
#'
#' @return a list of P values from tests across the populations, under LoF and
#'         functional tests.
test_enrichment_across_multiple_populations <- function(exac, biallelic_lof, biallelic_func, lof_func, cohort_n) {
    # define the ExAC different populations (AFR="African/African American",
    # EAS="East Asian", NFE="Non-Finnish European", SAS="South Asian")
    populations = names(cohort_n)
    biallelic_lof_combos = get_count_combinations(populations, biallelic_lof)
    biallelic_func_combos = get_count_combinations(populations, biallelic_func)
    lof_func_combos = get_count_combinations(populations, biallelic_lof + lof_func)
    
    p_values = list(lof=NA, functional=NA)
    p_values$biallelic_lof_p = sum_combo_tests(exac, cohort_n,
        , biallelic_lof_enrichment)
    p_values$biallelic_func_p = sum_combo_tests(exac, cohort_n,
        biallelic_func_combos, biallelic_func_enrichment)
    p_values$lof_func_p = sum_combo_tests(exac, cohort_n,
        lof_func_combos, lof_func_enrichment)
    
    return(p_values)
}

#' get a p-value that sums across different population possibilities
#'
#' @param exac list of frequency estimates for each ExAC population.
#' @param cohort_n list of number of probands in each population.
#' @param combos a dataframe of the possible count combinations for a functional
#'        type.
#' @param enrich_function function to test enrichment.
#' @export
#'
#' @return a p-value from testing for
sum_combo_tests <- function(exac, cohort_n, combos, enrich_function) {
    summed_p_value = 0
    for (pos in 1:nrow(combos)) {
        row_p_value = 1
        for (pop in names(cohort_n)) {
            # get the functional variants for the population from the row
            count = combos[[pop]][pos]
            
            # for populations where the count is not the highest in the row, we
            # want the probability of the population having that count families.
            # Adjust the count up by 1 (as the enrichment tests use count - 1),
            # and later calculate 1 - p-value to get the p-value for that exact
            # count.
            exact = FALSE
            if (count != max(combos[pos, ])) { count = count + 1 ; exact = TRUE }
            
            p_value = enrich_function(exac[[pop]], count, cohort_n[[pop]])
            
            # for counts that are not the largest in each row, get 1 - p-value
            # to obtain the p-value for that exact count.
            if (exact) {p_value = 1 - p_value}
            
            # the p-value for the row is the product of the p-values for every
            # population
            row_p_value = row_p_value * p_value
        }
        # the overall p-value is the sum of p-values for each row
        summed_p_value = summed_p_value + row_p_value
    }
    
    return(summed_p_value)
}

#' Get all the combinations of spreading the probands across populations.
#'
#' @param populations a vector of population names to be tested
#' @param count number of probands with inherited variants in the gene.
#' @export
#'
#' @return a list of count dataframes for each functional type
get_count_combinations <- function(populations, count) {
    # get a matrix of count combinations
    combos = expand.grid(rep(list(seq(0, count)), each=length(populations)))
    
    n_parity = ceiling(count/length(populations))
    # make sure each of the rows sums to the correct value, so that we only use
    # rows where the counts are dispersed correctly amongst the populations
    combos = combos[rowSums(combos) == count |
        (rowSums(combos) > count & apply(combos, 1, max) <= n_parity), ]
    
    # name the combinations by the population names
    names(combos) = populations
    
    return(combos)
}
