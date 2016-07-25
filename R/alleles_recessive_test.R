

#' test for enrichment of inherited variants
#'
#' @param counts list of counts of biallelically inherited genotypes in a gene for:
#'             biallelic_lof: inherited Lof/LoF variants.
#'             lof_func: ninherited Lof/Func variants.
#'             biallelic_func: inherited func/func variants.
#' @param variants dataframe of variants in gene
#' @param cohort_n number of probands in population.
#' @param autozygosity rate of autozygosity within the gene in the probands.
#' @export
#'
#' @return a list of P values from LoF and functional tests.
analyse_inherited_enrichment <- function(counts, variants, cohort_n=3072, autozygosity=0) {
    
    freqs = get_cumulative_frequencies(variants)
    
    if (!is.list(cohort_n) || !any(names(freqs) %in% names(cohort_n))) {
        p_values = enrichment_single_population(freqs, counts, sum(unlist(cohort_n)), autozygosity)
    } else {
        stopifnot(all(names(cohort_n) %in% names(freqs)))
        p_values = enrichment_multiple_populations(freqs, counts, cohort_n, autozygosity)
    }
    
    return(p_values)
}

#' test for enrichment of inherited variants
#'
#' @param freq list of cumulative frequencies of variation in a population for
#'        rare LoF variants, and rare functional variants.
#' @param counts list of counts of biallelically inherited genotypes in a gene for:
#'             biallelic_lof: inherited Lof/LoF variants.
#'             lof_func: ninherited Lof/Func variants.
#'             biallelic_func: inherited func/func variants.
#' @param cohort_n number of probands in population.
#' @param autozygosity rate of autozygosity in the cohort being investigated.
#' @export
#'
#' @examples
#' freq = list(lof=0.001, functional=0.1)
#' counts = list("biallelic_lof"=1, 'biallelic_func'=4, 'lof_func'=2)
#' enrichment_single_population(freq, counts, 1000)
#'
#' @return a list of P values from tests, under LoF and functional tests.
enrichment_single_population <- function(freq, counts, cohort_n, autozygosity=0) {
    
    # get the probability of getting more than or equal to the number of
    # observed inherited events
    freq$biallelic_lof_p = biallelic_lof_enrichment(freq, counts$biallelic_lof, cohort_n, autozygosity)
    freq$lof_func_p = lof_func_enrichment(freq, counts$biallelic_lof + counts$lof_func, cohort_n, autozygosity)
    freq$biallelic_func_p = biallelic_func_enrichment(freq, counts$biallelic_func, cohort_n, autozygosity)
    
    return(freq)
}

#' test for enrichment of biallelic LoF inherited variants
#'
#' @param freq list of cumulative frequencies of variation in a population for
#'        rare LoF variants, and rare functional variants.
#' @param count number of probands with inherited variants in the gene.
#' @param cohort_n number of probands in population.
#' @param autozygosity rate of autozygosity in the cohort being investigated.
#' @export
#'
#' @examples
#' freq = list(lof=0.001, functional=0.1)
#' count = 3
#' biallelic_lof_enrichment(freq, count, 1000)
#'
#' @return P-value from testing for biallelic LoF variants.
biallelic_lof_enrichment <- function(freq, count, cohort_n, autozygosity=0) {
    rate = (freq$lof ** 2) * (1 - autozygosity) + freq$lof * autozygosity
    return(stats::pbinom(count - 1, cohort_n, prob=rate, lower.tail=FALSE))
}

#' test enrichment of inherited biallelic LoF and compound heterozygous LoF/Func
#'
#' @param freq list of cumulative frequencies of variation in a population for
#'        rare LoF variants, and rare functional variants.
#' @param count number of probands with inherited variants in the gene.
#' @param cohort_n number of probands in population.
#' @param autozygosity rate of autozygosity in the cohort being investigated.
#' @export
#'
#' @examples
#' freq = list(lof=0.001, functional=0.1)
#' count = 3
#' lof_func_enrichment(freq, count, 1000)
#'
#' @return P-value from testing for biallelic LoF and Lof/Func variants.
lof_func_enrichment <- function(freq, count, cohort_n, autozygosity=0) {
    rate = (freq$lof ** 2) * (1 - autozygosity) +
        (freq$lof) * autozygosity +
        (2 * freq$lof * (1 - freq$lof) * freq$functional)
    return(stats::pbinom(count - 1, cohort_n, prob=rate, lower.tail=FALSE))
}

#' test enrichment of inherited biallelic functional variants
#'
#' @param freq list of cumulative frequencies of variation in a population for
#'        rare LoF variants, and rare functional variants.
#' @param count number of probands with inherited variants in the gene.
#' @param cohort_n number of probands in population.
#' @param autozygosity rate of autozygosity in the cohort being investigated.
#' @export
#'
#' @examples
#' freq = list(lof=0.001, functional=0.1)
#' count = 3
#' biallelic_func_enrichment(freq, count, 1000)
#'
#' @return P-value from testing for biallelic functional variants.
biallelic_func_enrichment <- function(freq, count, cohort_n, autozygosity=0) {
    rate = (freq$functional ** 2) * (1 - autozygosity) + freq$functional * autozygosity
    return(stats::pbinom(count - 1, cohort_n, prob=rate, lower.tail=FALSE))
}

#' test for enrichment of inherited variants in multiple populations
#'
#' @param freqs list of frequency estimates for each ExAC population.
#' @param counts list of counts of biallelically inherited genotypes in a gene for:
#'             biallelic_lof: inherited Lof/LoF variants.
#'             lof_func: ninherited Lof/Func variants.
#'             biallelic_func: inherited func/func variants.
#' @param cohort_n list of number of probands in each population.
#' @param autozygosity rate of autozygosity in the cohort being investigated.
#' @export
#'
#' @examples
#' freqs = list("AFR"=list("lof"=0.01, "functional"=0.1),
#'     "EAS"=list("lof"=0.02, "functional"=0.1))
#' counts = list("biallelic_lof"=4, 'biallelic_func'=6, 'lof_func'=3)
#' cohort_n = list("AFR"=50, "EAS"=150)
#' autozygous_rate = 0.005
#' enrichment_multiple_populations(freqs, counts, cohort_n, autozygous_rate)
#'
#' @return a list of P values from tests across the populations, under LoF and
#'         functional tests.
enrichment_multiple_populations <- function(freqs, counts, cohort_n, autozygosity=0) {
    # define the ExAC different populations (AFR="African/African American",
    # EAS="East Asian", NFE="Non-Finnish European", SAS="South Asian")
    populations = names(cohort_n)
    biallelic_lof_combos = get_count_combinations(populations, counts$biallelic_lof)
    biallelic_func_combos = get_count_combinations(populations, counts$biallelic_func)
    lof_func_combos = get_count_combinations(populations, counts$biallelic_lof + counts$lof_func)
    
    p_values = list(lof=NA, functional=NA)
    p_values$biallelic_lof_p = sum_combo_tests(freqs, cohort_n,
        biallelic_lof_combos, biallelic_lof_enrichment, autozygosity)
    p_values$lof_func_p = sum_combo_tests(freqs, cohort_n,
        lof_func_combos, lof_func_enrichment, autozygosity)
    p_values$biallelic_func_p = sum_combo_tests(freqs, cohort_n,
        biallelic_func_combos, biallelic_func_enrichment, autozygosity)
    
    return(p_values)
}

#' get a p-value that sums across different population possibilities
#'
#' @param freqs list of frequency estimates for each ExAC population.
#' @param cohort_n list of number of probands in each population.
#' @param combos a dataframe of the possible count combinations for a functional
#'        type.
#' @param enrich_function function to test enrichment.
#' @param autozygosity rate of autozygosity in the cohort being investigated.
#' @export
#'
#' @examples
#' freqs = list("A"=list("lof"=0.1, "functional"=0.1),
#'     "B"=list("lof"=0.1, "functional"=0.1))
#' cohort_n = list("A"=100, "B"=100)
#' combos = get_count_combinations(names(cohort_n), count=1)
#' sum_combo_tests(freqs, cohort_n, combos, biallelic_lof_enrichment)
#'
#' @return a p-value from testing for
sum_combo_tests <- function(freqs, cohort_n, combos, enrich_function, autozygosity=0) {
    if (length(combos) == 1 && is.na(combos)) { return(NA) }
    
    summed_p_value = 0
    for (pos in 1:nrow(combos)) {
        row_p_value = 1
        for (pop in names(cohort_n)) {
            # get the functional variants for the population from the row
            count = combos[[pop]][pos]
            
            # For populations where the count is not the highest in the row, we
            # want the probability of the population having that count families.
            # Otherwise we want the probability of having that count or greater.
            if (count == max(combos[pos, ])) {
                p_value = enrich_function(freqs[[pop]], count, cohort_n[[pop]], autozygosity)
            } else {
                inclusive = enrich_function(freqs[[pop]], count, cohort_n[[pop]], autozygosity)
                exclusive = enrich_function(freqs[[pop]], count + 1, cohort_n[[pop]], autozygosity)
                p_value = inclusive - exclusive
            }
            
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
#' @examples
#' populations = c("AFR", "NFE", "SAS")
#' count = 4
#' get_count_combinations(populations, count)
#'
#' @return a list of count dataframes for each functional type
get_count_combinations <- function(populations, count) {
    
    if (is.null(count) || is.na(count)) { return(NA)}
    
    stopifnot(count >= 0)
    
    # get a matrix of count combinations
    combos = expand.grid(rep(list(seq(0, count)), each=length(populations)))
    
    n_parity = ceiling(count/length(populations))
    # make sure each of the rows sums to the correct value, so that we only use
    # rows where the counts are dispersed correctly amongst the populations
    combos = combos[rowSums(combos) == count |
        (rowSums(combos) > count & apply(combos, 1, max) <= n_parity), ]
    
    # tidy up the dataframe so that it is of a standard format
    combos = data.frame(combos)
    row.names(combos) = 1:nrow(combos)
    
    # name the combinations by the population names
    names(combos) = populations
    
    return(combos)
}
