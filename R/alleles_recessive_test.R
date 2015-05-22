

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
        test_enrichment_across_multiple_populations(exac, biallelic_lof, biallelic_func, lof_func, cohort_n)
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
    populations = c("AFR", "EAS", "NFE", "SAS")
    combos = get_count_combinations(populations, biallelic_lof, biallelic_func, lof_func)
    biallelic_lof_combos = combos$biallelic_lof_combos
    biallelic_func_combos = combos$biallelic_func_combos
    lof_func_combos = combos$lof_func_combos
    
    biallelic_lof_p = sum_combo_tests(exac, populations, biallelic_lof_combos, "biallelic_lof_p")
    biallelic_func_p = sum_combo_tests(exac, populations, biallelic_func_combos, "biallelic_func_p")
    lof_func_p = sum_combo_tests(exac, populations, lof_func_combos, "lof_func_p")
    
    p_values = list(lof=NA, functional=NA, biallelic_lof_p=biallelic_lof_p,
        lof_func_p=lof_func_p, biallelic_func_p=biallelic_func_p)
    
    return(p_values)
}

#' get a p-value that sums across different population possibilities
#'
#' @param exac list of frequency estimates for each ExAC population.
#' @param populations a vector of population names to be tested.
#' @param combos a dataframe of the possible count combinations for a functional
#'        type.
#' @param cohort_n list of number of probands in each population.
#' @param p_name the name of the p-value to extract (e.g. "biallelic_lof_p").
#' @export
#'
#' @return a p-value from testing for
sum_combo_tests <- function(exac, cohort_n, populations, combos, p_name) {
    # define a list of counts for each of the functional types. I've arbitrarily
    # given them the names of the resulting p-values, but that is because we
    # pass in the name of the p-value that we want, so we can easily modify the
    # count using that name.
    counts = list(biallelic_lof_p=0, biallelic_func_p=0, lof_func_p=0)
    
    summed_p_value = 0
    for (pos in 1:nrow(combos)) {
        row_p_value = 1
        for (pop in populations) {
            # insert the count for the population into the correct functional
            # count position
            counts[[p_name]] = combos[[pop]][pos]
            
            # for populations with a count of zero, we want to exact probability
            # of the population having zero families. To do this, we adjust the
            # count to 1 (as "test_enrichment" uses count - 1), and later
            # calculate 1 - p-value to get the exact p-value.
            zero = FALSE
            if (counts[[p_name]] == 0) { counts[[p_name]] = 1 ; zero = TRUE }
            
            # define the function arguments
            args = list(freq=exac[[pop]],
                biallelic_lof=counts[["biallelic_lof_p"]],
                biallelic_func=counts[["biallelic_func_p"]],
                lof_func=counts[["lof_func_p"]],
                cohort_n=cohort_n[[pop]])
            
            # call the function using the arguments, then pull out the p-value
            # for the functional type we are looking at.
            p_values = do.call("test_enrichment", args)
            p_value = p_values[[p_name]]
            
            # for zero counts, get 1 - p-value to obtain the exact p-value.
            if (zero) {p_value = 1 - p_value}
            
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
#' @param biallelic_lof number of probands with inherited biallelic LoF variants
#'        in the gene.
#' @param biallelic_func number of probands with inherited biallelic functional
#'        variants in the gene.
#' @param lof_func number of probands with inherited Lof/Func variants in the gene.
#' @export
#'
#' @return a list of count dataframes for each functional type
get_count_combinations <- function(populations, biallelic_lof, biallelic_func, lof_func) {
    types = c("biallelic_lof", "biallelic_func", "lof_func")
    
    # get a matrix of count combinations for each functional type
    biallelic_lof_combos = expand.grid(rep(list(seq(0, biallelic_lof)), each=length(populations)))
    biallelic_func_combos = expand.grid(rep(list(seq(0, biallelic_func)), each=length(populations)))
    lof_func_combos = expand.grid(rep(list(seq(0, (lof_func + biallelic_lof))), each=length(populations)))
    
    # make sure each of the rows sums to the correct value, so that we only use
    # rows where the counts are dispersed correctly amongst the populations
    biallelic_lof_combos = biallelic_lof_combos[rowSums(biallelic_lof_combos) == biallelic_lof, ]
    biallelic_func_combos = biallelic_func_combos[rowSums(biallelic_func_combos) == biallelic_func, ]
    lof_func_combos = lof_func_combos[rowSums(lof_func_combos) == (lof_func + biallelic_lof), ]
    
    # check if any of the rows have all populations with counts greater than one.
    # if the row lacks this, then we need to add a final row where all the
    # counts are one, to cover the scenario
    if (!any(apply(biallelic_lof_combos, 1, function(x) all(x > 0)))) {
        biallelic_lof_combos = rbind(biallelic_lof_combos, rep(1, length(populations)))
    }
    if (!any(apply(biallelic_func_combos, 1, function(x) all(x > 0)))) {
        biallelic_func_combos = rbind(biallelic_func_combos, rep(1, length(populations)))
    }
    if (!any(apply(lof_func_combos, 1, function(x) all(x > 0)))) {
        lof_func_combos = rbind(lof_func_combos, rep(1, length(populations)))
    }
    
    # name the combinations by the population names
    names(biallelic_lof_combos) = populations
    names(biallelic_func_combos) = populations
    names(lof_func_combos) = populations
    
    values = list(biallelic_lof_combos=biallelic_lof_combos,
        biallelic_func_combos=biallelic_func_combos,
        lof_func_combos=lof_func_combos)
    
    return(values)
}
