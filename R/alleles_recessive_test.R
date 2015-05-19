

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
#' @param exac list of frequency estimates for each ExAC population
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
    
    # Since we check all combinations of the different functional counts, we
    # we end up with duplicates. Find out which rows are duplicates, so we
    # can exclude them later.
    lof_dups = duplicated(biallelic_lof_combos)
    func_dups = duplicated(biallelic_func_combos)
    lof_func_dups = duplicated(cbind(biallelic_lof_combos, lof_func_combos))
    
    # run through all the possible combinations of functional counts, and
    # get a p-value for each population in each row.
    all_p_values = list()
    for (pos in 1:nrow(biallelic_lof_combos)) {
        for (pop in populations) {
            lof_count = biallelic_lof_combos[[pop]][pos]
            func_count = biallelic_func_combos[[pop]][pos]
            lof_func_count = lof_func_combos[[pop]][pos]
            frequency = exac[[pop]]
            size = cohort_n[[pop]]
            
            p_values = test_enrichment(frequency, lof_count, func_count, lof_func_count, size)
            p_values$lof_n = lof_count
            p_values$func_count = func_count
            p_values$lof_func_n = lof_func_count
            p_values$cohort_n = size
            all_p_values[[pop]][[pos]] = p_values
        }
    }
    
    return(p_values)
}

#' Get all the combinations of spreading the probands across populations.
#'
#' @param populations a vector of population names to be tested
#' @param biallelic_lof number of probands with inherited biallelic LoF variants
#'       in the gene.
#' @param biallelic_func number of probands with inherited biallelic functional
#'        variants in the gene.
#' @param lof_func number of probands with inherited Lof/Func variants in the gene.
#' @export
#'
#' @return a list of count dataframes for each functional type
get_count_combinations <- function(populations, biallelic_lof, biallelic_func, lof_func) {
    types = c("biallelic_lof", "biallelic_func", "lof_func")
    combinations = expand.grid(rep(lapply(c(biallelic_lof, biallelic_func,
        lof_func), function(x) seq(0, x)), each=length(populations)))
    names(combinations) = paste(rep(populations, length(types)),
        sort(rep(types, length(populations))), sep=".")
    
    # split the combinations out by functional type, and rename the columns
    # to the corresponding population
    biallelic_lof_combos = combinations[, grepl("biallelic_lof", names(combinations))]
    biallelic_func_combos = combinations[, grepl("biallelic_func", names(combinations))]
    lof_func_combos = combinations[, grepl("lof_func", names(combinations))]
    names(biallelic_lof_combos) = populations
    names(biallelic_func_combos) = populations
    names(lof_func_combos) = populations
    
    # select only the combinations where the families are dispersed
    # correctly across the populations (each row should sum to the number of
    # families, no more, no less).
    use_rows = rowSums(biallelic_func_combos) == biallelic_func &
        rowSums(biallelic_lof_combos) == biallelic_lof &
        rowSums(lof_func_combos) == lof_func
    
    # drop the functional groups to the rows that are correct across all
    # groups.
    biallelic_lof_combos = biallelic_lof_combos[use_rows, ]
    biallelic_func_combos = biallelic_func_combos[use_rows, ]
    lof_func_combos = lof_func_combos[use_rows, ]
    
    values = list(biallelic_lof_combos=biallelic_lof_combos,
        biallelic_func_combos=biallelic_func_combos,
        lof_func_combos=lof_func_combos)
    
    return(values)
}
