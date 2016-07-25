
#' obtain cumulative allele frequency of rare lof and functional variants
#'
#' @param vars dataframe of variants (one row per allele), which includes the
#'        number of times that allele was observed within the population, as
#'        well as the total number of alleles in the population. Alternatively,
#'        this can be a list of dataframe, each for a different population
#         (e.g. list("EAS"=df(...), "SAS"=df(...))).
#' @export
#'
#' @return a list of loss of function cumulative frequency, and functional
#'         cumulative frequency. Alternatively, if the function was prvoided
#'         with a list of dataframe, return a list of frequency lists, named as
#'         per the input list.
#'
#' @examples
#' vars = read.table(header = TRUE, text = "
#'      AC  AN    CQ
#'      1   1000  missense_variant
#'      1   1000  stop_gained
#'      1   1000  stop_lost
#'      1   1000  synonymous_variant")
#' get_cumulative_frequencies(vars)
#'
#' vars2 = read.table(header = TRUE, text = "
#'      AC  AN    CQ
#'      1   1000  missense_variant
#'      1   1000  stop_gained
#'      1   1000  stop_lost
#'      1   1000  synonymous_variant")
#' var_list = list("first"=vars, "second"=vars2)
#' get_cumulative_frequencies(var_list)
get_cumulative_frequencies <- function(vars) {
    
    # define the VEP consequence types for loss of function and missense variants
    lof_cq = c("stop_gained", "splice_acceptor_variant", "splice_donor_variant",
        "frameshift_variant", "transcript_ablation")
    
    functional_cq = c("stop_lost", "initiator_codon_variant",
        "transcript_amplification", "inframe_insertion", "inframe_deletion",
           "missense_variant", "coding_sequence_variant")
    
    # if we have provided a list of dataframes, then run this function on each
    # of them in turn, and return a list of frequency lists.
    if (!is.data.frame(vars) & is.list(vars)) {
        return(lapply(vars, get_cumulative_frequencies))
    }
    
    vars$frequency = vars$AC/vars$AN
    
    # select the rare variants (possibly this should be done on a site basis,
    # rather than per allele?)
    vars = vars[vars$frequency < 0.01, ]
    
    # find the loss of function variants
    lof_vars = vars[vars$CQ %in% lof_cq, ]
    functional_vars = vars[vars$CQ %in% functional_cq, ]
    
    lof_freq = sum(lof_vars$frequency, na.rm=TRUE)
    functional_freq = sum(functional_vars$frequency, na.rm=TRUE)
    
    # What do we do if the frequency is zero? We won't be able to get meaningful
    # estimates of the enrichment of inherited variants. Estimate the frequency
    # as if the next individual to be included had a heterozygous genotype
    # for the consequence type. Calculate this using the site with the median
    # total alleles. The total allele count is adjusted by two, as if a new
    # biallelic individual had been included in the population.
    if (length(vars$AN[!is.na(vars$AN)]) != 0) {
        if (lof_freq == 0) { lof_freq = 1/(median(vars$AN, na.rm=TRUE) + 2) }
        if (functional_freq == 0) { functional_freq = 1/(stats::median(vars$AN, na.rm=TRUE) + 2) }
    } else {
        if (lof_freq == 0) { lof_freq = NA }
        if (functional_freq == 0) { functional_freq = NA }
    }
    
    frequencies = list(lof=lof_freq, functional=functional_freq)
    
    return(frequencies)
}
