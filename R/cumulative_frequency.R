
#' obtain cumulative allele frequency of rare lof and functional variants
#'
#' @param vars dataframe of variants (one row per allele), which includes the
#'     number of times that allele was observed within the population, as well as
#'     the total number of alleles in the population.
#' @export
#'
#' @return a list of loss of function cumulative frequency, and functional
#'     cumulative frequency.
get_cumulative_frequencies <- function(vars) {
    
    # define the VEP consequence types for loss of function and missense variants
    lof_cq = c("stop_gained", "splice_acceptor_variant", "splice_donor_variant",
        "frameshift_variant")
    
    vars$frequency = vars$AC/vars$AN
    
    # select the rare variants (possibly this should be done on a site basis,
    # rather than per allele?)
    vars = vars[vars$frequency < 0.01, ]
    
    # find the loss of function variants
    lof_vars = vars[vars$CQ %in% lof_cq, ]
    
    lof_freq = sum(lof_vars$frequency, na.rm=TRUE)
    functional_freq = sum(vars$frequency, na.rm=TRUE)
    
    # What do we do if the frequency is zero? We won't be able to get meaningful
    # estimates of the enrichment of inherited variants. Estimate the frequency
    # as if the next individual to be included had a heterozygous genotype
    # for the consequence type. Calculate this using the site with the fewest
    # total alleles. The total allele count is adjusted by two, as if a new
    # biallelic individual had been included in the population.
    if (lof_freq == 0) { lof_freq = 1/min(vars$AN) + 2 }
    if (functional_freq == 0) { functional_freq = 1/min(vars$AN) + 2 }
    
    frequencies = list(lof=lof_freq, functional=functional_freq)
    
    return(frequencies)
}
