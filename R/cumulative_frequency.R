
#' obtain cumulative allele frequency of rare lof and functional variants
#'
#' @param vars dataframe of variants (one row per allele), which includes the
#'        number of times that allele was observed within the population, as
#'        well as the total number of alleles in the population. Alternatively,
#'        this can be a list of dataframe, each for a different population
#         (e.g. list("EAS"=df(...), "SAS"=df(...))).
#' @param threshold minor allele frequency (MAF) threshold, we exclude variants
#'        with MAF values above or equal to this threshold. This needs to be
#'        matched to the thresh9old used during identification of the
#'        biallelically inherited genotypes.
#' @export
#'
#' @return a list of loss of function cumulative frequency, functional
#'         cumulative frequency and synonymous cumulative frequency.
#'         Alternatively, if the function was prvoided
#'         with a list of dataframe, return a list of frequency lists, named as
#'         per the input list.
#'
#' @examples
#' vars = read.table(header = TRUE, text = "
#'     CHROM  POS  REF  ALT  AC  AN    CQ
#'     1      1    A    G    1   1000  missense_variant
#'     1      2    G    C    1   1000  stop_gained
#'     1      3    T    A    1   1000  stop_lost
#'     1      4    G    T    1   1000  synonymous_variant")
#' get_cumulative_frequencies(vars)
#'
#' vars2 = read.table(header = TRUE, text = "
#'     CHROM  POS  REF  ALT  AC  AN    CQ
#'     1      1    A    G    1   1000  missense_variant
#'     1      2    G    C    1   1000  stop_gained
#'     1      3    T    A    1   1000  stop_lost
#'     1      4    G    T    1   1000  synonymous_variant")
#' var_list = list("first"=vars, "second"=vars2)
#' get_cumulative_frequencies(var_list)
#'
#' threshold = 0.005
#' get_cumulative_frequencies(var_list, threshold)
get_cumulative_frequencies <- function(vars, threshold=0.01) {
    
    # define the VEP consequence types for loss of function and missense variants
    lof_cq = c("stop_gained", "splice_acceptor_variant", "splice_donor_variant",
        "frameshift_variant", "transcript_ablation")
    
    functional_cq = c("stop_lost", "initiator_codon_variant",
        "transcript_amplification", "inframe_insertion", "inframe_deletion",
           "missense_variant", "coding_sequence_variant")
    
    silent_cq = c("synonymous_variant")
    
    # select the rare variants
    vars = remove_high_frequency_vars(vars, threshold)
    
    # if we have provided a list of dataframes, then run this function on each
    # of them in turn, and return a list of frequency lists.
    if (!is.data.frame(vars) & is.list(vars)) {
        return(lapply(vars, get_cumulative_frequencies))
    }
    
    vars$frequency = vars$AC/vars$AN
    
    # find the loss of function variants
    lof_vars = vars[vars$CQ %in% lof_cq, ]
    functional_vars = vars[vars$CQ %in% functional_cq, ]
    silent_vars = vars[vars$CQ %in% silent_cq, ]
    
    lof_freq = cumulative_frequency(lof_vars$frequency)
    functional_freq = cumulative_frequency(functional_vars$frequency)
    silent_freq = cumulative_frequency(silent_vars$frequency)
    
    # What do we do if the frequency is zero? We won't be able to get meaningful
    # estimates of the enrichment of inherited variants. Estimate the frequency
    # as if the next individual to be included had a heterozygous genotype
    # for the consequence type. Calculate this using the site with the median
    # total alleles. The total allele count is adjusted by two, as if a new
    # biallelic individual had been included in the population.
    if (length(vars$AN[!is.na(vars$AN)]) != 0) {
        if (lof_freq == 0) { lof_freq = 1/(stats::median(vars$AN, na.rm=TRUE) + 2) }
        if (functional_freq == 0) { functional_freq = 1/(stats::median(vars$AN, na.rm=TRUE) + 2) }
        if (silent_freq == 0) { silent_freq = 1/(stats::median(vars$AN, na.rm=TRUE) + 2) }
    } else {
        if (lof_freq == 0) { lof_freq = NA }
        if (functional_freq == 0) { functional_freq = NA }
        if (silent_freq == 0) { silent_freq = NA }
    }
    
    frequencies = list(lof=lof_freq, functional=functional_freq,
        synonymous=silent_freq)
    
    return(frequencies)
}

#' remove high frequency variants from the dataset
#'
#' Sometimes a variant is above the frequency threshold in one population, but
#' under it in another. We exclude variants from both populations in these cases.
#'
#' @param vars dataframe of variants (one row per allele), which includes the
#'        number of times that allele was observed within the population, as
#'        well as the total number of alleles in the population. Alternatively,
#'        this can be a list of dataframe, each for a different population
#         (e.g. list("EAS"=df(...), "SAS"=df(...))).
#' @param threshold minor allele frequency (MAF) threshold, we exclude variants
#'        with MAF values above or equal to this threshold. This needs to be
#'        matched to the threshold used during identification of the
#'        biallelically inherited genotypes.
#'
#' @return object with high frequency variants excluded.
#' @export
#'
#' @examples
#' vars = read.table(header = TRUE, text = "
#'     CHROM  POS  REF  ALT  AC  AN    CQ
#'     1      1    A    G    1   1000  missense_variant
#'     1      2    G    C    1   1000  stop_gained
#'     1      3    T    A    1   1000  stop_lost
#'     1      4    G    T    1   1000  synonymous_variant")
#' get_cumulative_frequencies(vars)
#'
#' vars2 = read.table(header = TRUE, text = "
#'     CHROM  POS  REF  ALT  AC  AN    CQ
#'     1      1    A    G    1   1000  missense_variant
#'     1      2    G    C    1   1000  stop_gained
#'     1      3    T    A    1   1000  stop_lost
#'     1      4    G    T    1   1000  synonymous_variant")
#' var_list = list("first"=vars, "second"=vars2)
#' threshold = 0.005
#' remove_high_frequency_vars(var_list, threshold)
remove_high_frequency_vars <- function(vars, threshold) {
    
    # if vars is a list of dataframes, look through them to exclude variants
    # above the threshold in any of the dataframes.
    if (!is.data.frame(vars) & is.list(vars)) {
        # add a frequency column to each population dataset
        freqs = lapply(vars, function(x) x[["AC"]]/x[["AN"]])
        vars = Map(cbind, vars, frequency=freqs)
        
        # get keys for all of the variants that fail the frequency threshold in
        # any population
        keys = lapply(vars, function(x)
            x[x[["frequency"]] > threshold, c("CHROM", "POS", "REF", "ALT")])
        
        # make a dataframe from the keys, so that we can later merge into the
        # original to pick up key matches
        keys = do.call("rbind", keys)
        keys = keys[!duplicated(keys), ]
        keys = na.omit(keys)
        if (nrow(keys) > 0) { keys$remove = TRUE
        } else { keys$remove = logical(0) }
        
        # merge and then remove variants
        vars = lapply(vars, function(x) merge(x, keys, by=c("CHROM", "POS", "REF", "ALT"), all.x=TRUE))
        vars = lapply(vars, function(x) x[is.na(x$remove), ])
        
        # remove unnecessary columns
        vars = lapply(vars, function(x) x[, c("CHROM", "POS", "REF", "ALT", "AC", "AN", "CQ")])
    } else {
        vars$frequency = vars$AC/vars$AN
        vars = vars[vars$frequency < threshold, ]
        
        vars$frequency = NULL
    }
    
    return(vars)
}

#' calculate cumulative allele frequency from a vector of allele frequencies
#'
#' Rather than simply summing the individual allele frequencies to get a
#' cumulative frequency, we need to take the chance of having two variants on
#' one chromosome into account. The chance of observing two variants on the same
#' chromosome is their frequencies multiplied together (e.g. p1 * p2), so the
#' probability of seeing them on independent chromosomes is p1 * (1 - p2).
#' This can be extended beyond two variants, to simply account for the
#' probability of the preceeding variants.
#'
#' @param probs vector of frequencies
#'
#' @return summed cumulative frequency.
#' @export
#'
#' @examples
#'
#' freqs = runif(100, min=1e-7, max=0.01)
#' cumulative_frequency(freqs)
cumulative_frequency <- function(probs) {
    
    if (length(probs) > 0) {
        probs = probs[!is.na(probs)]
    }
    
    total = 0
    for (p in probs) {
        total = total + p * (1 - total)
    }
    return(total)
}
