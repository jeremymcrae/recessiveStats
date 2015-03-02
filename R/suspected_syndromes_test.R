
#' open suspected syndromes
#'
#' @param path path to curated list of suspected syndromes per proband
#' @export
#'
#' @return list of probands, where each probandf has a vector of suspected
#'     syndromes.
open_suspected_syndromes <- function(path) {
    syndromes = read.table(path, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
    
    syndromes$syndrome = strsplit(syndromes$syndrome, ";")
    
    # convert the syndromes into a list, where the names of the list entries are
    # the proband IDs, so that we can extract syndromes per proband by proband
    # ID.
    syndromes = split(syndromes, syndromes$patient_id)
    syndromes = lapply(syndromes, function(x) as.vector(x[["syndrome"]])[[1]])
    
    return(syndromes)
}

#' calculate the similarity of suspected syndromes between two probands
#'
#' @param syndrome_rates dataframe of syndrome terms and the rates at which they
#'     are used within the DDD probands.
#' @param proband_1 vector of syndrome terms for a proband
#' @param proband_2 vector of syndrome terms for a proband
#' @export
#'
#' @return the similarity of the suspected syndromes for two syndromes, as rate
#'     between 0 and 1.
calculate_proband_similarity <- function(syndrome_rates, proband_1, proband_2) {
    
    proband_1 = unlist(proband_1)
    proband_2 = unlist(proband_2)
    
    in_both = as.vector(proband_1[proband_1 %in% proband_2])
    
    # if none of the syndromes match between the two probands, return the
    # highest possible rate, which is 1
    if (length(in_both) == 0) { return(1) }
    
    rates = syndrome_rates$rate[syndrome_rates$syndrome %in% in_both]
    rate = min(rates)
    
    return(rate)
}

#' calculate the similarity of suspected syndrome terms within a group of probands
#'
#' @param probands list of probands, with their suspected syndromes as vectors
#' @param syndrome_rates dataframe of syndrome terms and the rates at which they
#'     are used within the DDD probands.
#' @export
#'
#' @return summed similarity across the probands
calculate_gene_similarity <- function(probands, syndrome_rates) {
    
    similarity = c()
    for (x in 1:length(probands)) {
        for (y in 1:length(probands)) {
            if (x == y) { next }
            
            temp = calculate_proband_similarity(syndrome_rates, probands[x], probands[y])
            similarity = c(similarity, temp)
        }
    }
    
    similarity = -sum(log10(similarity))
    
    return(similarity)
}

#' get the chance of a group of probands sharing their syndrome terms
#'
#' @param syndromes list of syndrome terms per proband
#' @param probands vector of proband IDs for a gene
#' @param n_sims number of iterations to run
#' @export
#'
#' @return P-value for the chance that the probands share their syndrome terms
test_gene <- function(syndromes, probands, n_sims=10000) {
    
    # calculate the rate at which each syndrome term is used within the population
    term = unlist(syndromes)
    syndrome_rates = data.frame(table(term))
    syndrome_rates$syndrome = syndrome_rates$term
    syndrome_rates$rate = syndrome_rates$Freq/sum(syndrome_rates$Freq)
    
    # calculate the similarity between the observed probands
    proband_terms = syndromes[probands]
    observed_similarity = calculate_gene_similarity(proband_terms, syndrome_rates)
    
    distribution = rep(NA, n_sims)
    for (x in 1:n_sims) {
        sampled = sample(syndromes, length(probands))
        distribution[x] = calculate_gene_similarity(sampled, syndrome_rates)
    }
    
    distribution = sort(distribution)
    p_value = (sum(findInterval(x, observed_similarity)))/(n_sims + 1)
    
    if (p_value == 0) { p_value = 1/(n_sims + 1)}
    
    return(p_value)
}
