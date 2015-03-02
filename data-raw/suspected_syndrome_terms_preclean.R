# script to tidy up the suspected syndromes for use in finding similarity
# between suspected syndromes in groups of probands.

#' open the phenotype data,
#'
#' @param path path to the proband's phenotype file.
#' @export
#'
#' @return dataframe of phenotype data.
get_pheno_data <- function(path) {
    pheno = read.table(path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    # change the blank syndrome fields to NA
    pheno$syndrome[nchar(pheno$syndrome) == 0] = NA
    
    # remove all the rows with missing syndrome values
    pheno = pheno[!(is.na(pheno$syndrome)), ]
    
    return(pheno)
}

#' search for all the strings in the syndrome column
#'
#' @param syndromes vector of regex strings.
#' @param pheno dataframe of phenotypes for probands, containing a "syndrome"
#'     column.
#' @export
#'
#' @return dataframe, where each column is for a different regex string, and the
#'    rows correspond to the rows of th pheno dataframe
get_syndrome_matches <- function(syndromes, pheno) {
    all_matches = list()
    
    # convert spaces and hyphens to possible space and hyphen matches, since
    # some clinicians use spaces, others use hyphens, while others omit both
    syndromes = gsub("[- ]", "\\[ -\\]*", syndromes)
    
    # find all the matches to the syndrome regex strings
    all_matches = lapply(syndromes, grepl, x=tolower(pheno$syndrome))
    
    # convert the list to a dataframe
    all_matches = data.frame(all_matches)
    names(all_matches) = syndromes
    
    return(all_matches)
}

#' find which syndrome strings haven't been captured in the SYNDROMES list.
#'
#' These should mostly contain strings that are only seen once in These
#' probands, recurrent syndromes should be shifted into SYNDROMES in
#' syndromes_list.R
#'
#' @param all_matches dataframe of booleans for different syndrome regex strings
#' @param pheno dataframe of phenotypes for probands, containing a "syndrome"
#'     column.
#' @export
#'
#' @return dataframe
show_missing_syndromes <- function(all_matches, pheno) {
    
    has_match = rowSums(all_matches) > 0
    
    # find the syndromes that I haven't added to the list already
    missing = data.frame(table(pheno$syndrome[!has_match & ! is.na(pheno$syndrome)]))
    
    # show the most recurrent syndrome strings (these should all have 1-2
    # entries)
    print(head(missing[order(missing$Freq, decreasing=TRUE), ], 40))
    
    # how many probands have any information in the syndrome field, and how many
    # probands have at least one syndrome matched so far
    n_probands = length(pheno$syndrome) - sum(is.na(pheno$syndrome))
    n_with_match = sum(has_match)
    cat(paste("n probands:", n_probands, "\nwith match:", n_with_match, "\n", sep=" "))
    
    return(missing)
}

main <- function() {
    pheno = get_pheno_data(PHENOTYPES_PATH)
    all_matches = get_syndrome_matches(recessiveStats::SYNDROMES$regex, pheno)

    # rename the over-long syndrome text for one proband to something more readable
    pheno$syndrome[grepl("SOX3 causing X-linked isolated growth hormone deficiency with MR  IGHDIII caused isolated growth hormone deficiency with or without agammaglobuliaemia", pheno$syndrome)] = "X-linked with MR  isolated growth hormone deficiency/ Weil-Marchesani syndrome XLMR syndrome with seizures  hypogammaglobulinaemia"

    pheno$matched_terms = apply(all_matches, 1, function(x) recessiveStats:: SYNDROMES$name[as.vector(unlist(x))])
    pheno$matched_terms = sapply(pheno$matched_terms, paste, collapse=";")
    write.table(pheno[, c("patient_id", "syndrome", "matched_terms")], file="test.txt", sep="\t", row.names=FALSE, quote=FALSE)

    missing = show_missing_syndromes(all_matches, pheno)
}
