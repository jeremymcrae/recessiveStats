# script to tidy up the suspected syndromes for use in finding similarity
# between suspected syndromes in groups of probands.

library(gdata)

SYNDROMES_PATH = "data-raw/syndromes_cleaned.xlsx"
OUTPUT_PATH = "data-raw/syndromes_final.tsv"

#' open the phenotype data,
#'
#' @param path path to the proband's phenotype file.
#' @export
#'
#' @return dataframe of phenotype data.
get_syndrome_data <- function(path) {
    syndromes = read.xls(path, stringsAsFactors=FALSE)
    
    # remove all the inappropriate terms (defined in the excel table) from the
    # required probands
    syndromes$matched_terms[syndromes$remove != ""] = apply(syndromes[syndromes$remove != "", ], 1, function(x) gsub(paste("[;]+", x[["remove"]], "[;]+", sep=""), ";", x[["matched_terms"]]))
    
    syndromes$matched_terms[syndromes$remove != ""] = apply(syndromes[syndromes$remove != "", ], 1, function(x) gsub(paste("[;]*", x[["remove"]], "[;]*", sep=""), "", x[["matched_terms"]]))
    
    # there has to be an easier way to do this:
    syndromes$matched_terms[syndromes$matched_terms != "" & syndromes$add != ""] = paste(syndromes$matched_terms[syndromes$matched_terms != "" & syndromes$add != ""], syndromes$add[syndromes$matched_terms != "" & syndromes$add != ""], sep=";")
    
    syndromes$matched_terms[syndromes$matched_terms == "" & syndromes$add != ""] = syndromes$add[syndromes$matched_terms == "" & syndromes$add != ""]
    
    syndromes$syndrome = syndromes$matched_terms
    
    # change the blank syndrome fields to NA
    syndromes$syndrome[nchar(syndromes$syndrome) == 0] = NA
    
    # remove all the rows with missing syndrome values
    syndromes = syndromes[!(is.na(syndromes$syndrome)), ]
    
    syndromes = syndromes[, c("patient_id", "syndrome")]
    
    return(syndromes)
}


main <- function() {
    syndromes = get_syndrome_data(SYNDROMES_PATH)
    
    write.table(syndromes, file=OUTPUT_PATH, sep="\t", quote=FALSE, row.names=FALSE)
}

main()
