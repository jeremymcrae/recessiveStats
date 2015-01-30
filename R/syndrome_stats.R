
DATAFREEZE_DIR = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04"
PHENOTYPES_PATH = file.path(DATAFREEZE_DIR, "phenotypes_and_patient_info.txt")

# get the curated list of syndrome strings
source("data-raw/syndromes_list.R")

get_pheno_data <- function(path) {
    pheno = read.table(path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    # change the blank syndrome fields to NA
    pheno$syndrome[nchar(pheno$syndrome) == 0] = NA
    
    return(pheno)
}

# search for all the strings in the syndrome column
get_syndrome_matches <- function(syndromes, pheno) {
    all_matches = list()
    for (syndrome in syndromes) {
        syndrome = gsub("[- ]", "\\[ -\\]*", syndrome)
        matches = grepl(syndrome, tolower(pheno$syndrome))
        
        all_matches[[syndrome]] = matches
    }
    
    return(all_matches)
}

# find which syndrome strings haven't been captured in the SYNDROMES list.
# These should mostly contain strings that are only seen once in These
# probands, recurrent syndromes should be shifted into SYNDROMES in syndromes_list.R
show_missing_syndromes <- function(has_match, pheno) {
    # find the syndromes that I haven't added to the list already
    missing = data.frame(table(pheno$syndrome[!has_match & ! is.na(pheno$syndrome)]))
    
    # show the most recurrent syndrome strings (these should all have 1-2
    # entries)
    print(head(missing[order(a$Freq, decreasing=TRUE), ], 40))
    
    # how many probands have any information in the syndrome field
    print(length(pheno$syndrome) - sum(is.na(pheno$syndrome)))
    
    # how many probands have at least one syndrome matched so far
    print(sum(has_match))
    
    # # show the exact matches for each syndrome string. Note: this misses the
    # # strings that did NOT match the regex.
    # for (syndrome in names(all_matches)) {
    #     cat("\n")
    #     cat(syndrome)
    #     cat("\n")
    #     print(pheno$syndrome[all_matches[[syndrome]]])
    # }
    
    return(missing)
}

pheno = get_pheno_data(PHENOTYPES_PATH)
all_matches = get_syndrome_matches(SYNDROMES, pheno)
has_match = rowSums(data.frame(all_matches)) > 0

missing = show_missing_syndromes(has_match, pheno)
