
DATAFREEZE_DIR = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04"

PHENOTYPES_PATH = file.path(DATAFREEZE_DIR, "phenotypes_and_patient_info.txt")

pheno = read.table(PHENOTYPES_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# change the blank syndrome fields to NA
pheno$syndrome[nchar(pheno$syndrome) == 0] = NA


# get the curated list of synrome strings
source("data-raw/syndromes_list.R")

all_matches = list()
for (syndrome in SYNDROMES) {
    syndrome = gsub("[- ]", "\\[ -\\]*", syndrome)
    matches = grepl(syndrome, tolower(pheno$syndrome))
    
    all_matches[[syndrome]] = matches
}

has_match = rowSums(data.frame(all_matches)) > 0

# find the syndromes that I haven't added to the list already
a = data.frame(table(pheno$syndrome[!has_match & ! is.na(pheno$syndrome)]))
head(a[order(a$Freq, decreasing=TRUE), ], 40)

# how many probands have any information in the syndrome field
print(length(pheno$syndrome) - sum(is.na(pheno$syndrome)))

# how many probands have at least one syndrome matched so far
print(sum(has_match))

pheno$syndrome[grepl("wiedemann", tolower(pheno$syndrome))]


lapply(all_matches, sum)


for (syndrome in names(all_matches)) {
    cat("\n")
    cat(syndrome)
    cat("\n")
    print(pheno$syndrome[all_matches[[syndrome]]])
}
