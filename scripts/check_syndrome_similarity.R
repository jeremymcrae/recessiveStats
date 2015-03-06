
library(recessiveStats)
library(rjson)

decipher_ids_path = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/person_sanger_decipher.txt"
# probands_path = "/nfs/users/nfs_j/jm33/apps/hpo_similarity/data/probands_by_gene.json"
probands_path = "/nfs/users/nfs_j/jm33/apps/hpo_similarity/data/recessive_probands_by_gene.json"
suspected_syndromes_path = "/nfs/users/nfs_j/jm33/apps/recessiveStats/data-raw/syndromes_final.tsv"

decipher_ids = read.table(decipher_ids_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
genes = fromJSON(file=probands_path)

# swap all the DDD IDs to decipher IDs
genes = lapply(genes, function(x) unique(decipher_ids$decipher_id[decipher_ids$person_stable_id %in% x]))
syndromes = open_suspected_syndromes(suspected_syndromes_path)

for (gene in names(genes)) {
    probands = genes[[gene]]
    p_value = test_gene(syndromes, probands)
    
    if (!is.na(p_value)) {
        cat(paste(gene, "\t", p_value, "\n", sep=""))
    }
}
