# tally probands with recessive variants (as compound hets with
# loss-of-function, or functional events). The inital list of variants comes
# from clinical filtering output, run using all genes rather than restricting it
# to the known developmental disorder genes.

library(reshape)
library(jsonlite)

DATAFREEZE_DIR = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04"
TRIOS_PATH = file.path(DATAFREEZE_DIR, "family_relationships.txt")
PHENOTYPES_PATH = file.path(DATAFREEZE_DIR, "phenotypes_and_patient_info.txt")
SANGER_IDS_PATH = file.path(DATAFREEZE_DIR, "person_sanger_decipher.txt")
CLINICALFILTER_OUTPUT = "/nfs/users/nfs_j/jm33/clinical_reporting_output/clinical_reporting.2015-05-11.synonymous_variants.txt"
LIKELY_DIAGNOSED = "/lustre/scratch113/projects/ddd/users/jm33/ddd_likely_diagnosed.txt"

DATA_DIR = "/nfs/users/nfs_j/jm33/apps/recessiveStats/data-raw"
OUTPUT_COUNTS_PATH = file.path(DATA_DIR, "recessive_counts_by_gene.silent.txt")
OUTPUT_PROBANDS_JSON = file.path(DATA_DIR, "recessive_probands_by_gene.silent.json")

open_variants <- function() {
    diagnosed = read.table(LIKELY_DIAGNOSED, header=TRUE, stringsAsFactors=FALSE)
    variants = read.table(CLINICALFILTER_OUTPUT, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    # drop the likely diagnosed probands
    variants = variants[!(variants$proband %in% diagnosed$person_id), ]
    
    return(variants)
}

#' opens a file defining the DDD individuals, and their families
open_ddd_families <- function(path) {
    families = read.table(path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    families$affected = NULL
    families$path_to_vcf = NULL
    families$sex = NULL
    
    return(families)
}

open_ages <- function() {
    phenotypes = read.table(PHENOTYPES_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    ages = phenotypes[, c("patient_id", "decimal_age")]
    
    # make sure we have DDD IDs for each proband
    sanger_ids = read.table(SANGER_IDS_PATH, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    sanger_ids$sanger_id = NULL
    sanger_ids = sanger_ids[!duplicated(sanger_ids), ]
    
    ages = merge(ages, sanger_ids, by.x="patient_id", by.y="decipher_id", all.x=TRUE)
    
    return(ages)
}

#' remove variants that are shared between multiple probands of a family
get_independent_de_novos <- function(variants) {
    families = open_ddd_families(TRIOS_PATH)
    ages = open_ages()
    families = merge(families, ages, by.x="individual_id", by.y="person_stable_id", all.x=TRUE)
    
    # merge the family IDs with the variants
    variants = merge(variants, families, by.x="proband", by.y="individual_id", all.x=TRUE)
    
    # sort the variants by family ID, and age (with eldest siblings sorted first)
    # so that when we remove duplicates in families, we retain the eldest sibling
    variants = variants[order(variants[["family_id"]], variants[["decimal_age"]], decreasing=TRUE), ]
    
    # get a dataframe of genes per proband, where we have only one row per
    # proband per gene, so that later we can find families with more than one
    # proband without having to worry about each gene having multiple variants
    # from the compound hets.
    probands = variants[, c("proband", "family_id", "chrom", "gene")]
    probands = probands[!duplicated(probands),]
    
    family_genes = probands[, c("family_id", "chrom", "gene")]
    probands$duplicates = duplicated(family_genes)
    
    # restrict ourselves to the non-duplicates
    variants = merge(variants, probands, by=c("proband", "family_id", "chrom", "gene"), all.x=TRUE)
    variants = variants[!(variants$duplicates), ]
    
    return(variants)
}

count_for_uncertain_probands <- function(variants, uncertain) {
    # some compound het probands have two functional varinats in a gene, and one
    # lof variant. We can't distingush between between the proband being
    # lof/func, or func/func just from the numbers. We need to check the trio
    # genotypes for each variant, to identify which categories have different
    # trio genotypes.
    
    status = rep(NA, nrow(uncertain))
    for (row_num in 1:nrow(uncertain)) {
        row = uncertain[row_num, ]
        
        gene = variants[variants$proband == row$proband & variants$gene == row$gene, ]
        func_genotypes = gene$trio_genotype[gene$functional == "func"]
        lof_genotypes = gene$trio_genotype[gene$functional == "lof"]
        silent_genotypes = gene$trio_genotype[gene$functional == "silent"]
        
        # If we have two "functional" variants with differing genotypes, then
        # this could be a func/func category.
        if (length(table(silent_genotypes)) > 1) {
            status[row_num] = "silent_silent"
        }
        
        # If we have two "functional" variants with differing genotypes, then
        # this could be a func/func category.
        if (length(table(func_genotypes)) > 1) {
            status[row_num] = "func_func"
        }
        
        for (genotype in lof_genotypes) {
            if (length(func_genotypes[func_genotypes != genotype]) > 0) {
                status[row_num] = "lof_func"
            }
        }
        
        # If we have two "lof" variants with differing genotypes, this is a
        # lof/lof category.
        if (length(table(lof_genotypes)) > 1) {
            status[row_num] = "lof_lof"
        }
    }
    
    return(status)
}

# find the LoFs recessively inherited as variants on their own
count_homozygous_lofs_per_gene <- function(variants) {
    single_lofs = variants[grepl("single_variant", variants$result) &
        grepl("Biallelic", variants$inheritance) &
        variants$functional == "lof", ]
        
    # get the homozygous recessive single LoFs,
    single_lofs = single_lofs[sapply(strsplit(single_lofs$trio_genotype, "/"), "[", 1) == "2", ]
    
    # and only count one per gene (ie if a proband has two homozygous recessive
    # LoFs in a gene, we count that as a single proband)
    single_lofs = single_lofs[!duplicated(single_lofs[, c("proband", "chrom", "gene")]), ]
    
    single_lofs = single_lofs[, c("proband", "gene", "chrom")]
    single_lofs$category = "homozygous_lof"
    
    return(single_lofs)
}

# find the LoFs recessively inherited as variants on their own
count_homozygous_silent_per_gene <- function(variants) {
    single_lofs = variants[grepl("single_variant", variants$result) &
        grepl("Biallelic", variants$inheritance) &
        variants$functional == "silent", ]
        
    # get the homozygous recessive single LoFs,
    single_lofs = single_lofs[sapply(strsplit(single_lofs$trio_genotype, "/"), "[", 1) == "2", ]
    
    # and only count one per gene (ie if a proband has two homozygous recessive
    # LoFs in a gene, we count that as a single proband)
    single_lofs = single_lofs[!duplicated(single_lofs[, c("proband", "chrom", "gene")]), ]
    
    single_lofs = single_lofs[, c("proband", "gene", "chrom")]
    single_lofs$category = "homozygous_silent"
    
    return(single_lofs)
}

# find compound hets in genes
count_compound_hets <- function(variants) {
    # get only the compound het variants
    compound_hets = variants[grepl("compound_het", variants$result), ]
    
    # figure out whether each proband, for each gene is a lof/lof, a lof/func, or a
    # func/func
    probands = cast(compound_hets, proband + gene + chrom ~ functional, value="result", length)
    probands$category = NA
    probands$category[probands$func == 0 & probands$silent > 1] = "silent_silent"
    probands$category[probands$func == 0 & probands$lof > 1] = "lof_lof"
    probands$category[probands$func > 1 & probands$lof == 0] = "func_func"
    probands$category[probands$func == 1 & probands$lof == 1] = "lof_func"
    
    # account for the probands who have 2 or more "lof" variants, and 1 or more
    # "func" variants (or vice-versa).
    uncertain = probands[is.na(probands$category), ]
    probands$category[is.na(probands$category)] = count_for_uncertain_probands(compound_hets, uncertain)
    
    probands = probands[, c("proband", "gene", "chrom", "category")]
    
    return(probands)
}

count_probands_per_gene <- function(variants) {
    
    homozygous_synonymous = count_homozygous_silent_per_gene(variants)
    compound_hets = count_compound_hets(variants)
    
    probands = rbind(homozygous_synonymous, compound_hets)
    counts = cast(probands, gene + chrom ~ category, value="proband", length)
    
    counts[is.na(counts)] = 0
    counts$silent_sum = rowSums(counts[, c("silent_silent", "homozygous_silent")])
    
    # only include recurrent genes
    counts = counts[counts$silent_sum > 1, ]
    
    # drop the func/func probands, since there are far too many of them, and we
    # currently don't use them in analyses
    probands = probands[probands$category != "func_func", ]
    
    values = list(counts=counts, probands=probands)
    
    return(values)
}

set_functional_state <- function(variants) {
    lof_cq = c("stop_gained", "splice_acceptor_variant", "splice_donor_variant",
        "frameshift_variant")
    missense_cq = c("missense_variant", "initiator_codon_variant", "stop_lost",
        "inframe_deletion", "inframe_insertion", "splice_region_variant",
        "coding_sequence_variant")
    synonymous_cq = c("synonymous_variant")
    
    # define the two functional categories, loss-of-function versus functional
    variants$functional = NA
    variants$functional[grepl(paste(missense_cq, collapse="|"), variants$consequence)] = "func"
    variants$functional[grepl(paste(lof_cq, collapse="|"), variants$consequence)] = "lof"
    variants$functional[grepl(paste(synonymous_cq, collapse="|"), variants$consequence)] = "silent"
    
    return(variants)
}

main <- function() {
    variants = open_variants()
    variants = get_independent_de_novos(variants)
    variants = set_functional_state(variants)
    values = count_probands_per_gene(variants)
    counts = values$counts
    probands = values$probands
    
    proband_genes = sort(unique(probands$gene))
    probands_json = sapply(proband_genes, function(x) probands$proband[probands$gene == x])
    probands_json = toJSON(probands_json, pretty=TRUE)
    
    write.table(counts, file=OUTPUT_COUNTS_PATH, row.names=FALSE, quote=FALSE, sep="\t")
    
    fileConn = file(OUTPUT_PROBANDS_JSON)
    writeLines(probands_json, fileConn)
    close(fileConn)
}

main()
