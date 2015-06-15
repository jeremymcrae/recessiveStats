
library(Cairo)
library(dplyr)

AUTOZYGOSITY_DIR = "/nfs/users/nfs_j/jm33/apps/recessiveStats/data-raw/autozygosity"
KINSHIP_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/kinship_and_pca_trios.txt"
FAMILIES_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/family_relationships.txt"
SIGNIFICANT_PATH = "/nfs/users/nfs_j/jm33/apps/recessiveStats/results/recessive.combined_tests.ver2.txt"

#' identify the probands from consanguinous parents.
get_consanguinous_probands <- function(path) {
    cohort = read.table(path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    probands = cohort[cohort$king_kinship > 0, "proband_stable_id"]
    
    return(probands)
}

#' read the autozygosity output (which are individual files for each proband)
get_autozygous_regions <- function(path) {
    paths = Sys.glob(file.path(path, "*"))
    
    regions = vector("list", length(paths))
    
    for (pos in 1:length(paths)) {
        temp = read.table(paths[pos], sep="\t", header=TRUE,
            colClasses=c("character", "character", "numeric", "numeric"))
        regions[[pos]] = temp
    }
    
    regions = data.frame(rbind_all(regions))
    
    return(regions)
}

#' load a dataframe of family IDs andf individual IDs
get_families <- function(path) {
    families = read.table(path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    families = families[, c("family_id", "individual_id")]
    
    return(families)
}

#' for a list of probands, make sure we only count one proband from each family
get_independent_probands <- function(families, probands) {
    
    # find the family IDs for all the probands, then get the subset of probands
    # where the family ID is not a duplicate.
    family_ids = families$family_id[families$individual_id %in% probands]
    independent_probands = probands[!duplicated(family_ids)]
    
    return(independent_probands)
}

get_autozygosity_per_gene <- function(path, regions, families, subset=NULL) {
    genes = recessiveStats::gencode[recessiveStats::gencode$gene_type == "protein_coding", ]
    genes = genes[order(genes$gene), ]
    
    # Define the cohort size (we can't rely upon the number of unique probands
    # in the regions dataframe, since many probands will not have autozygous
    # regions).
    # probands_n = length(Sys.glob(file.path(path, "*")))
    probands_n = 3071
    
    # if we want to estimate the rates for consanguinous probands only, we take
    # a subset of the autozygous regions
    if (!is.null(subset)) {
        probands_n = length(subset)
        regions = regions[regions$sample_id %in% subset, ]
    }
    
    rates = data.frame(hgnc=genes$gene, chrom=genes$chr, count=NA, rate=NA)
    for (pos in 1:nrow(genes)) {
        # get the coordinates for the gene
        chrom = genes$chr[pos]
        start_pos = genes$start[pos]
        end_pos = genes$stop[pos]
        
        # find the autozygous regions that overlap the gene region
        in_chrom = regions[regions$chrom == chrom, ]
        overlapping = in_chrom[in_chrom$end_pos > start_pos
            & in_chrom$start_pos < end_pos, ]
        
        probands = unique(overlapping$sample_id)
        probands = get_independent_probands(families, probands)
        rates$count[pos] = length(probands)
        rates$rate[pos] = rates$count[pos]/probands_n
    }
    
    # remove the chrX rates, since all of the male probands appear autozygous
    rates = rates[rates$chrom != "X", ]
    
    return(rates)
}

plot_autozygosity <- function(all_probands, signif_path) {
    # exclude chrX, since males always appear autozygous, and I hadn't removed
    # them for this analysis run
    all_probands = all_probands[!all_probands$chrom == "X", ]
    
    # plot the autozygosity rates
    Cairo(file="autozygosity_rates.pdf", type="pdf", height=15, width=15, units="cm")
    hist(log10(all_probands$rate),
        xlab="Autozygosity rate in DDD probands",
        main="autozygosity rates for autosomal genes",
        xlim=c(-3.5, -1), las=1, cex=1.3, cex.lab=1.3, cex.axis=1.3, tck=-0.03)
    
    # also plot the rates for the most significantly associated genes, to show
    # that they fall within the typical spectrum of autozygosity
    signif_genes = read.table(signif_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    genes = head(signif_genes$hgnc, 8)
    for (gene in genes) {
        rate = all_probands$rate[all_probands$hgnc == gene]
        abline(v=log10(rate), col="gray")
    }
    
    dev.off()
}

families = get_families(FAMILIES_PATH)
regions_per_proband = get_autozygous_regions(AUTOZYGOSITY_DIR)
all_autozygous = get_autozygosity_per_gene(AUTOZYGOSITY_DIR, regions_per_proband, families)

write.table(all_autozygous, file="autozygosity.all_genes.all_probands.tsv", sep="\t", row.names=FALSE, quote=FALSE)

plot_autozygosity(all_autozygous, SIGNIFICANT_PATH)
