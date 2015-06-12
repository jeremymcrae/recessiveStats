
library(Cairo)
library(dplyr)

AUTOZYGOSITY_DIR = "/nfs/users/nfs_j/jm33/apps/recessiveStats/data-raw/autozygosity"
KINSHIP_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/kinship_and_pca_trios.txt"
FAMILIES_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/family_relationships.txt"

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
    
    if (any(duplicated(family_ids))) { print("some sibs found") }
    
    independent_probands = probands[!duplicated(family_ids)]
    
    return(independent_probands)
}

get_autozygosity_per_gene <- function(path, regions, families, subset=NULL) {
    genes = recessiveStats::gencode[recessiveStats::gencode$gene_type == "protein_coding", ]
    genes = genes[order(genes$gene), ]
    
    # Define the cohort size (we can't rely upon the number of unique probands
    # in the regions dataframe, since many probands will not have autozygous
    # regions).
    probands_n = length(Sys.glob(file.path(path, "*")))
    
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

plot_autozygosity <- function(all_probands, consang) {
    # exclude chrX, since males always appear autozygous, and I hadn't removed
    # them for this analysis run
    all_probands = all_probands[!all_probands$chrom == "X", ]
    all_density = density(all_probands$rate, na.rm=TRUE)
    consang = consang[!consang$chrom == "X", ]
    consang_density = density(consang$rate, na.rm=TRUE)

    max_y = max(all_density$y, consang_density$y)
    max_x = max(all_density$x, consang_density$x)
    min_x = min(all_density$x, consang_density$x)

    # plot the autozygosity rates
    Cairo(file="autozygosity_rates.pdf", type="pdf", height=15, width=15, units="cm")
    plot(all_density,
        xlab="Autozygosity rate in DDD probands",
        main="autozygosity rates for autosomal genes",
        xlim=c(min_x, max_x), ylim=c(0, max_y), las=1,
        cex=1.3, cex.lab=1.3, cex.axis=1.3, tck=-0.03)
    lines(consang_density, col="red")
    legend("topright", legend=c("all probands", "consanguinous probands"),
        col=c("black", "red"), bty="n", lwd=1)
    dev.off()
}

families = get_families(FAMILIES_PATH)
consanguinous = get_consanguinous_probands(KINSHIP_PATH)
regions_per_proband = get_autozygous_regions(AUTOZYGOSITY_DIR)
all_autozygous = get_autozygosity_per_gene(AUTOZYGOSITY_DIR, regions_per_proband, families)
consang_autozygous = get_autozygosity_per_gene(AUTOZYGOSITY_DIR, regions_per_proband, families, subset=consanguinous)

write.table(all_autozygous, file="autozygosity.all_genes.all_probands.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(consang_autozygous, file="autozygosity.all_genes.consanguinous_probands.tsv", sep="\t", row.names=FALSE, quote=FALSE)

plot_autozygosity(all_autozygous, consang_autozygous)
