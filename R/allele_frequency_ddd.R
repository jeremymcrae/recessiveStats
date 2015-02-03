
#' find the paths to the VCFs
#'
#' @param parents boolean for whether to restrict to parents only.
#' @param unaffected boolean for whether to restrict to unaffected only.
#' @export
#'
#' @return vector of paths to VCFs
get_vcf_paths <- function(parents=TRUE, unaffected=TRUE) {
    ped_path = file.path(DATAFREEZE_DIR, "family_relationships.txt")
    
    ped = read.table(ped_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    if (parents) { ped = ped[ped$dad_id == 0, ] }
    if (unaffected) { ped = ped[ped$affected == 1, ] }
    
    vcf = ped$path_to_vcf
    
    return(vcf)
}

#' gets variants for a sample in a genome region
#'
#' @param vcf_path path to VCF
#' @param chrom chromosome to get variants from
#' @param start lower bound of chromosomal region to extract variants from
#' @param end upper bound of chromosomal region to extract variants from
#' @export
#'
#' @return dataframe of variants
get_sample_variants <- function(vcf_path, chrom, start, end) {
    
    # extract variants within the region from the VCF
    vars = seqminer::readVCFToListByRange(fileName=vcf_path,
        range=paste(chrom, ":", start, "-", end, sep=""),
        annoType="",
        vcfColumn=c("CHROM", "POS"),
        vcfInfo=c("CQ", "HGNC"),
        vcfIndv="GT")
    
    # tidy up the genotype column
    vars[["GT"]] = as.vector(unlist(vars[["GT"]]))
    
    # convert the list of vectors to a dataframe
    if (length(vars[["POS"]]) > 0) {
        vars = data.frame(vars)
    }  else {
        # make a blank dataframe if we don't have any variants
        vars = data.frame("CHROM"=character(0), "POS"=character(0),
            "CQ"=character(0), "HGNC"=character(0), "GT"=character(0),
            "sampleId"=character(0))
    }
    
    return(vars)
}

#' loads the variants for a given gene from source VCFs
#'
#' @param vcfs vector of vcf paths
#' @param hgnc hgnc symbol as character string
#' @param chrom chromosome as character string
#' @export
#'
#' @return dataframe of variants
get_variants_for_gene <- function(vcfs, hgnc, chrom) {
    
    vcfs = get_vcf_paths()
    
    # find the gene coordinates, sometimes we have multiple genes with the same
    # HGNC symbol on different chromosomes, so we also need to restrict by
    # chromosome
    rows = recessiveStats::gencode[recessiveStats::gencode$gene == hgnc &&
        recessiveStats::gencode$chr == chrom, ]
    stopifnot(nrow(rows) == 1)
    
    all_vars = vector("list", length=length(vcfs))
    pb = txtProgressBar(min=0, max=length(vcfs), style=3)
    for (i in 1:length(vcfs)) {
        path = vcfs[[i]]
        all_vars[[i]] = get_sample_variants(path,
            chrom=rows$chr,
            start=rows$start,
            end=rows$stop)
        
        i = i + 1
        if (i %% 20 == 0) { setTxtProgressBar(pb, i) }
    }
    
    
}

get_functional_vars <- function(vars) {
    # define the VEP consequence types for loss of function and missense variants
    lof_cq = c("stop_gained", "splice_acceptor_variant", "splice_donor_variant",
        "frameshift_variant")
    missense_cq = c("missense_variant", "initiator_codon_variant", "stop_lost",
        "inframe_deletion", "inframe_insertion", "splice_region_variant")
}
