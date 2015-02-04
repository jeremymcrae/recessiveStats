
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

#' find where the required data lies within the file
#'
#' Seek in a file to the section containing the data that we want. We start
#' with a reasonable guess about where in the file the data is located, but
#' we seek up or downstream of that (depending on whether the data at that
#' position is before or after the data that we want), then bisect to pin
#' down where to start loading data from.
#'
#' @param path path to the datafile
#' @param start_pos start position of the variant
#' @export
#' @return the approximate byte offset in the file (to within a few hundred
#'    bytes) where the required data starts
find_file_pos <- function(path, start_pos) {
    
    # seek to the correct part of the file
    file_con = file(path, open="r")
    
    # make sure we don't use negative seek positions, or go beyond the file size
    lower = 0
    upper = file.info(path)$size
    
    while (upper - lower > 100) {
        mid_point = lower + (upper - lower)/2
        prev_pos = seek(file_con, mid_point, origin="start")
        
        read = scan(file_con, what=character(), skip=1, sep="\t", nmax=2, quiet=TRUE)
        temp_start = as.numeric(read[2])
        
        # narrow the boundaries to half the previous size
        if (temp_start > start_pos) {upper = mid_point}
        if (temp_start <= start_pos) {lower = mid_point}
    }
    
    close(file_con)
    
    # make sure the seek position is within the file
    seek_position = max(lower, 0)
    seek_position = min(seek_position, file.info(path)$size)
    
    return(seek_position)
}

#' loads DDD VEP annotations lines for variants
#'
#' @param chrom chromosome as character string
#' @param start start of gene range
#' @param end end of gene range
#' @export
#'
#' @return dataframe of VEP annotations
get_lines_from_vep_vcf <- function(chrom, start, end) {
    vep_dir = "/lustre/scratch113/projects/ddd/users/ddd/test_msVCF_annotate_vep/vep_annotated_vcfs/results_vcfs"
    vep_path = file.path(vep_dir, paste(chrom, ".txt", sep=""))
    
    # get the byte offsets for where the data lies within the file (this is a
    # workaround to speed access to the data, while not having tabix indexes).
    file_start = max(find_file_pos(vep_path, start) - 10000, 0)
    file_end = min(find_file_pos(vep_path, end) + 10000, file.info(vep_path)$size)
    
    # load the lines from the file (assume the lines are about 400 bytes long)
    line_length = 400
    con = file(vep_path, open="r")
    seek(con, file_start, origin="start")
    lines_n = (file_end - file_start) / line_length
    vep = read.table(con, skip=1, sep="\t", nrows=lines_n, stringsAsFactors=FALSE)
    names(vep) = c("chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format")
    close(con)
    
    # make sure we have captured variants that lie outside the gene sites,
    # otherwise we don't know if we have missed any variants inside the gene range
    stopifnot(vep$pos[1] < start)
    stopifnot(vep$pos[nrow(vep)] > end)
    
    return(vep)
}

#' gets the DDD VEP annotations for variants
#'
#' @param chrom chromosome as character string
#' @param start start of gene range
#' @param end end of gene range
#' @export
#'
#' @return dataframe of VEP annotations
get_ddd_vep_annotations <- function(chrom, start, end) {
    
    vep = get_lines_from_vep_vcf(chrom, start, end)
    vep = vep[vep$pos > start & vep$pos < end, ]
    
    vep$HGNC = NA
    vep$CQ = NA
    for (pos in 1:nrow(vep)) {
        row = vep[pos, ]
        info = unlist(strsplit(row$info, ";"))
        # some rows lack a "HGNC" field, but do have a "HGNC_ALL" field, so we
        # shall use this instead
        if (!any(grepl("HGNC=", info))) { info = gsub("HGNC_ALL", "HGNC", info) }
        
        # make sure there is a HGNC field (I'm not sure if this will be used,
        # since we are looking at variants within gene regions)
        if (!any(grepl("HGNC=", info))) { info = c(info, "HGNC=") }
        
        # get the HGNC and consequence fields from the INFO
        vep$HGNC[pos] = unlist(strsplit(info[grepl("HGNC=", info)], "="))[2]
        vep$CQ[pos] = unlist(strsplit(info[grepl("CQ=", info)], "="))[2]
    }
    
    # trim the unnecesary columns, it'll make a cleaner merge later
    vep$info = NULL
    vep$qual = NULL
    vep$filter = NULL
    vep$format = NULL
    vep$id = NULL
    
    return(vep)
}

#' loads the variants for a given gene from source VCFs
#'
#' @param hgnc hgnc symbol as character string
#' @param chrom chromosome as character string
#' @export
#'
#' @return dataframe of variants
get_ddd_variants_for_gene <- function(hgnc, chrom) {
    
    # vcfs = get_vcf_paths()
    
    vcfs_dir = "/lustre/scratch114/projects/ddd/release/20140912/final"
    vcf_path = Sys.glob(file.path(vcfs_dir, paste(chrom, "\\:1-*.vcf.gz", sep="")))
    
    # find the gene coordinates, sometimes we have multiple genes with the same
    # HGNC symbol on different chromosomes, so we also need to restrict by
    # chromosome
    rows = recessiveStats::gencode[recessiveStats::gencode$gene == hgnc &
        recessiveStats::gencode$chr == chrom, ]
    stopifnot(nrow(rows) == 1)
    
    start=rows$start
    end=rows$stop
    
    # extract variants within the region from the VCF
    vars = seqminer::readVCFToListByRange(fileName=vcf_path,
        range=paste(chrom, ":", start, "-", end, sep=""),
        annoType="",
        vcfColumn=c("CHROM", "POS", "REF", "ALT"),
        vcfInfo=c(),
        vcfIndv=c("GT"))
    
    vep = get_ddd_vep_annotations(chrom, start, end)
    
}
