# functions to assess autozygosity rates within a gene
# Currently, these functions depend heavily upon:
#     a) bcftools (version 1.3)
#     b) VCF (actually, BCF) file for probands used to discover biallelically
#        inherited variants.
# It should be straightforward to adapt the code to define a different set of
# probands to analyse, but the bcftools requirement is much stricter. You could
# define the path to the bcftools executable data-raw/constants.R, rerun the
# constants.R code, then reinstall the package.

#' identify the genome regions where a proband has runs-of-homozygosity
#'
#' @param bcf_path path to a BCF for the full genome across all probands
#' @param proband ID for the proband to analyse
#' @param chrom chromosome to analyse
#' @param map genetic recombination map for the given chromosome
#' @export
#'
#' @return dataframe of proband ID, chromosome, start and end coordinates.
check_sample_autozygosity_genome_wide <- function(bcf_path, proband, chrom, map) {
    
    coords = NULL
    for (chrom in as.character(seq(1, 22))) {
        chrom_map = gsub('XXX', chrom, map)
        
        chrom_coords = check_autozygosity_in_chrom(bcf_path, proband, chrom, map)
        
        if (is.null(coords)) {
            coords = chrom_coords
        } else {
            coords = rbind(coords, chrom_coords)
        }
    }
    
    return(coords)
}

#' identify chromosomal regions where a proband has runs-of-homozygosity
#'
#' @param bcf_path path to a BCF for the full genome across all probands
#' @param proband ID for the proband to analyse
#' @param chrom chromosome to analyse
#' @param map genetic recombination map for the given chromosome
#' @export
#'
#' @return dataframe of proband ID, chromosome, start and end coordinates.
check_autozygosity_in_chrom <- function(bcf_path, proband, chrom, map) {
    command = BCFTOOLS
    args = c("roh",
        "--estimate-AF", "-",
        "--GTs-only", 30,
        "--skip-indels",
        "--regions", chrom,
        "--genetic-map", map,
        "--sample", proband,
        bcf_path)
    
    # analyse ROH and get results for the proband
    roh_output = system2(command, args, stdout=TRUE, stderr=FALSE)
    roh = try(utils::read.table(text=roh_output, sep="\t"), silent=TRUE)
    
    # if the command doesn't return any output (such as for probands who only
    # have CNV calls, which would not make it through to the merged VCF, which
    # only contains biallelic SNVs/indels)
    if (class(roh) == "try-error") {
        roh = data.frame("chrom"=character(0), "pos"=numeric(0),
            "state"=numeric(0), "score"=numeric(0))
    }
    
    names(roh) = c("chrom", "pos", "state", "score")
    
    # if there aren't any ROH regions predicted in the output, then return a
    # blank dataframe
    if (!any(roh$state == 1)) {
        return(data.frame("chrom"=character(0), "start_pos"=numeric(0),
            "end_pos"=numeric(0)))
    }
    
    # identify the start and end positions of ROH ranges
    row_positions = which(roh$state == 1)
    ends = which(diff(row_positions) > 1)
    starts = c(row_positions[1], row_positions[ends + 1])
    ends = c(row_positions[ends], row_positions[length(row_positions)])
    
    # define a table of ROH ranges
    coords = data.frame("chrom"=roh$chrom[starts], "start_pos"=roh$pos[starts],
        "end_pos"=roh$pos[ends])
    
    return(coords)
}
