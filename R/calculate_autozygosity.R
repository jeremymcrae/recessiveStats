# functions to assess autozygosity rates within a gene
# Currently, these functions depend heavily upon:
#     a) bcftools (version 1.0)
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
#' @export
#'
#' @return dataframe of proband ID, chromosome, start and end coordinates.
check_sample_autozygosity_genome_wide <- function(bcf_path, proband) {
    command = BCFTOOLS
    args = c("roh",
        "--biallelic-sites",
        "--estimate-AF", "all",
        "--GTs-only", 30,
        "--skip-indels",
        "--samples", proband,
        bcf_path)
    
    # analyse ROH and get results for the proband
    roh_output = system2(command, args, stdout=TRUE, stderr=FALSE)
    roh = read.table(text=roh_output, sep="\t")
    names(roh) = c("sample_id", "chrom", "pos", "p_value", "state")
    
    # if there aren't any ROH regions predicted in the output, then return a
    # blank dataframe
    if (!any(roh$state == 1)) {
        return(data.frame(sample_id=character(0), chrom=character(0),
            start_pos=numeric(0),end_pos=numeric(0)))
    }
    
    # identify the start and end positions of ROH ranges
    row_positions = which(roh$state == 1)
    ends = which(diff(row_positions) > 1)
    starts = c(row_positions[1], row_positions[ends + 1])
    ends = c(row_positions[ends], row_positions[length(row_positions)])
    
    # define a table of ROH ranges
    coords = data.frame(sample_id=roh$sample_id[starts],
        chrom=roh$chrom[starts], start_pos=roh$pos[starts],
        end_pos=roh$pos[ends])
    
    return(coords)
}
