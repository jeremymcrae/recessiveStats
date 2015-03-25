
#' get the genomic coordinates for a gene from an internal gencode dataframe
#'
#' @param hgnc HGNC symbol for the gene that we want data for eg "ARID1B"
#' @param chrom chromosome string e.g. "6"
#' @export
#'
#' @return data frame of genomic cooridantes
get_gene_coordinates <- function(hgnc, chrom) {
    # find the gene coordinates, sometimes we have multiple genes with the same
    # HGNC symbol on different chromosomes, so we also need to restrict by
    # chromosome
    rows = recessiveStats::gencode[recessiveStats::gencode$gene == hgnc &
        recessiveStats::gencode$chr == chrom, ]
    
    # give reasonable error messages if we cannot find a single gene range
    if (nrow(rows) == 0) {
        stop(paste("Unable to find a gene range for ", hgnc, " on chrom ",
            chrom, ". Check that the gene is on the correct chromosome.", sep=""))
    } else if (nrow(rows) > 1) {
        stop(paste("Multiple gene ranges for ", hgnc, " on chrom ",
            chrom, ". Check that the gene symbol?", sep=""))
    }
    
    return(rows)
}
