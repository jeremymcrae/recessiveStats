

#' get and parse ExAC variants within a gene
#'
#' @param hgnc HGNC symbol for the gene that we want data for eg "ARID1B"
#' @param chrom chromosome string e.g. "6"
#' @export
#'
#' @return data frame of variants in gene
get_exac_variants_for_gene <- function(hgnc, chrom) {
    
    # find the gene coordinates, sometimes we have multiple genes with the same
    # HGNC symbol on different chromosomes, so we also need to restrict by
    # chromosome
    rows = recessiveStats::gencode[recessiveStats::gencode$gene == hgnc &
        recessiveStats::gencode$chr == chrom, ]
    stopifnot(nrow(rows) == 1)
    
    start=rows$start
    end=rows$stop
    
    # extract variants within the region from the VCF
    vars = seqminer::readVCFToListByRange(fileName=EXAC_PATH,
        range=paste(chrom, ":", start, "-", end, sep=""),
        annoType="",
        vcfColumn=c("CHROM", "POS", "REF", "ALT"),
        vcfInfo=c("AC", "AN", "CSQ"),
        vcfIndv=c())
    
    vars$sampleId = NULL
    
    vars = data.frame(vars, stringsAsFactors=FALSE)
    
    vars$CQ = apply(vars, 1, parse_vep_output, hgnc=hgnc)
    vars$CSQ = NULL
    
    vars = standardise_multiple_alt_variants(vars)
    
    return(vars)
}
