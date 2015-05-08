

#' get and parse ExAC variants within a gene
#'
#' @param hgnc HGNC symbol for the gene that we want data for eg "ARID1B"
#' @param chrom chromosome string e.g. "6"
#' @param check_last_base boolean for whether to check if last base non-lofs can
#'        be LoF.
#' @export
#'
#' @return data frame of variants in gene
get_exac_variants_for_gene <- function(hgnc, chrom, check_last_base=FALSE) {
    
    rows = get_gene_coordinates(hgnc, chrom)
    start=rows$start
    end=rows$stop
    
    if (!file.exists(get("EXAC", envir=exacPathEnv))) {
        stop(paste("Cannot find the ExAC file. Check that the file exists at: ",
            get("EXAC", envir=exacPathEnv), ", or obtain the ExAC datasets, ",
            "then run set_exac_path(YOUR_PATH) before trying this function ",
            "again.", sep=""))
    }
    
    # define the populations of interest within the ExAC dataset
    populations = list(AFR="African/African American", AMR="American",
        EAS="East Asian", FIN="Finnish", NFE="Non-Finnish European",
        SAS="South Asian")
    
    # extract variants within the region from the VCF
    vars = seqminer::readVCFToListByRange(fileName=get("EXAC", envir=exacPathEnv),
        range=paste(chrom, ":", start, "-", end, sep=""),
        annoType="",
        vcfColumn=c("CHROM", "POS", "REF", "ALT"),
        vcfInfo=c("AC", "AN", paste("AC_", names(populations), sep=""),
            paste("AN_", names(populations), sep=""), "CSQ"),
        vcfIndv=c())
    
    # convert the vars list to a data frame, then extract the necessary info
    # from the VEP data field.
    vars$sampleId = NULL
    vars = data.frame(vars, stringsAsFactors=FALSE)
    vars$CQ = apply(vars, 1, parse_vep_output, hgnc=hgnc)
    vars$CSQ = NULL
    
    if (check_last_base) {
        ends = get_exon_ends(hgnc)
        exon_ends = ends$all_ends
        strand = ends$strand
        vars$CQ = apply(vars, 1, check_for_last_base_in_exon, exon_ends=exon_ends, strand=strand)
    }
    
    vars = get_list_of_population_frequencies(vars, populations)
    
    return(groups)
}

#' create a list of allele frequency dataframes for each population.
#'
#' @param vars data frame of variants in gene, with allele counts and totals for
#'        several ExAC sub-groups.
#' @param populations list of populations (e.g. list(EAS="East Asian",
#'        FIN="Finnish")). Perhaps this should be a simple vector of names, but
#'        this way I have more information about what the populations are.
#' @export
#'
#' @return list of dataframes, one for each population. The dataframes contain
#'         allele count, allele number and VEP consequence strings for each
#'         allele in the input dataframe. The list is named for each of the
#'         populations.
get_list_of_population_frequencies <- function(vars, populations) {
    
    # define a blank list to place the formatted dataframes into
    groups = vector(mode="list", length=length(populations))
    names(groups) = names(populations)
    
    for (population in names(populations)) {
        # define the allele count and total alleles for the population
        vars$AC = vars[[paste("AC_", population, sep="")]]
        vars$AN = vars[[paste("AN_", population, sep="")]]
        
        # split out multiple alleles, drop all unecessary columns, and add the
        # dataframe to the list
        by_allele = standardise_multiple_alt_variants(vars)
        by_allele = by_allele[, c("CHROM", "POS", "REF", "ALT", "AC", "AN", "CQ")]
        groups[[population]] = by_allele
    }
    
    return(groups)
}
