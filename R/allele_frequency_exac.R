
# consequence list, as sorted at http://www.ensembl.org/info/genome/variation/predicted_data.html
consequences = c("transcript_ablation", "splice_donor_variant",
    "splice_acceptor_variant", "stop_gained", "frameshift_variant", "stop_lost",
    "initiator_codon_variant", "transcript_amplification", "inframe_insertion",
    "inframe_deletion", "missense_variant", "splice_region_variant",
    "incomplete_terminal_codon_variant", "stop_retained_variant",
    "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant",
    "5_prime_UTR_variant", "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant", "intron_variant",
    "NMD_transcript_variant", "non_coding_transcript_variant",
    "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation",
    "TFBS_amplification", "TF_binding_site_variant",
    "regulatory_region_ablation", "regulatory_region_amplification",
    "regulatory_region_variant", "feature_elongation", "feature_truncation",
    "intergenic_variant")
severity = data.frame(consequence=consequences, rank=seq(1:length(consequences)),
stringsAsFactors=FALSE)

#' parses VEP output
#'
#' @param data single row dataframe containing allele column and VEP prediction
#'     column.
#' @param hgnc HGNC symbol for the gene that we are currently investigating.
#' @export
#'
#' @return vep consequence string for the most severe consequence for the gene
#'     of interest. If non of the predicted consequences lie within the gene of
#'    interest, return "NA". If the variant has multiple alleles, return a
#'    comma-separated list of VEP consequences for the alleles.
parse_vep_output <- function(data, hgnc) {
    
    # Sometimes the VEP prediction is a comma-separated list of VEP predictions.
    # This can occur for two reasons, when we have predictions for the different
    # transcripts that the variant lies within, and secondly, if the variant has
    # multiple alleles, as we need sepaarte predicitons for all the different
    # alleles.
    vep_predictions = unlist(strsplit(data[["CSQ"]], ","))
    
    consequences = list()
    for (prediction in vep_predictions) {
        prediction = unlist(strsplit(prediction, "\\|"))
        
        allele = prediction[1]
        symbol = prediction[15]
        cq = prediction[5]
        
        # make sure we have an entry for the allele. If non of the predictions
        # for an allele match the required HGNC symbol, then we will have a zero
        # length vector, which should give NA when we look for the most severe
        # consequence.
        if (!allele %in% names(consequences)) { consequences[[allele]] = c() }
        
        if (symbol != hgnc) { next }
        
        cq = unlist(strsplit(cq, "&"))
        consequences[[allele]] = c(consequences[[allele]], cq)
    }
    
    # get the most severe consequence for each allele
    cq = paste(sapply(consequences, get_most_severe_consequence), collapse=",")
    
    return(cq)
}

#' get the most severe VEP consequence form a list of consequences
#'
#' @param consequences vector of VEP consequence strings.
#' @export
#'
#' @return most severe VEP consequence string, or NA if the consequences vector
#'     has no entries.
get_most_severe_consequence <- function(consequences) {
    
    most_severe = NA
    best_severity = NA
    
    ranks = severity$rank[severity$consequence %in% cq]
    
    for (cq in consequences) {
        
        temp_severity = severity$rank[severity$consequence == cq]
        if (is.na(best_severity) | temp_severity < best_severity) {
            best_severity = temp_severity
            most_severe = cq
        }
    }
    
    return(most_severe)
}

#' tidy up a variant dataset, so that multiple alt variants are in standardised rows
#'
#' Make a new dataframe for the multiple alt variants, where each row is
#' for a single alt allele, so alleles at a single site span multiple rows.
#'
#' @param vars dataframe of variants, with consequence, alt allele, and allele
#'     count columns.
#' @export
#'
#' @return dataframe of variants with mutlple alt alleles spread across separate
#'     rows.
standardise_multiple_alt_variants <- function(vars) {
    single_vars = vars[!(grepl(",", vars$CQ)), ]
    multi_vars = vars[grepl(",", vars$CQ), ]
    
    # spread the multiple alts into individual rows
    new_vars = multi_vars[0, ]
    for (i in 1:nrow(multi_vars)) {
        row = multi_vars[i, ]
        alts = unlist(strsplit(row$ALT, ","))
        
        for (j in 1:length(alts)) {
            new_row = row
            new_row$ALT = alts[j]
            new_row$CQ = unlist(strsplit(row$CQ, ","))[j]
            new_row$AC = unlist(strsplit(row$AC, ","))[j]
            new_vars = rbind(new_vars, new_row)
        }
    }
    
    vars = rbind(single_vars, new_vars)
    
    return(vars)
}

#' exclude nonfunctional variants, or nonfunctional alleles of variants
#'
#' @param vars dataframe of variants, with consequence, alt allele, and allele
#'     count columns.
#' @export
#'
#' @return dataframe of varinats with nonfunctional variants/alleles trimmed
remove_nonfunctional_variants <- function(vars) {
    
    # define the VEP consequence types for loss of function and missense variants
    lof_cq = c("stop_gained", "splice_acceptor_variant", "splice_donor_variant",
        "frameshift_variant")
    missense_cq = c("missense_variant", "initiator_codon_variant", "stop_lost",
        "inframe_deletion", "inframe_insertion", "splice_region_variant")
    
    vars = vars[vars$CQ %in% c(lof_cq, missense_cq), ]
    
    return(vars)
}

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
    vars = remove_nonfunctional_variants(vars)
    
    return(vars)
}
