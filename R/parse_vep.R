
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
#' @param variant single row dataframe containing allele column and VEP prediction
#'     column.
#' @param hgnc HGNC symbol for the gene that we are currently investigating.
#' @export
#'
#' @return vep consequence string for the most severe consequence for the gene
#'     of interest. If none of the predicted consequences lie within the gene of
#'    interest, return "NA". If the variant has multiple alleles, return a
#'    comma-separated list of VEP consequences for the alleles.
parse_vep_output <- function(variant, hgnc) {
    
    # Sometimes the VEP prediction is a comma-separated list of VEP predictions.
    # This can occur for two reasons, when we have predictions for the different
    # transcripts that the variant lies within, and secondly, if the variant has
    # multiple alleles, as we need sepaarte predicitons for all the different
    # alleles.
    vep_predictions = unlist(strsplit(variant[["CSQ"]], ","))
    
    # make sure we have a consequence entry for each allele. If none of the
    # predictions for an allele match the required HGNC symbol, then we will
    # have a zero length vector, which should give NA when we look for the most
    # severe consequence.
    consequences = list()
    alleles = unlist(strsplit(variant[["ALT"]], ","))
    
    for (prediction in vep_predictions) {
        prediction = unlist(strsplit(prediction, "\\|"))
        
        allele = prediction[1]
        symbol = prediction[15]
        cq = prediction[5]
        
        if (symbol != hgnc) { next }
        
        cq = unlist(strsplit(cq, "&"))
        consequences[[allele]] = c(consequences[[allele]], cq)
    }
    
    # sort the consequences by the alt alleles from the VCF field
    # TODO: some indels cometimes have different allele names, and different
    # allele orders to the alt allele codes. We should fix these, but some are
    # not simple.
    if (length(consequences) > 1) {
        if (all(sort(names(consequences)) == sort(alleles))) {
            consequences = consequences[alleles]
        }
    }
    
    # if we haven't found any consequences, retun NA
    if (length(consequences) == 0) { return(NA) }
    
    # get the most severe consequence for each allele
    cq = paste(sapply(consequences, get_most_severe_consequence), collapse=",")
    
    return(cq)
}

#' get the most severe VEP consequence form a list of consequences
#'
#' @param consequences vector of VEP consequence strings. Can be an empty vector
#' @export
#'
#' @return most severe VEP consequence string, or NA if the consequences vector
#'     is empty.
get_most_severe_consequence <- function(consequences) {
    
    # convert the consequence strings to their rank (in terms of severity)
    ranks = severity$rank[severity$consequence %in% consequences]
    
    if (length(ranks) == 0) { return(NA) }
    
    most_severe = severity$consequence[severity$rank == min(ranks)]
    
    return(most_severe)
}
