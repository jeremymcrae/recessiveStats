
#' tidy up a variant dataset, so that multiple alt variants are in standardised rows
#'
#' Make a new dataframe for the multiple alt variants, where each row is
#' for a single alt allele, so alleles at a single site span multiple rows.
#'
#' @param vars dataframe of variants, with consequence, alt allele, and allele
#'     count columns.
#' @param include_hgnc whether the dataframe also includes a HGNC column that
#'     needs to be standardised.
#' @export
#'
#' @return dataframe of variants with mutlple alt alleles spread across separate
#'     rows.
standardise_multiple_alt_variants <- function(vars, include_hgnc=FALSE) {
    single_vars = vars[!(grepl(",", vars$CQ)), ]
    multi_vars = vars[grepl(",", vars$CQ), ]
    
    # spread the multiple alts into individual rows
    new_vars = multi_vars[0, ]
    if (nrow(multi_vars) > 0) {
        for (i in 1:nrow(multi_vars)) {
            row = multi_vars[i, ]
            alts = unlist(strsplit(row$ALT, ","))
            
            for (j in 1:length(alts)) {
                new_row = row
                new_row$ALT = alts[j]
                new_row$CQ = unlist(strsplit(row$CQ, ","))[j]
                new_row$AC = unlist(strsplit(row$AC, ","))[j]
                
                if (include_hgnc) {
                    new_row$HGNC = unlist(strsplit(row$HGNC, ","))[j]
                }
                
                new_vars = rbind(new_vars, new_row)
            }
        }
    }
    
    vars = rbind(single_vars, new_vars)
    
    # make sure the allele counts and totals are numeric, so that later
    # calculations work seamlessly
    vars$AN = as.numeric(as.character(vars$AN))
    vars$AC = as.numeric(as.character(vars$AC))
    
    return(vars)
}

#' exclude nonfunctional variants, or nonfunctional alleles of variants
#'
#' @param vars dataframe of variants, with consequence, alt allele, and allele
#'     count columns.
#' @export
#'
#' @return dataframe of variants with nonfunctional variants/alleles trimmed
remove_nonfunctional_variants <- function(vars) {
    
    # define the VEP consequence types for loss of function and missense variants
    lof_cq = c("stop_gained", "splice_acceptor_variant", "splice_donor_variant",
        "frameshift_variant", "transcript_ablation")
    
    functional_cq = c("stop_lost", "initiator_codon_variant",
        "transcript_amplification", "inframe_insertion", "inframe_deletion",
           "missense_variant", "coding_sequence_variant")
    
    vars = vars[vars$CQ %in% c(lof_cq, functional_cq), ]
    
    return(vars)
}
