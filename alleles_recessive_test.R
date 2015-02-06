

#' test for enrichment of inherited variants
#'
#' @param hgnc HGNC symbol for a gene.
#' @param chrom chromosome that the gene is on.
#' @param lof_lof number of inherited Lof/LoF events.
#' @param lof_func number of inherited Lof/Func events.
#' @export
#'
#' @return a vector of P values from tests using the DDD population, the ExAC
#'     population, under LoF and functional tests.
analyse_inherited_enrichment <- function(hgnc, chrom, lof_lof, lof_func) {
    ddd = get_ddd_variants_for_gene(hgnc, chrom)
    exac = get_exac_variants_for_gene(hgnc, chrom)
    
    ddd = get_cumulative_frequencies(ddd)
    exac = get_cumulative_frequencies(exac)
    
    ddd$lof_rate = ddd$lof * ddd$lof
    ddd$lof_func_rate = ddd$lof * ddd$functional
    exac$lof_rate = exac$lof * ddd$lof
    exac$lof_func_rate = exac$lof * exac$functional
    
    # get the probability of o getting more than or equal to the number of
    # observed inherited events
    ddd_lof_p = pbinom(lof_lof - 1, 4297, prob=ddd$lof_rate, lower.tail=FALSE)
    exac_lof_p = pbinom(lof_lof - 1, 4297, prob=exac$lof_rate, lower.tail=FALSE)
    ddd_func_p = pbinom(lof_func - 1, 4297, prob=ddd$lof_func_rate, lower.tail=FALSE)
    exac_func_p = pbinom(lof_func - 1, 4297, prob=exac$lof_func_rate, lower.tail=FALSE)
    
    p_values = c(ddd_lof_p, exac_lof_p, ddd_func_p, exac_func_p)
    
    return(p_values)
}
