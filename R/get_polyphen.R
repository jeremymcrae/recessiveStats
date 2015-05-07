
#' get a polyphen prediction for a variant
#'
#' @param chrom chromosome that the variant is on (e.g. "8")
#' @param start start nucleotide position of the variant
#' @param end end nucleotide position of the variant
#' @param alt alternate allele code (e.g. "C")
#' @export
#'
#' @return the polyphen prediction for a variant.
get_polyphen <- function(chrom, start, end, alt) {

    # construct the VEP REST url for the variant
    base_url = "http://grch37.rest.ensembl.org/"
    ext = "/vep/human/region/"
    ext = paste(ext, chrom, ":", start, ":", end, "/", alt, sep="")
    url = paste(base_url, ext)

    # get the VEP data for the variant
    request = request_from_ensembl(url)
    request = rjson::fromJSON(request)[[1]]

    score = 0
    term = NA
    for (transcript in request$transcript_consequences) {
        # ignore non-missense consequences, since they won't have a polyphen prediction
        if (transcript$consequence_terms != "missense_variant") { next }

        if ("polyphen_prediction" %in% names(transcript)) {
            if (transcript$polyphen_score > score) {
                score = transcript$polyphen_score
                term = transcript$polyphen_prediction
            }
        }
    }

    return(term)
}
