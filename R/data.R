#' table of syndromes, with regexes to match in the probands recorded syndromes
#'
#' @format A data frame with two columns variables: \code{regex} and \code{name}
"SYNDROMES"

#' table of genecode gene coordinates
#'
#' Tom Fitzgerald provided this table, where the positions were derived from a
#' GENCODE data release.
#'
#' @source url(http://www.gencodegenes.org/releases/19.html)
#'
#' @format A data frame with two columns variables: \code{chr}, \code{start},
#'     \code{stop}, \code{gene_id}, \code{transcript_id}, \code{gene},
#'     \code{status}, \code{level} and \code{gene_type}.
"gencode"
