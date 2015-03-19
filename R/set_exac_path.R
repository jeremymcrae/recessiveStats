

exacPathEnv = new.env()
assign("EXAC", "/lustre/scratch113/resources/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz", envir=exacPathEnv)


#' sets the path to the ExAC dataset files
#'
#' @param path path to the ExAC VCF.
#' @export
set_exac_path <- function(path) {
    assign("EXAC", path, envir=exacPathEnv)
}
