# functions to assess autozygosity rates within a gene
# Currently, these functions depend heavily upon:
#     a) bcftools (version 1.0)
#     b) files for DDD probands (and a file listing the undiagnosed probands)
# It should be straightforward to adapt the code to define a different set of
# probands to analyse, but the bcftools requirement is much stricter. You could
# define the path to the bcftools executable data-raw/constants.R, rerun the
# constants.R code, then reinstall the package.

#' identify DDD probands to analyse
#'
#' @param chrom chromosome of interest (this will only be important if we have
#'        different sets of probands for each chromosome file, which shouldn't
#'        be the case, but we might as well be sure.)
#' @param vcf_dir directory containing multi-sample VCFs per-chromosome for the
#'        population of interest.
#' @export
#'
#' @examples
#' get_undiagnosed_sanger_ids("14")
get_undiagnosed_sanger_ids <- function(chrom, vcf_dir=DDD_VCFS_DIR) {
    # get the probands from the DDD study
    ddd = get_ddd_cohort(parents=FALSE, unaffected=FALSE)
    ddd = ddd[!duplicated(ddd), ]
    ddd = ddd[ddd$dad_id != 0, ]
    
    # remove the diagnosed probands
    diagnosed = "/lustre/scratch113/projects/ddd/users/jm33/ddd_likely_diagnosed.txt"
    diagnosed = read.table(diagnosed, header=TRUE)
    ddd = ddd[!ddd$individual_id %in% diagnosed$person_id, ]
    
    # if the gene is on chrX, we can only estimate autozygosity for females, 
    # since males are hemizygous for chrX.
    if (chrom == "X") { ddd = ddd[ddd$sex == "F", ] }
    
    # identify the chromosome-specific VCF to extract from
    vcf_path = Sys.glob(file.path(vcf_dir, paste(chrom, "\\:1-*.vcf.gz", sep="")))
    
    # extract the sample IDs from the multi-sample VCF
    command = BCFTOOLS
    args = c("query",
        "-l", vcf_path)
    samples = system2(command, args, stdout=TRUE)
    
    proband_ids = samples[samples %in% ddd$sanger_id]
    
    return(proband_ids)
}

#' create a BCF for a small genome region, for a limited number of probands
#'
#' @param chrom chromosome that we wish to analyse
#' @param start_pos nucleotide position that starts the genome region
#' @param end_pos nucleotide position that end the genome region
#' @param vcf_dir directory containing multi-sample VCFs per-chromosome for the
#'        population of interest.
#' @export
#'
#' @return path to the newly created BCF
extract_region <- function(chrom, start_pos, end_pos, vcf_dir) {
    # identify the VCF to extract from. The chromosome specific VCF should be
    # named as: CHROM:1-CHROM_BP_END.vcf.gz e.g. 1:1-249250621.vcf.gz
    vcf_path = Sys.glob(file.path(vcf_dir, paste(chrom, "\\:1-*.vcf.gz", sep="")))
    
    # identify the samples to select, and where to put the resulting file
    probands = get_undiagnosed_sanger_ids(chrom, vcf_dir)
    temp_bcf = tempfile(fileext=".bcf")
    
    command = BCFTOOLS
    args = c("view",
        "--samples", paste(probands, collapse=","),
        "--regions", paste(chrom, ":", start_pos, "-", end_pos, sep=""),
        "--output-type", "b",
        "--output-file", temp_bcf,
        vcf_path)
    
    system2(command, args)
    
    return(temp_bcf)
}

#' identify whether a proband is autozygous within a genome region
#'
#' @param bcf_path path to a BCF for the genome region
#' @param chrom chromosome that we wish to analyse
#' @param start_pos nucleotide position that starts the genome region
#' @param end_pos nucleotide position that end the genome region
#' @param proband ID for the proband to analyse
#' @export
#'
#' @return True/False for whether the sample is autozygous within the region
check_sample_autozygosity <- function(bcf_path, chrom, start_pos, end_pos, proband) {
    
    command = BCFTOOLS
    args = c("roh",
        "--biallelic-sites",
        "--estimate-AF", "all",
        "--GTs-only", 30,
        "--skip-indels",
        "--samples", proband,
        bcf_path)
    
    # # for bcftools version 1.2, which for an undetermined reason does not
    # # find a known ROH region for a control sample. I suspect I need to either
    # # rejig the command/options, or amend the VCF.
    # command = "/software/hgi/pkglocal/bcftools-1.2/bin/bcftools"
    # args = c("roh",
    #     "--sample", proband,
    #     "--estimate-AF", "-",
    #     "--skip-indels",
    #     "--GTs-only", 30,
    #     bcf_path)
    
    # analyse ROH and get results for the proband
    roh_output = system2(command, args, stdout=TRUE, stderr=FALSE)
    roh = read.table(text=roh_output, sep="\t")
    names(roh) = c("sample_id", "chrom", "pos", "p_value", "state")
    
    # # for bcftools version 1.2
    # names(roh) = c("chrom", "pos", "state", "quality")
    
    # subset down to the variants that lie within the required region
    roh = roh[roh$pos >= start_pos & roh$pos <= end_pos, ]
    
    # find whether any site shows evidence of autozygosity
    autozygous = any(roh$state != 0)
    
    return(autozygous)
}

#' identify whether a proband is autozygous within a genome region
#'
#' @param hgnc HGNC symbol for the gene of interest
#' @param chrom chromosome that we wish to analyse
#' @param vcf_dir directory containing multi-sample VCFs per-chromosome for the
#'        population of interest.
#' @export
#'
#' @return True/False for whether the sample is autozygous within the region
get_autozygous_rate <- function(hgnc, chrom, vcf_dir=DDD_VCFS_DIR) {
    rows = get_gene_coordinates(hgnc, chrom)
    
    start=rows$start
    end=rows$stop
    offset = 2000000
    
    probands = get_undiagnosed_sanger_ids(chrom, vcf_dir)
    temp_bcf = extract_region(chrom, max(start - offset, 1), end + offset, vcf_dir)
    
    autozygous = sapply(probands, function(x) check_sample_autozygosity(temp_bcf, chrom, start, end, x))
    autozygous_rate = sum(autozygous)/length(autozygous)
    
    unlink(temp_bcf)
    
    return(autozygous_rate)
}
