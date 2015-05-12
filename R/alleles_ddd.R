
#' get the details of the DDD cohort
#'
#' @param parents boolean for whether to restrict to parents only.
#' @param unaffected boolean for whether to restrict to unaffected only.
#' @export
#'
#' @return vector of paths to VCFs
get_ddd_cohort <- function(parents=TRUE, unaffected=TRUE) {
    ped_path = file.path(DATAFREEZE_DIR, "family_relationships.txt")
    sanger_id_path = file.path(DATAFREEZE_DIR, "person_sanger_decipher.txt")
    
    ped = read.table(ped_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    sanger_ids = read.table(sanger_id_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    if (parents) { ped = ped[ped$dad_id == 0, ] }
    if (unaffected) { ped = ped[ped$affected == 1, ] }
    
    ped = merge(ped, sanger_ids, by.x="individual_id", by.y="person_stable_id", all.x=TRUE)
    
    return(ped)
}

#' find where the required data lies within the file
#'
#' Seek in a file to the section containing the data that we want. We start
#' with a reasonable guess about where in the file the data is located, but
#' we seek up or downstream of that (depending on whether the data at that
#' position is before or after the data that we want), then bisect to pin
#' down where to start loading data from.
#'
#' @param path path to the datafile
#' @param start_pos start position of the variant
#' @export
#' @return the approximate byte offset in the file (to within a few hundred
#'    bytes) where the required data starts
find_file_pos <- function(path, start_pos) {
    
    # seek to the correct part of the file
    file_con = file(path, open="r")
    
    # make sure we don't use negative seek positions, or go beyond the file size
    lower = 0
    upper = file.info(path)$size
    
    while (upper - lower > 100) {
        mid_point = lower + (upper - lower)/2
        prev_pos = seek(file_con, mid_point, origin="start")
        
        read = scan(file_con, what=character(), skip=1, sep="\t", nmax=2, quiet=TRUE)
        # if we have gone back into the VCF heaeader, then shift the mid point
        # so we get out of it again
        if (grepl("^#", read[1])) {upper = mid_point; next}
        temp_start = as.numeric(read[2])
        
        # narrow the boundaries to half the previous size
        if (temp_start > start_pos) {upper = mid_point}
        if (temp_start <= start_pos) {lower = mid_point}
    }
    
    close(file_con)
    
    # make sure the seek position is within the file
    seek_position = max(lower, 0)
    seek_position = min(seek_position, file.info(path)$size)
    
    return(seek_position)
}

#' loads DDD VEP annotations lines for variants
#'
#' @param chrom chromosome as character string
#' @param start start of gene range
#' @param end end of gene range
#' @export
#'
#' @return dataframe of VEP annotations
get_lines_from_vep_vcf <- function(chrom, start, end) {
    vep_dir = "/lustre/scratch113/projects/ddd/users/ddd/test_msVCF_annotate_vep/vep_annotated_vcfs/results_vcfs"
    vep_path = file.path(vep_dir, paste(chrom, ".txt", sep=""))
    
    # get the byte offsets for where the data lies within the file (this is a
    # workaround to speed access to the data, while not having tabix indexes).
    file_start = max(find_file_pos(vep_path, start) - 10000, 0)
    file_end = min(find_file_pos(vep_path, end) + 10000, file.info(vep_path)$size)
    
    # load the lines from the file (assume the lines are about 400 bytes long)
    line_length = 400
    con = file(vep_path, open="r")
    seek(con, file_start, origin="start")
    lines_n = (file_end - file_start) / line_length
    vep = read.table(con, skip=1, sep="\t", nrows=lines_n, stringsAsFactors=FALSE)
    names(vep) = c("chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format")
    close(con)
    
    # make sure we have captured variants that lie outside the gene sites,
    # otherwise we don't know if we have missed any variants inside the gene range
    stopifnot(vep$pos[1] < start)
    stopifnot(vep$pos[nrow(vep)] > end)
    
    return(vep)
}

#' gets the DDD VEP annotations for variants
#'
#' @param chrom chromosome as character string
#' @param start start of gene range
#' @param end end of gene range
#' @export
#'
#' @return dataframe of VEP annotations
get_ddd_vep_annotations <- function(chrom, start, end) {
    
    vep = get_lines_from_vep_vcf(chrom, start, end)
    vep = vep[vep$pos > start & vep$pos < end, ]
    
    vep$HGNC = NA
    vep$CQ = NA
    for (pos in 1:nrow(vep)) {
        row = vep[pos, ]
        info = unlist(strsplit(row$info, ";"))
        # some rows lack a "HGNC" field, but do have a "HGNC_ALL" field, so we
        # shall use this instead
        if (!any(grepl("HGNC=", info))) { info = gsub("HGNC_ALL", "HGNC", info) }
        
        # make sure there is a HGNC field (I'm not sure if this will be used,
        # since we are looking at variants within gene regions)
        if (!any(grepl("HGNC=", info))) { info = c(info, "HGNC=") }
        
        # get the HGNC and consequence fields from the INFO
        vep$HGNC[pos] = unlist(strsplit(info[grepl("HGNC=", info)], "="))[2]
        vep$CQ[pos] = unlist(strsplit(info[grepl("CQ=", info)], "="))[2]
    }
    
    # trim the unnecesary columns, it'll make a cleaner merge later
    vep$info = NULL
    vep$qual = NULL
    vep$filter = NULL
    vep$format = NULL
    vep$id = NULL
    
    return(vep)
}

#' reformat to hemizygous chrX genotypes in a single column
#'
#' The column has previously been vetted as for a male, on chromosome X, not
#' within a pseudoautosomal region.
#'
#' @param column column of genotypes
#' @export
#'
#' @return column of variants, where male allosomal chrX genotypes have been
#'     converted to hemizygous genotypes (e.g. "0/0" -> "0").
get_hemizygous <- function(column) {
    
    # recode the NA values to avoid splitting or excluding NAs later
    column[column == "./."] = "a/b"
    
    alleles = strsplit(column, "/")
    first = sapply(alleles, "[", 1)
    second = sapply(alleles, "[", 2)
    
    # Male chrX genotypes are assigned one of their alleles. We pick the first
    # when both are identical, since that will convert "0/0" to "0" and "1/1" to
    # "1", which suits the later purpose of getting the frequency of the alleles
    # in the population. Reinsert the NA values.
    column[first == second] = first[first == second]
    column[column == "a/b"] = NA
    
    return(column)
}

#' reformat chrX genotypes for males to hemizygous
#'
#' Female chrX genotypes are untouched, as are genotypes on autosomal
#' chromosomes. Males with heterozygous chrX genotypes are left unchanged, as we
#' cannot be sure which genotype is correct. Possibly they should be removed. We
#' also do not alter genotypes within pseduo-autosomal regions.
#'
#' @param vars seqminer::readVCFToListByRange output, list of values, which
#'     includes a genotype matrix as "GT", and sample IDs as sampleId.
#' @param geno dataframe of biallelic genotypes, one column per individual
#' @param ddd_parents information about DDD individuals, including sex
#' @export
#'
#' @return dataframe of variants, where male allosomal chrX genotypes have been
#'     converted to hemizygous genotypes (e.g. "0/0" -> "0").
reformat_chrX_genotypes <- function(vars, geno, ddd_parents) {
    
    # define the pseudoautosomal regions
    x_par = data.frame(matrix(c(60001, 2699520, 154930290, 155260560,
        88456802, 92375509), ncol=2, byrow=TRUE))
    names(x_par) = c("start", "end")
    
    # don't alter non-chrX variants
    if (vars$CHROM[1] != "X") {
        geno[geno == "./."] = NA
        return(geno)
    }
    
    # don't alter genes overlapping the pseudoautosomal regions
    start_pos = vars$POS[1]
    end_pos = vars$POS[length(vars$POS)]
    if (any(start_pos < x_par$end & end_pos > x_par$start)) { return(geno) }
    
    # convert biallelic male genotypes on chrX to hemizygous genotypes
    is_male = names(geno) %in% ddd_parents$sanger_id[ddd_parents$sex == "M"]
    geno[, is_male] = apply(geno[, is_male], 2, get_hemizygous)
    
    # convert the female null genotype codes to NA, since it's faster to do the
    # males separately. This seems to be awkwardly constructed, perhaps there
    # is a simpler way to express this.
    female = geno[, !is_male]
    female[female == "./."] = NA
    geno[, !is_male] = female
    
    return(geno)
}

#' count alleles for a variant
#'
#' @param variant single row of data frame of VCF genotypes
#' @export
#'
#' @return dataframe of variants
tally_alleles <- function(variant) {
    alt_column = which(names(variant) == "ALT")
    alts = unlist(strsplit(variant[alt_column], ","))
    
    # drop the alt column, so the variant only contains genotype entries
    variant = variant[-alt_column]
    
    # drop NA genotypes, and split the genotypes into alleles, then count the
    # different alleles used in the genotypes
    variant = variant[!is.na(variant)]
    variant = unlist(strsplit(variant, "/"))
    counts = table(variant)
    
    # get the non-ref allele counts
    alt_positions = names(counts)[names(counts) != 0]
    all_alts = 1:length(alts)
    missing = all_alts[!(all_alts %in% alt_positions)]
    non_ref = counts[names(counts) != 0]
    
    if (length(missing) > 0) {
        new = rep(0, length(missing))
        names(new) = missing
        non_ref = c(non_ref, new)
        non_ref = non_ref[order(names(non_ref))]
    }
    
    values = c(paste(non_ref, collapse=","), sum(counts))
    
    return(values)
}


#' quickly strips out the variants which are obviously nonfunctional
#'
#' @param vars seqminer::readVCFToListByRange output, list of values, which
#'     includes a genotype matrix as "GT", and sample IDs as sampleId.
#' @param vep dataframe of VEP consequences for the vars entry
#' @export
#'
#' @return dataframe of variants
remove_easy_nonfunctional <- function(vars, vep) {
    
    nonfunc_cq = c("intron_variant", "5_prime_UTR_variant", "3_prime_UTR_variant")
    
    nonfunc = vep$CQ %in% nonfunc_cq
    
    # remove all the nonfunctional variants from the sections of the vars list
    vars$CHROM = vars$CHROM[!nonfunc]
    vars$POS = vars$POS[!nonfunc]
    vars$REF = vars$REF[!nonfunc]
    vars$ALT = vars$ALT[!nonfunc]
    vars$GT = vars$GT[, !nonfunc]
    
    return(vars)
}

#' reformats the genotypes dataset
#'
#' @param vars seqminer::readVCFToListByRange output, list of values, which
#'     includes a genotype matrix as "GT", and sample IDs as sampleId.
#' @param vep dataframe of VEP consequences for the vars entry
#' @param probands vector of proband IDs, or NA, or a zero-length vector.
#' @export
#'
#' @return dataframe of variants
convert_genotypes <- function(vars, vep, probands) {
    vars = remove_easy_nonfunctional(vars, vep)
    
    cat("converting DDD genotypes\n")
    # convert the genotype matrix to a dataframe where each column is for a
    # separate individual
    geno = data.frame(t(vars$GT), stringsAsFactors=FALSE)
    names(geno) = vars$sampleId
    
    # restrict the genotypes to individuals who are unaffected parents
    ddd = get_ddd_cohort(parents=FALSE, unaffected=FALSE)
    ddd_parents = ddd[ddd$dad_id == 0 & ddd$affected == 1, ]
    geno = geno[, names(geno) %in% ddd_parents$sanger_id]
    
    # remove any parents of the probands, since they could skew the allele
    # frequencies
    if (!is.null(probands)) {
        proband_parents = ddd[ddd$individual_id %in% probands, ]
        geno = geno[, !(names(geno) %in% proband_parents$sanger_id)]
    }
    
    geno = reformat_chrX_genotypes(vars, geno, ddd_parents)
    geno$ALT = vars$ALT
    
    cat("tallying DDD allele counts\n")
    counts = apply(geno, 1, tally_alleles)
    
    vars$AC = counts[1, ]
    vars$AN = counts[2, ]
    vars$GT = NULL
    vars$sampleId = NULL
    
    vars = data.frame(vars, stringsAsFactors=FALSE)
    
    return(vars)
}

#' loads the variants for a given gene from source VCFs
#'
#' @param hgnc hgnc symbol as character string
#' @param chrom chromosome as character string
#' @param probands vector of IDs for probands for the gene, or NULL
#' @param check_last_base boolean for whether to check if last base non-lofs can
#'     be LoF.
#' @export
#'
#' @return dataframe of variants
get_ddd_variants_for_gene <- function(hgnc, chrom, probands=NULL, check_last_base=FALSE) {
    
    vcf_path = Sys.glob(file.path(DDD_VCFS_DIR, paste(chrom, "\\:1-*.vcf.gz", sep="")))
    
    rows = get_gene_coordinates(hgnc, chrom)
    
    start=rows$start
    end=rows$stop
    
    cat(paste("loading DDD vcfs for ", hgnc, "\n", sep = ""))
    # extract variants within the region from the VCF
    vars = seqminer::readVCFToListByRange(fileName=vcf_path,
        range=paste(chrom, ":", start, "-", end, sep=""),
        annoType="",
        vcfColumn=c("CHROM", "POS", "REF", "ALT"),
        vcfInfo=c(),
        vcfIndv=c("GT"))
    
    vep = get_ddd_vep_annotations(chrom, start, end)
    vars = convert_genotypes(vars, vep, probands)
    
    vars = merge(vars, vep, by.x=c("CHROM", "POS", "REF", "ALT"),
        by.y=c("chrom", "pos", "ref", "alt"), all.x=TRUE)
        
    vars = standardise_multiple_alt_variants(vars, include_hgnc=TRUE)
    
    if (check_last_base) {
        ends = get_exon_ends(hgnc)
        exon_ends = ends$all_ends
        strand = ends$strand
        vars$CQ = apply(vars, 1, check_for_last_base_in_exon, exon_ends=exon_ends, strand=strand)
    }
    
    # remove alleles with none observed in the unaffected DDD parents, and
    # alleles not in the required gene
    vars = vars[vars$AC > 0 & vars$HGNC == hgnc, ]
    
    # remove the HGNC column, to match the output for the ExAC functions
    vars$HGNC = NULL
        
    return(vars)
}
