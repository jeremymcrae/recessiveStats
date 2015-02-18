

#' modifies the consequence where the variant is at the last base in the exon
#'
#' @param variant row of dataframe
#' @param exon_ends vector of exon end positions for the variants gene.
#' @export
#'
#' @return dataframe where the consequence has been modified if appropriate.
check_for_last_base_in_exon <- function(variant, exon_ends) {
    
    allowed = c("missense_variant", "synonymous_variant")
    
    # ignore variants where the consequence isn't missense or synonymous
    if (!variant[["CQ"]] %in% allowed) { return(variant) }
    
    # we don't modify variants that don't have G as reference allele
    if (variant[["REF"]] != "G") { return(variant) }
    
    if (variant[["POS"]] %in% exon_ends) {
        variant[["CQ"]] = "splice_donor_variant"
        print(paste(variant[["CHROM"]], ":", variant[["POS"]], " is non-LoF ref allele G at last base of exon", sep=""))
    }
    
    return(variant)
}

#' gets all exon end positions for an HGNC symbol
#'
#' @param hgnc_symbol HGNC symbol for gene of interest.
#' @export
#'
#' @return list of vectors of exon end chromosome positions.
get_exon_ends <- function(hgnc_symbol) {
    base_url = "http://grch37.rest.ensembl.org/"
    
    gene_ids = get_gene_ids_for_hgnc(base_url, hgnc_symbol)
    transcript_ids = get_transcript_ids_for_ensembl_gene_ids(base_url, gene_ids, hgnc_symbol)
    
    all_ends = c()
    for (transcript_id in transcript_ids) {
        ext = paste("/overlap/id/", transcript_id, "?feature=exon", sep="")
        url = paste(base_url, ext, sep="")
        json = request_from_ensembl(url)
        json = rjson::fromJSON(json)
        
        # restrict to the exons for the transcript
        exons = json[sapply(json, function(x) x[["Parent"]] == transcript_id)]
        
        # get the final position of each exon
        ends = sapply(exons, function(x) x[["end"]])
        if (exons[[1]]$strand == -1) {
            ends = sapply(exons, function(x) x[["start"]])
        }
        
        all_ends = c(all_ends, ends)
    }
    
    # we only need a unique set of exon end positions, since we simply match
    # against these positions
    all_ends = unique(all_ends)
    
    return(all_ends)
}

#' gets the gene IDs for HGNC symbols
#'
#' @param base_url url for REST server.
#' @param hgnc_symbol HGNC symbol for gene of interest.
#' @export
#'
#' @return vector of ensembl gene IDs.
get_gene_ids_for_hgnc <- function(base_url, hgnc_symbol) {
    ext = paste("/xrefs/symbol/homo_sapiens/", hgnc_symbol, sep="")
    url = paste(base_url, ext, sep="")
    
    json = request_from_ensembl(url)
    json = rjson::fromJSON(json)
    
    genes = json[sapply(json, function(x) x[["type"]] == "gene")]
    gene_ids = sapply(genes, function(x) x[["id"]])
    
    return(gene_ids)
}

#' gets the transcript IDs for ensembl gene IDs
#'
#' @param base_url url for REST server.
#' @param gene_ids vector of gene IDs.
#' @param hgnc HGNC symbol for gene of interest.
#' @export
#'
#' @return vector of transcript IDs.
get_transcript_ids_for_ensembl_gene_ids <- function(base_url, gene_ids, hgnc) {
    chroms = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
        "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
        "X", "Y")
    hgnc_symbols = hgnc
    
    transcript_ids = c()
    for (gene_id in gene_ids) {
        ext = paste("/overlap/id/", gene_id, "?feature=transcript", sep="")
        url = paste(base_url, ext, sep="")
        json = request_from_ensembl(url)
        json = rjson::fromJSON(json)
        
        for (item in json) {
            # ignore transcripts not on the standard chromosomes
            # (non-default chroms fail to map the known de novo variants
            # to the gene location
            if (item[["Parent"]] != gene_id |
                !(item[["seq_region_name"]] %in% chroms) |
                all(!grepl(paste(hgnc_symbols, collapse="|"), item[["external_name"]]))) {
                next
            }
            
            transcript_ids = c(transcript_ids, item[["id"]])
        }
    }
    
    return(transcript_ids)
}

timeEnv = new.env()
assign("initial_time", Sys.time(), envir=timeEnv)

#' make a URL request to the Ensembl service
#'
#' @param url string of URL to be accessed
#' @param tries number of attempts that have been mader to access the URL
#'
#' @export
#' @return a character string, typically json encoded
request_from_ensembl <- function(url, tries=0) {
    
    # strip any possible whitespace from the url, since imported datasets
    # occasionally contain whitespace in the values used to construct the URL.
    url = gsub(" ", "", url)
    
    # check that we are not requesting urls fater than that allowed by Ensembl,
    # sleep until the period expires
    current_time = Sys.time()
    diff = 0.067 - as.numeric(current_time - get("initial_time", current_time, envir=timeEnv))
    if (diff > 0) { Sys.sleep(diff) }
    
    # set the previous request time to that of the current request
    assign("initial_time", current_time, envir=timeEnv)
    
    # cut out after making 5 attempts to access the url
    tries = tries + 1
    stopifnot(tries <= 5)
    
    request = httr::GET(url)
    
    # handle the possible http request return status codes, such as when the
    # server is unavailable, when we have made too many requests, or requested
    # an impossible URL.
    if (request$status_code == 503) { # server down
        Sys.sleep(30)
        return(request_from_ensembl(url, tries))
    } else if (request$status_code == 429) { # too frequent requests
        reset_time = as.numeric(request$headers$`x-ratelimit-reset`)
        Sys.sleep(reset_time)
        return(request_from_ensembl(url, tries))
    } else if (request$status_code == 400) { # bad url
        msg = paste("bad url request: ", url, sep = "")
        stop(msg)
    } else if (request$status_code != 200) { # other status errors
        msg = paste("unknown status: ", request$status$code, " at: ", url, sep = "")
        stop(msg)
    }
    
    return(intToUtf8(request$content))
}
