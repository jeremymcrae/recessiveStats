

library(recessiveStats)

genes = recessiveStats::gencode$gene[recessiveStats::gencode$gene_type == "protein_coding"]

get_bjobs <- function() {
    bjobs = system2("bjobs", stdout=TRUE)
    
    # if there aren't any running jobs, return a blank dataframe
    if (length(bjobs) == 1) {
        return(data.frame(
            JOBID=character(0),
            USER=character(0),
            STAT=character(0),
            QUEUE=character(0),
            FROM_HOST=character(0),
            EXEC_HOST=character(0),
            JOB_NAME=character(0),
            SUBMIT_TIME=character(0)))
    }
    
    tmp = tempfile()
    cat(file=tmp, bjobs, sep = "\n")
    bjobs = read.fwf(file=tmp, widths=c(8, 8, 6, 11, 12, 12, 11, 12), strip.white=TRUE)
    unlink(tmp)
    
    # rename the dataframe, since the file doesn't work if we try to use the
    # header
    names(bjobs) = as.vector(unlist(bjobs[1, ]))
    bjobs = bjobs[-1, ]
    
    return(bjobs)
}

get_running_jobs <- function() {
    bjobs = get_bjobs()
    return(bjobs = bjobs[bjobs$STAT == "RUN", ])
}

all_genes = list()
for (hgnc in sort(genes)) {
    chrom = recessiveStats::gencode$chr[recessiveStats::gencode$gene == hgnc &
        recessiveStats::gencode$gene_type == "protein_coding"]
    
    while (nrow(get_bjobs()) > 500) { Sys.sleep(30) }
    
    command = "bsub"
    args = c("-o", "./results/autozygosity.all_genes.consanguinous_probands.bjob_output.txt",
        "-q", "normal",
        "-R", "\"select[mem>100] rusage[mem=100]\"",
        "-M", 100,
        "bash", "-c",
        "\"/software/R-3.1.2/bin/Rscript", "./scripts/autozygosity/gene_autozygosity.R",
        "--hgnc", hgnc,
        "--chrom", chrom,
        "--consang",
        "\"")
    
    system2(command, args)
    Sys.sleep(0.5)
}
