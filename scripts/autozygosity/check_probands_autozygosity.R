#### script to identify ROH regions for all probands

library(recessiveStats)
library(argparse)

BCF_PATH = "/lustre/scratch113/projects/ddd/users/jm33/ddd_4k.bcftools.bcf"
OUTPUT_DIR = "data-raw/autozygosity"
SCRIPT_PATH = "scripts/autozygosity/proband_autozygosity.R"
RSCRIPT_BINARY = "/software/R-3.2.2/bin/Rscript"

get_options <- function() {
    parser = ArgumentParser()
    parser$add_argument("--bcf", default=BCF_PATH, help="Path to bcf to analyse.")
    parser$add_argument("--diagnosed",
        help="Include path to table listing probands with diagnoses, if you want to exclude these probands.")
    parser$add_argument("--script", default=SCRIPT_PATH, help="Path to R script to run.")
    parser$add_argument("--rbinary", default=RSCRIPT_BINARY, help="Path to Rscript binary.")
    parser$add_argument("--output-folder", default=OUTPUT_DIR, help="Path to put output files into.")

    args = parser$parse_args()
    
    return(args)
}

get_bjobs <- function() {
    args = c("-o", "\"JOBID", "USER", "STAT", "QUEUE", "JOB_NAME", "delimiter=';'\"")
    bjobs = system2("bjobs", args, stdout=TRUE)
    
    # if there aren't any running jobs, return a blank dataframe
    if (length(bjobs) == 1) {
        return(data.frame(JOBID=character(0), USER=character(0), STAT=character(0),
            QUEUE=character(0), JOB_NAME=character(0)))
    }
    
    bjobs = read.table(text=bjobs, sep=";", header=TRUE)
    
    return(bjobs)
}

get_running_jobs <- function() {
    bjobs = get_bjobs()
    return(bjobs = bjobs[bjobs$STAT == "RUN", ])
}

get_proband_ids <- function() {
    ddd = recessiveStats::get_ddd_cohort(parents=FALSE, unaffected=FALSE)
    probands = ddd[ddd$dad_id != 0 & ddd$mum_id != 0, ]
    
    return(sort(unique(probands$individual_id)))
}

get_undiagnosed_sanger_ids <- function(diagnosed_path) {
    probands = get_proband_ids()
    
    # remove the probands with diagnoses, or likely diagnostic variants.
    diagnosed = read.table(diagnosed_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    probands = probands[!probands %in% diagnosed$person_id]
    
    return(probands)
}

main <- function() {
    
    args = get_options()
    
    if (is.null(args$diagnosed)) {
        probands = get_proband_ids()
    } else {
        probands = get_undiagnosed_sanger_ids(args$diagnosed)
    }
    
    for (proband in probands) {
        
        while (nrow(get_bjobs()) > 100) { Sys.sleep(30) }
        
        output_path = file.path(args$output_folder, proband)
        
        command = "bsub"
        arguments = c("-o", "get_autozygosity.bjob",
            "bash", "-c",
            "\"", args$rbinary, args$script,
            "--proband", proband,
            "--bcf", args$bcf,
            "--output", output_path,
            "\"")
        
        system2(command, arguments)
        Sys.sleep(0.25)
    }
}

main()
