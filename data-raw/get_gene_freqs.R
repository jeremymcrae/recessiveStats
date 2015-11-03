# extract the ExAC cumulative LoF and functional frequencies for a single gene

# # Get frequencies for all genes by running this as a LSF job array:
# bsub \
#   -J freqs[1-20350]%100 \
#   -o results/frequencies/gene.%I.bjob \
#   -R "select[mem>15000] rusage[mem=15000]" \
#   -M 15000 \
#   bash -c "/software/R-3.2.2/bin/Rscript scripts/get_gene_freqs.R \
#      --position \$LSB_JOBINDEX \
#      --output results/frequencies/gene.\$LSB_JOBINDEX.txt"

# # After the jobs complete, create a single file with all gene data:
# head -n 1 results/frequencies/gene.1.txt > results/exac_frequencies.txt
# tail -n +2 -q results/frequencies/gene.*.txt | sort >> results/exac_frequencies.txt
# rm results/frequencies/gene.*.txt

library(recessiveStats)
library(argparse)

get_options <- function() {
    parser = ArgumentParser()
    parser$add_argument("--position", type="integer",
        help="row number of gene to check (from 1 to ~20400).")
    parser$add_argument("--output", help="where to put the output.")
    args = parser$parse_args()
    
    return(args)
}

main <- function() {
    args = get_options()
    
    # restrict outselves to the ~20000 protein coding genes
    genes = recessiveStats::gencode[recessiveStats::gencode$gene_type == "protein_coding", ]
    
    # find which gene we want to get the frequencies for
    hgnc = sort(unique(genes$gene))[args$position]
    chrom = genes$chr[genes$gene == hgnc]
    
    # don't process genes with NA (which could happen if we pass a gene position
    # that doesn't fit within the gene dataframe).
    if (is.na(hgnc)) { stop(0) }
    
    # extract the variants for the ExAC populations
    exac = try(get_exac_variants_for_gene(hgnc, chrom), silent=TRUE)
    
    exac_freqs = list()
    if (class(exac) != "try-error") {
        exac_freqs = get_cumulative_frequencies(exac)
    }
    
    # make sure we include all the populations from ExAC
    missing = list(lof=NA, functional=NA, synonymous=NA)
    for (pop in c("AFR", "AMR", "EAS", "FIN", "NFE", "SAS")) {
        if (!pop %in% names(exac_freqs)) {
            exac_freqs[[pop]] = missing
        }
    }
    
    # write the data to a file
    freqs = data.frame(hgnc=hgnc, chrom=chrom, exac_freqs)
    write.table(freqs, file=args$output, sep="\t", row.names=FALSE, quote=FALSE)
}

main()
