# prepare the gencode dataset for the package

gencode_path = "data-raw/all_gencode_genes_v19.txt"
gencode = read.table(gencode_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# strip out a few columns, to reduce saved file size
gencode$gene_id = NULL
gencode$transcript_id = NULL

# standardise the chromosome strings
gencode$chr = gsub("chr", "", gencode$chr)

save(gencode, file="data/gencode.rda", compress="xz")
