# unit testing for the recessiveStats functions

library(recessiveStats)
library(testthat)

context("VEP parsing checks")

test_that("parse_vep_output output is correct", {
    
    # check a variant where none of the alleleic or transcript consequences
    # are for the required gene
    hgnc = "SAMD11"
    csq = paste("A|ENSG00000187634|ENST00000420190|Transcript|",
        "splice_donor_variant||||||rs375781587|1||1|SAMD11|HGNC|28706|",
        "protein_coding|||ENSP00000411579||Q5SV95_HUMAN&I7FV93_HUMAN&",
        "A6PWC8_HUMAN|UPI000155D47C||||3/6||ENST00000420190.1:c.254+1G>A|||||||",
        "A:0.000228|A:0||||||||INTRON_SIZE:702|||HC", sep="")
    variant = list("CSQ"=csq, "ALT"="A")
    expect_identical(parse_vep_output(variant, hgnc), "splice_donor_variant")
    
    # define a blank consequence vector of the same length as actual VEP output
    # vectors. That way we can simply use the blank vector, and alter the fields
    # we want, rather than continually defining the long VEP strings
    csq = rep("", 49)
    
    # check a variant with a single allele in a single transcript matching the
    # desired gene
    transcript = csq
    transcript[c(1,5,15)] = c("A", "splice_donor_variant", "SAMD11")
    variant = list("CSQ"=paste(transcript, collapse="|"), "ALT"="A")
    expect_identical(parse_vep_output(variant, hgnc), "splice_donor_variant")
    
    # check a variant where none of the alleleic or transcript consequences
    # are for the required gene
    hgnc = "TEST"
    variant = list("CSQ"=paste(transcript, collapse="|"), "ALT"="A")
    expect_identical(parse_vep_output(variant, hgnc), NA)
    
    # check a variant with multiple conseuqnces in a single transcript
    hgnc = "SAMD11"
    transcript[c(1,5,15)] = c("A", "splice_region_variant&synonymous_variant", "SAMD11")
    variant = list("CSQ"=paste(transcript, collapse="|"), "ALT"="A")
    expect_identical(parse_vep_output(variant, hgnc), "splice_region_variant")
    
    # check a variant with mutlple consequences across difference transcripts
    # for the same gene (but where they are all for the same allele).
    transcript_1 = csq
    transcript_2 = csq
    transcript_1[c(1,5,15)] = c("A", "missense_variant", "SAMD11")
    transcript_2[c(1,5,15)] = c("A", "splice_region_variant&synonymous_variant", "SAMD11")
    variant = list("CSQ"=paste(paste(transcript_1, collapse="|"),
        paste(transcript_2, collapse="|"), sep=","), "ALT"="A")
    expect_identical(parse_vep_output(variant, hgnc), "missense_variant")
    
    # check a variant with multple alleles
    transcript_1[c(1,5,15)] = c("A", "missense_variant", "SAMD11")
    transcript_2[c(1,5,15)] = c("G", "splice_region_variant&synonymous_variant", "SAMD11")
    variant = list("CSQ"=paste(paste(transcript_1, collapse="|"),
        paste(transcript_2, collapse="|"), sep=","), "ALT"="A,G")
    expect_identical(parse_vep_output(variant, hgnc), "missense_variant,splice_region_variant")
    
    # check that the consequences order matches the ALT allele order
    variant = list("CSQ"=paste(paste(transcript_1, collapse="|"),
        paste(transcript_2, collapse="|"), sep=","), "ALT"="G,A")
    expect_identical(parse_vep_output(variant, hgnc), "splice_region_variant,missense_variant")
})

test_that("get_most_severe_consequence output is correct", {
    # check two LoF consequences
    cq = c("transcript_ablation", "splice_donor_variant")
    expect_identical(get_most_severe_consequence(cq), "transcript_ablation")
    
    # check a Lof and a functional consequence
    cq = c("splice_donor_variant", "inframe_deletion")
    expect_identical(get_most_severe_consequence(cq), "splice_donor_variant")
    
    # check two functional consequences
    cq = c("missense_variant", "inframe_deletion")
    expect_identical(get_most_severe_consequence(cq), "inframe_deletion")
    
    # check an empty vector, although I don't think we should end up here, since
    # have a different check for empty vectors in he function that calls
    # get_most_severe_consequence
    cq = c()
    expect_identical(get_most_severe_consequence(cq), NA)
})
