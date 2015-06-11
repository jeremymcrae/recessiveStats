# unit testing for the recessiveStats functions

library(recessiveStats)
library(testthat)

context("Tidy variants checks")

test_that("standardise_multiple_alt_variants output is correct", {
    
    # check that a variant without any alt alleles is processed appropriately
    vars = read.table(header=TRUE, text="
        CHROM  POS  REF  ALT  AC  AN  CQ
        1      1    A    G    1   5   stop_lost",
        colClasses=c("character","numeric","character","character","character",
            "character","character"))
    
    result = read.table(header=TRUE, text="
        CHROM  POS  REF  ALT  AC  AN  CQ
        1      1    A    G    1   5   stop_lost",
        colClasses=c("character","numeric","character","character","numeric",
            "numeric","character"))
    expect_identical(standardise_multiple_alt_variants(vars), result)
    
    # check that a variant with multiple alt alleles is split out correctly
    vars = read.table(header=TRUE, text="
        CHROM  POS  REF  ALT   AC   AN  CQ
        1      1    A    G,C   1,2  5   stop_lost,missense_variant",
        colClasses=c("character","numeric","character","character","character",
            "character","character"))
    
    result = read.table(header=TRUE, text="
        CHROM  POS  REF  ALT  AC  AN  CQ
        1      1    A    G    1   5   stop_lost
        1      1    A    C    2   5   missense_variant",
        colClasses=c("character","numeric","character","character","numeric",
            "numeric","character"))
    expect_identical(standardise_multiple_alt_variants(vars), result)
    
    # check that HGNC symbols are handled correctly without the include_hgnc flag
    vars = read.table(header=TRUE, text="
        CHROM  POS  REF  ALT   AC    AN  CQ                          HGNC
        1      1    A    G     1     5   missense_variant            TEST
        1      2    A    G,C   1,2   5   stop_lost,missense_variant  TEST",
        colClasses=c("character","numeric","character","character","character",
            "character","character"))
    
    result = read.table(header=TRUE, text="
        CHROM  POS  REF  ALT  AC  AN  CQ                HGNC
        1      1    A    G    1   5   missense_variant  TEST
        1      2    A    G    1   5   stop_lost         TEST
        1      2    A    C    2   5   missense_variant  TEST",
        colClasses=c("character","numeric","character","character","numeric",
            "numeric","character"))
    vars = standardise_multiple_alt_variants(vars)
    row.names(vars) = 1:3 # make sure the row names match the expectation
    expect_identical(vars, result)
    
    # check that that HGNC symbols are handled correctly when using the
    # include_hgnc flag
    vars = read.table(header=TRUE, text="
        CHROM  POS  REF  ALT   AC    AN  CQ                          HGNC
        1      1    A    G     1     5   missense_variant            TEST
        1      2    A    G,C   1,2   5   stop_lost,missense_variant  TEST,TEST2",
        colClasses=c("character","numeric","character","character","character",
            "character","character"))
    
    result = read.table(header=TRUE, text="
        CHROM  POS  REF  ALT  AC  AN  CQ                HGNC
        1      1    A    G    1   5   missense_variant  TEST
        1      2    A    G    1   5   stop_lost         TEST
        1      2    A    C    2   5   missense_variant  TEST2",
        colClasses=c("character","numeric","character","character","numeric",
            "numeric","character"))
    vars = standardise_multiple_alt_variants(vars, include_hgnc=TRUE)
    row.names(vars) = 1:3 # make sure the row names match the expectation
    expect_identical(vars, result)
})

test_that("remove_nonfunctional_variants output is correct", {
    
    # make an easy table with one variant that will not be excluded
    vars = read.table(header=TRUE, text="
        CHROM POS REF ALT CQ
        1     1   A   G   stop_lost", stringsAsFactors=FALSE)
    expect_identical(remove_nonfunctional_variants(vars), vars)
    
    # now define a table where at least one variant is dropped
    vars = read.table(header=TRUE, text="
        CHROM POS REF ALT CQ
        1     1   A   G   stop_lost
        1     2   A   C   intron_variant", stringsAsFactors=FALSE)
    
    result = read.table(header=TRUE, text="
        CHROM POS REF ALT CQ
        1     1   A   G   stop_lost", stringsAsFactors=FALSE)
    
    expect_identical(remove_nonfunctional_variants(vars), result)
    
    # define a table where we have both LoF and non-LoF fcuntioal variants
    vars = read.table(header=TRUE, text="
        CHROM POS REF ALT CQ
        1     1   A   G   stop_lost
        1     3   A   T   missense_variant
        1     2   A   C   intron_variant", stringsAsFactors=FALSE)
        
    result = read.table(header=TRUE, text="
        CHROM POS REF ALT CQ
        1     1   A   G   stop_lost
        1     3   A   T   missense_variant", stringsAsFactors=FALSE)
    expect_identical(remove_nonfunctional_variants(vars), result)
})
