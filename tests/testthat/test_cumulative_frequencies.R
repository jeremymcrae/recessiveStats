# unit testing for the recessiveStats functions

library(recessiveStats)
library(testthat)

context("Cumulative frequency checks")

test_that("get_cumulative_frequencies output is correct for simplest table", {
    vars = read.table(header = TRUE, text = "
        AC AN CQ
        1 1000  missense_variant
        1 1000  stop_gained
        ")
    
    expect_equal(get_cumulative_frequencies(vars),
        list(lof=0.001, functional=0.001, synonymous=0.000998004))
})

test_that("correct cumulative frequencies with larger table", {
    vars = read.table(header = TRUE, text = "
        CHROM  POS  REF  ALT  AC AN CQ
        1      1    A    G    1 1000  missense_variant
        1      2    G    C    1 1000  stop_gained
        1      3    T    A    1 1000  stop_lost
        1      4    G    T    1 1000  synonymous_variant
        ")
    
    expect_equal(get_cumulative_frequencies(vars),
        list(lof=0.001, functional=0.001999, synonymous=0.001))
})

test_that("correct cumulative frequencies with higher MAF", {
    vars = read.table(header = TRUE, text = "
        CHROM  POS  REF  ALT  AC AN    CQ
        1      1    A    G    1  1000  missense_variant
        1      2    G    C    1  1000  stop_gained
        1      3    T    A    20 1000  stop_lost
        1      4    G    T    1  1000  synonymous_variant
        ")
    
    expect_equal(get_cumulative_frequencies(vars),
        list(lof=0.001, functional=0.001, synonymous=0.001))
})

test_that("correct cumulative frequencies without rare functional variants", {
    vars = read.table(header = TRUE, text = "
        CHROM  POS  REF  ALT  AC AN    CQ
        1      1    A    G    20 1000  missense_variant
        1      2    G    C    1  500   stop_gained
        1      3    T    A    20 1000  stop_lost
        1      4    G    T    1  1000  synonymous_variant
        ")
    
    # if we don't have any rare functional vars, the number is determined from
    # the lowest allele number (plus two).
    expect_equal(get_cumulative_frequencies(vars),
        list(lof=0.002, functional=0.00132978723404255, synonymous=0.001))
})

test_that("correct cumulative frequencies when we lack a rare variants", {
    vars = read.table(header = TRUE, text = "
        CHROM  POS  REF  ALT  AC AN    CQ
        1      1    A    G    20 1000  missense_variant
        1      2    G    C    20 500   stop_gained
        1      3    T    A    20 1000  stop_lost
        1      4    G    T    20 1000  synonymous_variant
        ")
    
    # if we don't have any rare variants, then we get NA values
    expect_equal(get_cumulative_frequencies(vars),
        list(lof=NA, functional=NA, synonymous=NA))
})

test_that("correct cumulative frequencies when we test a list of dataframes", {
    vars1 = read.table(header = TRUE, text = "
        CHROM  POS  REF  ALT  AC AN   CQ
        1      1    A    G    1 1000  missense_variant
        1      2    G    C    1 1000  stop_gained
        1      3    G    T    1 1000  synonymous_variant
        ")
    
    vars2 = read.table(header = TRUE, text = "
        CHROM  POS  REF  ALT  AC AN   CQ
        1      1    A    G    1 1000  missense_variant
        1      2    G    C    1 1000  stop_gained
        1      3    G    T    1 1000  synonymous_variant
        ")
    vars = list("first"=vars1, "second"=vars2)
    
    expect_equal(get_cumulative_frequencies(vars),
        list("first"=list(lof=0.001, functional=0.001, synonymous=0.001),
            "second"=list(lof=0.001, functional=0.001, synonymous=0.001)))
})

test_that("correct cumulative frequencies when one pop is above the threshold", {
    vars1 = read.table(header = TRUE, text = "
        CHROM  POS  REF  ALT  AC AN   CQ
        1      1    A    G    1 1000  missense_variant
        1      2    A    G    1 1000  missense_variant
        1      3    G    C    1 1000  stop_gained
        1      4    G    T    1 1000  synonymous_variant
        ")
    
    vars2 = read.table(header = TRUE, text = "
        CHROM  POS  REF  ALT  AC AN    CQ
        1      1    A    G    1  1000  missense_variant
        1      2    A    G    50 1000  missense_variant
        1      3    G    C    1  1000  stop_gained
        1      4    G    T    1  1000  synonymous_variant
        ")
    vars = list("first"=vars1, "second"=vars2)
    
    expect_equal(get_cumulative_frequencies(vars),
        list("first"=list(lof=0.001, functional=0.001, synonymous=0.001),
            "second"=list(lof=0.001, functional=0.001, synonymous=0.001)))
})

test_that("correct cumulative frequencies when we change the frequency threshold", {
    vars = read.table(header = TRUE, text = "
        AC AN CQ
        1 1000  missense_variant
        1 1000  stop_gained
        8 1000  stop_gained
        1 1000  synonymous_variant
        ")
    
    expect_equal(get_cumulative_frequencies(vars, threshold=0.005),
        list(lof=0.001, functional=0.001, synonymous=0.001))
    
    # expect higher frequencies when we use a higher threshold
    expect_equal(get_cumulative_frequencies(vars, threshold=0.01),
        list(lof=0.008992, functional=0.001, synonymous=0.001))
})

test_that("remove_high_frequency_vars is correct for easiest cases", {
    vars1 = read.table(header = TRUE, text = "
        CHROM  POS  REF  ALT  AC AN   CQ
        1      1    A    G    1 1000  missense_variant
        1      2    A    G    1 1000  missense_variant
        1      3    G    C    1 1000  stop_gained
        1      4    G    T    1 1000  synonymous_variant
        ")
    
    vars2 = read.table(header = TRUE, text = "
        CHROM  POS  REF  ALT  AC AN    CQ
        1      1    A    G    1  1000  missense_variant
        1      2    A    G    10 1000  missense_variant
        1      3    G    C    1  1000  stop_gained
        1      4    G    T    1  1000  synonymous_variant
        ")
    
    expect_equal(remove_high_frequency_vars(vars1, 0.01), vars1)
    expect_equal(remove_high_frequency_vars(vars2, 0.01), vars2[-2, ])
    
    vars = list("first"=vars1, "second"=vars2)
    expected = list("first"=vars1, "second"=vars2)
    expect_equal(remove_high_frequency_vars(vars, 0.02), expected)
})

test_that("remove_high_frequency_vars is correct when one pop is above the threshold", {
    vars1 = read.table(header = TRUE, text = "
        CHROM  POS  REF  ALT  AC AN   CQ
        1      1    A    G    1 1000  missense_variant
        1      2    A    G    1 1000  missense_variant
        1      3    G    C    1 1000  stop_gained
        1      4    G    T    1 1000  synonymous_variant
        ")
    
    vars2 = read.table(header = TRUE, text = "
        CHROM  POS  REF  ALT  AC AN    CQ
        1      1    A    G    1  1000  missense_variant
        1      2    A    G    9 1000  missense_variant
        1      3    G    C    1  1000  stop_gained
        1      4    G    T    1  1000  synonymous_variant
        ")
    vars = list("first"=vars1, "second"=vars2)
    
    expected = list("first"=vars1[-2, ], "second"=vars2[-2, ])
    
    expect_equal(remove_high_frequency_vars(vars, 0.005),
        expected)
})


test_that("correct cumulative frequency caclulation", {
    
    freqs = rep(0.005, 500)
    expect_equal(cumulative_frequency(freqs), 0.9184281386)
    
    # corect expectation for empty vector
    freqs = c()
    expect_equal(cumulative_frequency(freqs), 0)
    
    # correct expectation for zero value
    freqs = 0
    expect_equal(cumulative_frequency(freqs), 0)
    
    # correct expectation for vector with NA
    freqs = c(NA, 0.01)
    expect_equal(cumulative_frequency(freqs), 0.01)
    
    # correct expectation for vector with NA only
    freqs = c(NA)
    expect_equal(cumulative_frequency(freqs), 0)
})
