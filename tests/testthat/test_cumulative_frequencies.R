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
        AC AN CQ
        1 1000  missense_variant
        1 1000  stop_gained
        1 1000  stop_lost
        1 1000  synonymous_variant
        ")
    
    expect_equal(get_cumulative_frequencies(vars),
        list(lof=0.001, functional=0.002, synonymous=0.001))
})

test_that("correct cumulative frequencies with higher MAF", {
    vars = read.table(header = TRUE, text = "
        AC AN CQ
        1 1000  missense_variant
        1 1000  stop_gained
        20 1000  stop_lost
        1 1000  synonymous_variant
        ")
    
    expect_equal(get_cumulative_frequencies(vars),
        list(lof=0.001, functional=0.001, synonymous=0.001))
})

test_that("correct cumulative frequencies without rare functional variants", {
    vars = read.table(header = TRUE, text = "
        AC AN CQ
        20 800  missense_variant
        1 500  stop_gained
        20 800  stop_lost
        1 1000  synonymous_variant
        ")
    
    # if we don't have any rare functional vars, the number is determined from
    # the lowest allele number (plus two).
    expect_equal(get_cumulative_frequencies(vars),
        list(lof=0.002, functional=0.00132978723404255, synonymous=0.001))
})

test_that("correct cumulative frequencies when we lack a rare variants", {
    vars = read.table(header = TRUE, text = "
        AC AN CQ
        20 1000  missense_variant
        20 500  stop_gained
        20 1000  stop_lost
        20 1000  synonymous_variant
        ")
    
    # if we don't have any rare variants, then we get NA values
    expect_equal(get_cumulative_frequencies(vars),
        list(lof=NA, functional=NA, synonymous=NA))
})

test_that("correct cumulative frequencies when we test a list of dataframes", {
    vars1 = read.table(header = TRUE, text = "
        AC AN CQ
        1 1000  missense_variant
        1 1000  stop_gained
        ")
    
    vars2 = read.table(header = TRUE, text = "
        AC AN CQ
        1 1000  missense_variant
        1 1000  stop_gained
        ")
    vars = list("first"=vars1, "second"=vars2)
    
    expect_equal(get_cumulative_frequencies(vars),
        list("first"=list(lof=0.001, functional=0.001, synonymous=0.000998004),
            "second"=list(lof=0.001, functional=0.001, synonymous=0.000998004)))
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
        list(lof=0.009, functional=0.001, synonymous=0.001))
})
