# unit testing for the recessiveStats functions

library(recessiveStats)
library(testthat)

context("Recessive statistic checks")

test_that("test_enrichment output is correct", {
    
    freq = list(lof=0.001, functional=0.1)
    
    result = list(lof=0.001, functional=0.1,
        biallelic_lof_p=0.000999500666126,
        lof_func_p=0.001158649871933,
        biallelic_func_p=0.989927345227986)
    
    expect_equal(test_enrichment(freq, 1, 4, 2, 1000), result)
})

test_that("test_enrichment output is correct for NA input", {
    
    freq = list(lof=NA, functional=NA)
    
    result = list(lof=NA, functional=NA,
        biallelic_lof_p=as.numeric(NA),
        lof_func_p=as.numeric(NA),
        biallelic_func_p=as.numeric(NA))
    
    expect_identical(test_enrichment(freq, 1, 4, 2, 1000), result)
})
