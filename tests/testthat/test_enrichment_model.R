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

context("Population combination generation checks")

test_that("get_count_combinations is correct", {
    
    populations = c("A", "B")
    count = 2
    
    result = read.table(header=TRUE, text="
        A B
        2 0
        1 1
        0 2",
        colClasses=c("integer", "integer"))
    
    expect_equal(get_count_combinations(populations, count), result)
    expect_error(get_count_combinations(populations, -1), "count >= 0 is not TRUE")
})

test_that("get_count_combinations is correct for count of zero", {
    
    populations = c("A", "B")
    count = 0
    
    result = read.table(header=TRUE, text="
        A B
        0 0",
        colClasses=c("integer", "integer"))
    
    expect_equal(get_count_combinations(populations, count), result)
})

test_that("get_count_combinations is correct for more cohorts", {
    
    populations = c("A", "B", "C")
    count = 4
    
    result = read.table(header=TRUE, text="
        A B C
        4 0 0
        3 1 0
        2 2 0
        1 3 0
        0 4 0
        3 0 1
        2 1 1
        1 2 1
        2 2 1
        0 3 1
        2 0 2
        1 1 2
        2 1 2
        0 2 2
        1 2 2
        2 2 2
        1 0 3
        0 1 3
        0 0 4",
        colClasses=c("integer", "integer"))
    
    expect_equal(get_count_combinations(populations, count), result)
})

context("Enrichment tests return the correct P-values")

test_that("sum_combo_tests is correct", {
    
    exac = list("A"=list("lof"=0.1, "functional"=0.1),
        "B"=list("lof"=0.1, "functional"=0.1))
    cohort_n = list("A"=100, "B"=100)
    combos = get_count_combinations(names(cohort_n), count=1)
    
    expect_equal(sum_combo_tests(exac, cohort_n, combos,
        biallelic_lof_enrichment), 0.866020325142038)
    
    # and test for a larger number of families
    combos = get_count_combinations(names(cohort_n), count=7)
    expect_equal(sum_combo_tests(exac, cohort_n, combos,
        biallelic_lof_enrichment), 0.00429554232885157)
    
    # check that summing across populations where the frequencies are the same
    # gives the same answer as for the single population p-value, even if the
    # population sizes differ between ethnicities.
    cohort_n = list("A"=20, "B"=180)
    expect_equal(sum_combo_tests(exac, cohort_n, combos, biallelic_lof_enrichment),
        pbinom(6, 200, 0.01, lower.tail=FALSE))
})

test_that("single enrichment tests are correct", {
    
    freq = list("lof"=0.1, "functional"=0.2)
    cohort_n = 100
    count = 2
    
    # make sure the calculation uses the binomial model we expect
    expect_equal(biallelic_lof_enrichment(freq, count, cohort_n),
        pbinom(1, 100, 0.01, lower.tail=FALSE))
    
    # make sure the calculation gives us the exact value we expect
    expect_equal(biallelic_lof_enrichment(freq, count, cohort_n),
        0.264238021077044)
    
    # a count of zero should have a p-value of 1
    expect_equal(biallelic_lof_enrichment(freq, 0, cohort_n), 1)
    
    # check the lof/func enrichment calculation
    expect_equal(lof_func_enrichment(freq, count, cohort_n),
        0.947531944275668)
    
    # check the functional enrichment calculation
    expect_equal(biallelic_func_enrichment(freq, count, cohort_n),
        0.91283668331261)
    
    # check that NA counts return NA values
    expect_identical(biallelic_func_enrichment(freq, NA, cohort_n),
        as.numeric(NA))
})
