#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")
library(testthat)
script <- sub(".*=", "", commandArgs()[4])
source(paste0(dirname(script), "/", substr(basename(script),6, nchar(basename(script))-2), "_functions.r", sep=""))

test_that("get_spearman_corr no NA", {
	numbers1 <- c(-5, .9, 0)
	numbers2 <- c(-3, 1.9, 12)
	testpairs <- data.frame("brain1" = numbers1, "brain2" = numbers2)

	result <- get_spearman_corr(testpairs, "brain")

	expect_equal(result, .5)
})

test_that("test_get_spearman_corr with NA", {
	numbers1 <- c(-5, .9, 0, NA, 3, NA)
	numbers2 <- c(-3, 1.9, 12, 6, NA, NA)
	testpairs <- data.frame("brain1" = numbers1, "brain2" = numbers2)

	result <- get_spearman_corr(testpairs, "brain")

	expect_equal(result, .5)
})

test_that("test_get_spearman_corr with all NA column", {
	numbers1 <- c(-5, .9, 0, NA, 3, NA)
	numbers2 <- c(NA, NA, NA, NA, NA, NA)
	testpairs <- data.frame("brain1" = numbers1, "brain2" = numbers2)

	result <- get_spearman_corr(testpairs, "brain")

	expect_equal(is.na(result), TRUE)
})

test_that("test_get_spearman_corr with dataframe with no rows", {
	testpairs <- data.frame("brain1" = c(), "brain2" = c())

	expect_error(get_spearman_corr(testpairs, "brain"))
})

test_that("test_get_spearman_corr with dataframe with no columns", {
	testpairs <- data.frame()

	expect_error(get_spearman_corr(testpairs, "brain"))
})

test_that("test_get_spearman_corr with dataframe with no propperly named columns", {
	numbers1 <- c(-5, .9, 0)
	numbers2 <- c(-3, 1.9, 12)
	testpairs <- data.frame("brain3" = numbers1, "brain5" = numbers2)

	expect_error(get_spearman_corr(testpairs, "brain"))
})

test_that("test_sum_union_of_patterns", {
	tissues <- data.frame("Name" = c("brain", "liver"))
	tpmthresh <- 10
	spnum <- 1
	OG <- "OG_1"
	testpairs <- data.frame(c(1, 30), c(5, 2), c(OG, OG))
	colnames(testpairs) <- c("brain1", "liver1", "OG")

	result <- sum_union_of_patterns(OG, spnum, testpairs, tissues, tpmthresh)

	expect_equal(result, 1)
})

test_that("test_absolute_difference_union_of_patterns", {
	tissues <- data.frame("Name" = c("brain", "liver"))
	tpmthresh <- 10
	OG <- "OG_1"
	testpairs <- data.frame(c(1, 30), c(1, 5), c(60, 80), c(7, 33), c(OG, OG))
	colnames(testpairs) <- c("brain1", "liver1", "brain2", "liver2", "OG")

	result <- absolute_difference_union_of_patterns(OG, testpairs, tissues, tpmthresh)

	expect_equal(result, 1)
})


