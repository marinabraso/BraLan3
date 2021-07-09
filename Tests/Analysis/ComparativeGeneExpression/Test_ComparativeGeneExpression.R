#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")
library(testthat)
script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,6, nchar(script)-12), "_functions.R", sep=""))

test_that("test", {
	a <- 6
	b <- 7

	result <- a*b

	expect_equal(result, 42)
})


test_that("test_get_spearman_corr", {
	tissue <- "brain"
	rnorm <- rnorm(100)
	testpairs <- data.frame(paste0(tissue,"1")=rnorm, paste0(tissue,"2")=rnorm*10)
	print(head(testpairs))

	result <- GetSpearmanCorr(testpairs, tissue)

	expect_equal(result, 1)
})




