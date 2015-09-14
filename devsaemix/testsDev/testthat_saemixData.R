library(saemix)
library(testthat)
context("Creating a SaemixData object")

test_that("Errors in creating a SaemixData object", {
  x1<-saemixData()

  expect_that(x1, equals("Creation of saemixData failed"))
})


test_that("Successful creation of a SaemixData object", {
  data(theo.saemix)
  x1<-saemixData(name.data=theo.saemix,header=T, name.predictors=c("Dose","Time"),verbose=F)
  
  expect_that(class(x1)[[1]], equals("SaemixData")) 
  expect_that(length(x1@name.predictors), equals(2))
})
