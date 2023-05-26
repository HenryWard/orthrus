context("Inner scoring functions")
library(orthrus)

test_that("fix_residual_length works on valid input", {
  a <- c(1,2)
  b <- c(3,4)
  expanded <- expand.grid(a, b)
  result <- expanded[1:2,]
  expect_equal(fix_residual_length(expanded, 2), result)
})

test_that("fix_residual_length works on valid expanding input", {
  a <- c(1,2)
  b <- c(3,4)
  expanded <- expand.grid(a, b)
  to_append <- data.frame(matrix(NA, nrow = 1, ncol = 2))
  result <- rbind(expanded, to_append)
  expect_equal(fix_residual_length(expanded, 5), result)
})

test_that("fix_residual_length works for no expected change", {
  a <- c(1,2)
  b <- c(3,4)
  expanded <- expand.grid(a, b)
  result <- expanded
  expect_equal(fix_residual_length(expanded, 4), result)
})