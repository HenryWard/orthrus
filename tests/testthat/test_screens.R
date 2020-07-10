context("Screen utility functions")
library(chymeraR)

test_that("add_screen fails on invalid input", {
  expect_error(add_screen(name = NULL, replicates = c("Test1", "Test2")))
  expect_error(add_screen(name = NA, replicates = c("Test1", "Test2")))
  expect_error(add_screen(name = 3, replicates = c("Test1", "Test2")))
  expect_error(add_screen(name = "Test", replicates = c(NA, "Test2")))
  expect_error(add_screen(name = "Test", replicates = NULL))
  expect_error(add_screen(name = "Test", replicates = NA))
  expect_error(add_screen(name = "Test", replicates = 3))
  expect_error(add_screen(name = "Test", replicates = c("Test1", "Test2"),
                          target_coverage = NULL))
  expect_error(add_screen(name = "Test", replicates = c("Test1", "Test2"),
                          target_coverage = NA))
  expect_error(add_screen(name = "Test", replicates = c("Test1", "Test2"),
                          target_coverage = "hello"))
})

test_that("add_screen works on valid input", {
  reps <- c("Test1", "Test2")
  screen <- list()
  screen_list <- list()
  screen[["replicates"]] <- reps
  screen[["target_coverage"]] <- 200
  screen_list[["Test"]] <- screen
  expect_equal(add_screen(name = "Test", replicates = reps), screen_list)
})

test_that("remove_screen fails on invalid input", {
  reps <- c("Test1", "Test2")
  screen_list <- add_screen(name = "Test", replicates = reps)
  expect_error(remove_screen(screen_list, NULL))
  expect_error(remove_screen(screen_list, NA))
  expect_error(remove_screen(screen_list, c("Test", "Hello")))
  expect_warning(remove_screen(screen_list, "NotHere"))
})

test_that("remove_screen works on valid input", {
  reps <- c("Test1", "Test2")
  screen_list <- add_screen(name = "Test", replicates = reps)
  expect_equal(remove_screen(screen_list, "Test"), list())
})


