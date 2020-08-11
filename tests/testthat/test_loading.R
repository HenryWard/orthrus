context("Data loading functions")
library(orthrus)

# We test a plotting function because it's one of two possible functions called at the
# beginning of the pipeline
test_that("plot_screen_reads fails on invalid input", {
  screens <- add_screen(name = "HAP1_T0", replicates = "DoesNotExist")
  df <- chymera_paralog
  expect_error(plot_screen_reads(df, screens, "output"))
})

test_that("normalize_screens fails on invalid input", {
  screens <- add_screen(name = "HAP1_T0", replicates = "DoesNotExist")
  df <- chymera_paralog
  expect_error(normalize_screens(df, screens))
})

