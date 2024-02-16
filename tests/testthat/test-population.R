test_that("custom families are named", {
  expect_equal(ols_with_error(rnorm)$family, "ols_with_error")
  expect_equal(custom_family(rnorm, function(x) x)$family, "custom_family")
})

test_that("sampling without error matches linear predictor", {
  toy <- population(
    x = predictor("rnorm"),
    y = response(1 + 2 * x, error_scale = 0)
  )

  d <- toy |>
    sample_x(n = 10) |>
    sample_y()

  expect_equal(nrow(d), 10)
  expect_equal(d$y, 1 + 2 * d$x)
})

test_that("ols_with_error reports errors", {
  bad_err_fn <- function(n) { 1 }

  fam <- ols_with_error(bad_err_fn)

  expect_error(fam$simulate(NULL, 1, data.frame(), 1:10),
               class = "regressinator_error_length")
})

test_that("response errors on missing error_scale", {
  expect_error(response(4 + 2 * x), class = "regressinator_error_scale")
})

test_that("population_predictors gets all predictor names", {
  foo <- population(
    x1 = predictor("rnorm"),
    x2 = predictor("rnorm"),
    y = response(x1 + x2, error_scale = 1)
  )

  expect_named(population_predictors(foo), c("x1", "x2"), ignore.order = TRUE)
})
