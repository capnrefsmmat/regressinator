test_that("diagnose_model rejects invalid fn", {
  fn <- function(fit) "duck"

  expect_error(diagnose_model(lm(dist ~ speed, data = cars), fn),
               class = "regressinator_diagnostic_class")
})

test_that("sampling_distribution rejects invalid fn", {
  fn <- function(fit) "duck"
  data <- sample_x(
    population(x = predictor("rnorm"),
               y = response(2, error_scale = 1)),
    n = 10)

  expect_error(sampling_distribution(lm(y ~ x, data = data), data, fn),
               class = "regressinator_diagnostic_class")
})
