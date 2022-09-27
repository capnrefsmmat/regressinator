test_that("model_lineup rejects invalid fn", {
  fn <- function(fit) "duck"

  expect_error(model_lineup(lm(dist ~ speed, data = cars), fn),
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

test_that("parametric_boot_distribution produces new response data", {
  # there was a bug that it would put the simulated data in a column named "y",
  # even if the response variable is not named "y"
  fit <- lm(dist ~ speed, data = cars)

  boot_dist <- parametric_boot_distribution(fit, fn = tidy)

  coefs <- boot_dist |>
    group_by(term) |>
    summarize(sd = sd(estimate))

  expect_true(all(coefs$sd > 0))
})
