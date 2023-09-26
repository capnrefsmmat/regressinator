test_that("model_lineup rejects invalid fn", {
  fn <- function(fit) "duck"

  expect_error(model_lineup(lm(dist ~ speed, data = cars), fn),
               class = "regressinator_diagnostic_class")
})

test_that("sampling_distribution rejects invalid fn", {
  fn <- function(fit) "duck"
  data <- sample_y(
    sample_x(
      population(x = predictor("rnorm"),
                 y = response(2, error_scale = 1)),
      n = 10)
  )

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

test_that("sampling_distribution resamples x when fixed_x = FALSE", {
  # there was a bug that prevented fixed_x = FALSE from working correctly,
  # because sampling_distribution() passed sample_x() the sample, not its parent
  # population

  pop <- population(
    x = predictor("rnorm", mean = 0, sd = 1),
    y = response(4 + 2 * x, family = gaussian(), error_scale = 1.0)
  )

  samp <- pop |>
    sample_x(n = 100) |>
    sample_y()

  fit <- lm(y ~ x, data = samp)

  samples <- sampling_distribution(fit, samp, fn = function(fit) {
    data.frame(x = model.frame(fit)$x)
  }, nsim = 2, fixed_x = FALSE)

  # most of the Xs should not be unique. Technically all, but of course I don't
  # trust floats to be unique like real numbers
  expect_gt(length(unique(samples$x)), 150)
})
