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

test_that("parametric_boot_distribution handles transformed predictors", {
  # Using I(x^2) in a model formula means model.matrix() returns an `I(x^2)`
  # column instead of an `x` column; need to fetch the original predictors to
  # successfully use `update()`.

  fit <- lm(mpg ~ I(hp^2), data = mtcars)

  out <- model_lineup(fit)

  expect_equal(nrow(out), 20 * nrow(mtcars))

  # Ensure matrix columns, such as those provided from splines, also work
  fit <- lm(mpg ~ splines::ns(hp, df = 3), data = mtcars)

  out <- model_lineup(fit)

  expect_equal(nrow(out), 20 * nrow(mtcars))
})

test_that("parametric_boot_distribution includes all columns from the data", {
  linear_pop <- population(
    x1 = predictor("rnorm", mean = 4, sd = 10),
    x2 = predictor("runif", min = 0, max = 10),
    y = response(
      0.7 + 2.2 * x1 - 0.2 * x2,
      family = gaussian(),
      error_scale = 1.5
    )
  )

  sample <- linear_pop |>
    sample_x(n = 50) |>
    sample_y()

  # notice null_fit does not contain x2
  fit <- lm(y ~ x1 + x2, data = sample)
  null_fit <- lm(y ~ x1, data = sample)

  # we want the simulated datasets to contain x2 so fit can be updated
  null_dist <- parametric_boot_distribution(null_fit, fit, data = sample,
                                            nsim = 100)

  expect_equal(nrow(null_dist), 3 * 100)
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
