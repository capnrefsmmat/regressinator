suppressMessages(library(dplyr))

test_that("partial_residuals() produces correct amount of data", {
  fit <- lm(mpg ~ cyl + disp + hp, data = mtcars)

  out <- partial_residuals(fit)

  expect_equal(nrow(out), nrow(mtcars) * 3)
  expect_setequal(unique(out$.predictor_name),
                  c("cyl", "disp", "hp"))

  # tidyselect syntax
  out <- partial_residuals(fit, c(disp, hp))

  expect_equal(nrow(out), nrow(mtcars) * 2)
  expect_setequal(unique(out$.predictor_name),
                  c("disp", "hp"))
})

test_that("partial_residuals() works on models fit to population samples", {
  # partial_residuals() had a bug that caused it to rely on drop=TRUE behavior
  # in data frames; if given a tibble or population sample, it would
  # accidentally produce list columns instead of the correct columns
  pop <- population(
    x = predictor("rnorm"),
    y = response(x, error_scale = 1)
  )

  samp <- pop |>
    sample_x(n = 100) |>
    sample_y()

  fit <- lm(y ~ x, data = samp)

  out <- partial_residuals(fit)

  expect_type(out$.predictor_value, "double") # not a tibble
  expect_setequal(names(out),
                  c("x", ".predictor_name", ".predictor_value",
                    ".predictor_effect", ".partial_resid"))
})

test_that("partial_residuals() omits factors", {
  mtcars$cylinders <- factor(mtcars$cyl)

  fit <- lm(mpg ~ cylinders * disp + hp, data = mtcars)

  out <- partial_residuals(fit)

  expect_equal(nrow(out), nrow(mtcars) * 2)
  expect_setequal(unique(out$.predictor_name),
                  c("disp", "hp"))
})

test_that("partial_residuals() rejects factor() in formulas", {
  fit <- lm(mpg ~ factor(cyl) * disp + hp, data = mtcars)

  expect_error(partial_residuals(fit),
               class = "regressinator_transmutation_factor")
})

test_that("partial_residuals() gives correct results for GLMs", {
  # Can compare to residuals(x, type = "partial") when predictors enter directly
  # as regressors with no transformations; there may be a constant offset,
  # however, so center each

  fit <- glm(cyl ~ drat + wt, family = poisson, data = mtcars)
  pr <- residuals(fit, type = "partial")

  out <- partial_residuals(fit)

  drat_pr <- pr[, "drat"]
  drat_out <- out |> filter(.predictor_name == "drat") |> pull(.partial_resid)

  expect_equal(unname(drat_pr - mean(drat_pr)),
               drat_out - mean(drat_out))

  wt_pr <- pr[, "wt"]
  wt_out <- out |> filter(.predictor_name == "wt") |> pull(.partial_resid)

  expect_equal(unname(wt_pr - mean(wt_pr)),
               wt_out - mean(wt_out))
})

test_that("partial_residuals() handles offsets", {
  fit <- lm(mpg ~ hp, offset = qsec, data = mtcars)

  out <- partial_residuals(fit) |>
    pull(.partial_resid)

  pr <- residuals(fit, type = "partial")
  pr <- pr[, "hp"] - mean(pr[, "hp"])

  # as in previous test, centered values should match
  expect_equal(out - mean(out), unname(pr))
})

test_that("binned_residuals() produces correct amount of data", {
  fit <- lm(mpg ~ hp + qsec, data = mtcars)

  out <- binned_residuals(fit, breaks = 5)

  expect_equal(nrow(out), 5 * 2)
  expect_setequal(unique(out$predictor_name), c("hp", "qsec"))
})

test_that("binned_residuals() omits factors", {
  mtcars$cylinders <- factor(mtcars$cyl)

  fit <- lm(mpg ~ cylinders * disp + hp, data = mtcars)

  out <- binned_residuals(fit, breaks = 5)

  expect_equal(nrow(out), 5 * 2)
  expect_setequal(unique(out$predictor_name),
                  c("disp", "hp"))
})

test_that("binned_residuals() rejects factor() in formulas", {
  fit <- lm(mpg ~ factor(cyl) * disp + hp, data = mtcars)

  expect_error(binned_residuals(fit, breaks = 5),
               class = "regressinator_transmutation_factor")
})

test_that("augment_longer() produces correct amount of data", {
  fit <- lm(mpg ~ cyl + disp + hp, data = mtcars)

  out <- augment_longer(fit)

  expect_equal(nrow(out), nrow(mtcars) * 3)
  expect_setequal(unique(out$.predictor_name),
                  c("cyl", "disp", "hp"))
})

test_that("augment_longer() omits factors", {
  mtcars$cylinders <- factor(mtcars$cyl)

  fit <- lm(mpg ~ cylinders * disp + hp, data = mtcars)

  out <- augment_longer(fit)

  expect_equal(nrow(out), nrow(mtcars) * 2)
  expect_setequal(unique(out$.predictor_name),
                  c("disp", "hp"))
})

test_that("augment_longer() keeps factors if there are no numerics", {
  mtcars$cylinders <- factor(mtcars$cyl)
  mtcars$engine_shape <- factor(mtcars$vs)

  fit <- lm(mpg ~ cylinders + engine_shape, data = mtcars)

  out <- augment_longer(fit)

  expect_equal(nrow(out), nrow(mtcars) * 2)
})
