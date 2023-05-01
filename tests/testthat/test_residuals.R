library(dplyr)

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
