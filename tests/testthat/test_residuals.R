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
