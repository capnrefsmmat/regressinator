library(mvtnorm)

test_that("design_x.population checks predictor names", {
  test_pop <- population(
    x = predictor("rnorm"),
    z = predictor("rmvnorm", mean = 0:2), # produces z1 through z3
    y = response(x + z1 + z2, error_scale = 1)
  )

  # Designed predictors must be part of the population
  bad_design <- data.frame(foo = 1:10)
  expect_error(design_x(test_pop, bad_design),
               class = "regressinator_design_preds")

  # Can't design the response
  resp_design <- data.frame(y = 1:10)
  expect_error(design_x(test_pop, resp_design),
               class = "regressinator_design_preds")

  # Multivariate predictors: z becomes z1, z2, z3. No way to know without
  # sampling it first, but it should be accepted
  okay_design <- data.frame(z1 = 1:10)
  expect_no_error(design_x(test_pop, okay_design))
})
