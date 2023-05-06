test_that("binomial samples have adjustable size", {
  pop <- population(x = predictor("rnorm"),
                    y = response(0, family = binomial()))
  samples <- pop |> sample_x(n = 100) |> sample_y()

  expect_equal(range(samples$y), c(0, 1))

  pop <- population(x = predictor("rnorm"),
                    y = response(0, family = binomial(), size = 2))
  samples <- pop |> sample_x(n = 100) |> sample_y()

  expect_equal(range(samples$y), c(0, 2))
})

test_that("binomial family rejects invalid size", {
  pop <- population(x = predictor("rnorm"),
                    y = response(0, family = binomial(), size = c(2, 3)))
  samples <- pop |> sample_x(n = 100)

  expect_error(sample_y(samples))

  pop <- population(x = predictor("rnorm"),
                    y = response(0, family = binomial(), size = 2.5))
  samples <- pop |> sample_x(n = 100)

  expect_error(sample_y(samples))
})

test_that("sample_y() throws classed errors", {
  pop <- population(
    x = predictor("rnorm"),
    y = response(foo, error_scale = 1)
  )

  expect_error(pop |> sample_x(10) |> sample_y(),
               class = "regressinator_eval_response")

  pop <- population(
    x = predictor("rnorm"),
    y = response(x, error_scale = x2 + 2)
  )

  expect_error(pop |> sample_x(10) |> sample_y(),
               class = "regressinator_eval_error_scale")

  pop <- population(
    x = predictor("rnorm"),
    y = response(x, family = binomial(), size = foo)
  )

  expect_error(pop |> sample_x(10) |> sample_y(),
               class = "regressinator_eval_size")
})
