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

