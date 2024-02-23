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

test_that("sample_x() throws classed error", {
  # invalid distribution function
  pop <- population(x = predictor("nonexistent_fn", mean = 0),
                    y = response(x, error_scale = 1))

  expect_error(sample_x(pop, 10),
               class = "regressinator_sample_dist")
})

test_that("sample_x() names multivariate predictors", {
  # unnamed multivariate
  runnamed <- function(n) {
    cbind(1:n, 1:n)
  }
  pop <- population(x = predictor(runnamed))

  expect_named(sample_x(pop, 10), c("x1", "x2"))

  # named multivariate
  rnamed <- function(n) {
    cbind(a = 1:n, b = 1:n)
  }
  pop <- population(x = predictor(rnamed))

  expect_named(sample_x(pop, 10), c("xa", "xb"))
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

test_that("response expr evaluated in right environment", {
  # Refer to local variables
  slope <- 2.5
  intercept <- 1.0
  pop <- population(
    x = predictor("rnorm"),
    y = response(intercept + slope * x, error_scale = 1)
  )

  expect_no_error(pop |> sample_x(10) |> sample_y())

  # refer to other environments
  foo <- function() {
    slope1 <- 2.5
    intercept1 <- 1.0
    return(population(
      x = predictor("rnorm"),
      y = response(intercept1 + slope1 * x, error_scale = 1)
    ))
  }

  expect_no_error(foo() |> sample_x(10) |> sample_y())
})

test_that("error_scale argument evaluated in right environment", {
  # Refer to prior predictors
  pop <- population(
    x = predictor("rnorm"),
    y = response(x, error_scale = x**2)
  )

  expect_no_error(pop |> sample_x(10) |> sample_y())

  # Refer to local variables
  sigma <- 2
  pop <- population(
    x = predictor("rnorm"),
    y = response(x, error_scale = sigma)
  )

  expect_no_error(pop |> sample_x(10) |> sample_y())
})

test_that("binomial response size argument evaluated in right environment", {
  # Should be able to refer to other predictors
  pop <- population(
    x = predictor("rpois", lambda = 10),
    y = response(x / 10, family = binomial(), size = x)
  )

  expect_no_error(pop |>
                    sample_x(10) |>
                    sample_y())

  # Should also be able to refer to local variables
  s <- 14
  pop <- population(
    x = predictor("rnorm"),
    y = response(x, family = binomial(), size = s)
  )

  expect_no_error(out <- pop |>
                    sample_x(10) |>
                    sample_y())
  expect_true(all(out$y <= s))

  # And variables in different environments
  foo <- function() {
    q <- 20
    return(population(
      x = predictor("rnorm"),
      y = response(x, family = binomial(), size = q)
    ))
  }

  expect_no_error(out <- foo() |>
                    sample_x(10) |>
                    sample_y())
  expect_true(all(out$y <= 20))
})
