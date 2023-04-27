test_that("in_interaction handles all formula cases", {
  skip("in_interaction implementation is incomplete")

  expect_true(in_interaction(y ~ x1 * x2, "x2"))
  expect_false(in_interaction(y ~ x1 + x2, "x2"))
  expect_true(in_interaction(y ~ x1 + I(x1 * x2), "x2"))

  expect_false(in_interaction(y ~ poly(x1, 2) + x2, "x1"))
  expect_false(in_interaction(y ~ poly(x1, 2) + x2, "x2"))

  expect_true(in_interaction(y ~ x1 * x2 - 1, "x1"))

  expect_true(in_interaction(y ~ x1 + x2 + x1:x2, "x1"))

  expect_true(in_interaction(y ~ (x1 + x2 + x3)^3, "x1"))
  expect_false(in_interaction(y ~ (x1 + x2 + x3)^3 + x4, "x4"))
})

test_that("prototype_for produces sensible data frames", {
  expect_equal(
    prototype_for(
      data.frame(x1 = 1:2,
                 x2 = 3:4),
      "x1"
    ),
    data.frame(x1 = 1:2,
               x2 = 0)
  )

  foo <- as.factor(c("a", "b"))
  expect_equal(
    prototype_for(
      data.frame(x1 = 1:2,
                 x2 = foo,
                 x3 = c(TRUE, FALSE),
                 x4 = c("b", "a")),
      "x1"
    ),
    data.frame(x1 = 1:2,
               x2 = foo[1],
               x3 = FALSE,
               x4 = "a")
  )
})

test_that("drop_factors handles all kinds of factors", {
  foo <- data.frame(
    foo = 1:4,
    bar = c(TRUE, FALSE, TRUE, FALSE),
    baz = c("Ducks", "Geese", "Penguins", "Walruses"),
    spam = 2:5
  )

  expected <- data.frame(foo = 1:4, spam = 2:5)

  expect_equal(suppressMessages(drop_factors(foo)), expected)

  # ensure result isn't a vector when there's only one remaining column
  expect_equal(suppressMessages(drop_factors(foo[, -4])),
               expected[, -2, drop = FALSE])
})

test_that("detect_transmutation() rejects factor() calls", {
  expect_no_error(detect_transmutation(foo ~ bar))

  expect_error(detect_transmutation(foo ~ factor(bar)),
               class = "regressinator_transmutation_factor")

  expect_error(detect_transmutation(foo ~ factor(bar) + baz),
               class = "regressinator_transmutation_factor")

  expect_error(detect_transmutation(foo ~ factor(bar) * baz),
               class = "regressinator_transmutation_factor")

  expect_error(detect_transmutation(foo ~ ham + bar + factor(bar):baz),
               class = "regressinator_transmutation_factor")

  expect_no_error(detect_transmutation(foo ~ ham + bar + bar:baz))
})
