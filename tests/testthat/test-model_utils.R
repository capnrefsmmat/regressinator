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
