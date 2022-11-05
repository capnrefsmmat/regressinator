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
