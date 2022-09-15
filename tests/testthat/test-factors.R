test_that("rfactor runs with default probabilities", {
  # a previous bug caused this to fail
  expect_equal(length(rfactor(10, c("duck", "goose"))), 10)
})

test_that("by_level calculates correctly", {
  foo <- factor(c("spam", "ham", "spam", "ducks"))
  res <- c("spam" = 4, "ham" = 10, "spam" = 4, "ducks" = 16.7)

  expect_equal(by_level(foo, spam = 4, ham = 10, ducks = 16.7), res)
  expect_equal(by_level(foo, c(spam = 4, ham = 10, ducks = 16.7)), res)
})

test_that("by_level rejects invalid arguments", {
  expect_error(by_level(1:10, 1, 2, 3), class = "regressinator_by_level_arg")

  expect_warning(by_level(letters, a = 2, b = 4),
                 class = "regressinator_by_level_missing_level")
})
