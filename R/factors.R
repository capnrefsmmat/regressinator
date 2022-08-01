#' Draw random values from a factor variable
#'
#' To specify the population distribution of a factor variable, specify the
#' probability for each of its factor levels. When drawn from the population,
#' factor levels are drawn with replacement according to their probability.
#'
#' @param n Number of values to draw
#' @param levels Character vector specifying the levels for the factor
#' @param prob Vector specifying the probability for each factor level
#' @return Sample of `n` values from `levels`, drawn in proportion to their
#'   probabilities. By default, levels are equally likely.
#' @examples
#' rfactor(5, c("foo", "bar", "baz"), c(0.4, 0.3, 0.3))
#' @export
rfactor <- function(n, levels, prob = 1 / length(levels)) {
  return(factor(sample(levels, size = n, prob = prob, replace = TRUE),
                levels))
}

#' Convert factor levels to numeric values
#'
#' Replace each entry in a vector with its corresponding numeric value, for
#' instance to use a factor variable to specify intercepts for different groups
#' in a regression model.
#'
#' @param x Vector of factor values
#' @param ... Mapping from factor levels to values. Can be provided either as a
#'   series of named arguments, whose names correspond to factor levels, or as a
#'   single named vector.
#' @return Named vector of same length as `x`, with values replaced with those
#'   specified. Names are the original factor level name.
#' @examples
#' foo <- factor(c("spam", "ham", "spam", "ducks"))
#'
#' by_level(foo, spam = 4, ham = 10, ducks = 16.7)
#'
#' by_level(foo, c("spam" = 4, "ham" = 10, "ducks" = 16.7))
#' @export
by_level <- function(x, ...) {
  levels <- list(...)

  if (is.null(names(levels))) {
    if (length(levels) == 1) {
      # we've been passed a single named vector, not a list of associations.
      # just use that named vector
      levels <- levels[[1]]
    } else {
      cli_abort("by_level() should be provided either a single named vector or a sequence of named arguments")
    }
  } else {
    levels <- unlist(levels)
  }

  observed_levels <- unique(x)
  if (!all(observed_levels %in% names(levels))) {
    undef_levels <- paste0(setdiff(observed_levels, names(levels)),
                           collapse = ", ")

    cli_warn(c("not all factor levels were assigned a value in by_level()",
               "*" = "levels without a value: {undef_levels}",
               "i" = "levels without a value will be given value NA"))
  }

  return(levels[as.character(x)])
}
