#' Draw a data frame of X values from the specified population.
#'
#' @param population Population, as defined by `population()`
#' @param n Number of observations to draw from the population
#' @return Data frame of `n` rows, with columns matching the variables specified
#'   in the population
#' @importFrom assertthat assert_that
#' @export
sample_x <- function(population, n) {
  assert_that(inherits(population, "population"))

  sampled_data <- lapply(
    population$predictors,
    function(var) {
      args <- within(var, rm("dist"))
      args$n <- n
      do.call(var$dist, args)
    })

  # sampled_data is now a named list. names are variable names, entries are
  # vectors (for univariate predictors) or matrices (for multivariate
  # predictors). to get this into a data frame, we must split out the matrices
  for (predictor in names(sampled_data)) {
    if (is.matrix(sampled_data[[predictor]])) {
      pred_data <- sampled_data[[predictor]]
      sampled_data[[predictor]] <- NULL

      for (col in seq_len(ncol(pred_data))) {
        col_name <- paste0(predictor, col)

        sampled_data[[col_name]] <- pred_data[, col]
      }
    }
  }

  return(structure(
    as.data.frame(sampled_data),
    population = population,
    class = c("population_sample", "data.frame")))
}

parent_population <- function(sample) {
  attr(sample, "population")
}

#' Augment the X values with sampled response values
#'
#' Given the data frame of `xs`, use the true relationship in the population to
#' draw random Y values corresponding to those X values.
#'
#' @param xs Sample X values drawn from the population, as obtained from
#'   `sample_x()`
#' @return Data frame of `xs` with additional column named `y`
#' @importFrom cli cli_abort
#' @importFrom stats rbinom rpois rnorm
#' @importFrom assertthat assert_that
#' @examples
#' # A population with a simple linear relationship
#' pop <- population(
#'   0.7 + 2.2 * x1 - 0.2 * x2,
#'   predictors = list(
#'       x1 = list(dist = "rnorm", mean = 4, sd = 10),
#'       x2 = list(dist = "runif", min = 0, max = 10)
#'   ),
#'   error_scale = 1.0)
#'
#' pop |>
#'   sample_x(10) |>
#'   sample_y()
#' @export
sample_y <- function(xs) {
  assert_that(inherits(xs, "population_sample"))

  n <- nrow(xs)
  population <- parent_population(xs)

  # on the response scale
  y_resp <- population$family$linkinv(eval(population$response, envir = xs))

  family_name <- population$family$family

  if (family_name == "gaussian") {
    y_resp <- rnorm(n, mean = y_resp, sd = 1.0) *
      eval(population$error_scale, envir = xs)
  } else if (family_name == "ols_with_error") {
    y_resp <- y_resp +
      population$family$simulate(NULL, 1, env = xs, ftd = rep(0, n)) *
      eval(population$error_scale, envir = xs)
  } else if (family_name == "binomial") {
    y_resp <- rbinom(n, size = 1, prob = y_resp)
  } else if (family_name == "poisson") {
    y_resp <- rpois(n, lambda = y_resp)
  } else if (family_name == "custom_family") {
    y_resp <- population$family$simulate(NULL, 1, env = xs, ftd = y_resp)
  } else {
    cli_abort(c("Unable to simulate from population family",
                "*" = "Population family is {family_name}",
                "i" = "Supported families are gaussian, ols_with_error, binomial, and poisson"))
  }

  xs$y <- y_resp

  return(xs)
}

#' @export
print.population_sample <- function(x, ...) {
  cat("\nSample of ", nrow(x), " observations from", sep = "")

  print(parent_population(x))

  cat("\n")

  NextMethod("print")
}
