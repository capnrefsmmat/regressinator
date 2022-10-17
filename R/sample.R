#' Draw a data frame from the specified population.
#'
#' Sampling is split into two steps, for predictors and for response variables,
#' to allow users to choose which to simulate. `sample_x()` will only sample
#' predictor variables, and `sample_y()` will augment a data frame of predictors
#' with columns for response variables, overwriting any already present. Hence
#' one can use `sample_y()` as part of a simulation with fixed predictors, for
#' instance.
#'
#' @param population Population, as defined by `population()`.
#' @param n Number of observations to draw from the population.
#' @return Data frame (tibble) of `n` rows, with columns matching the variables
#'   specified in the population.
#' @importFrom tibble as_tibble
#' @export
sample_x <- function(population, n) {
  if (!inherits(population, "population")) {
    cli_abort("population argument must be a population defined with `population()`")
  }

  predictors <- Filter(
    function(v) { inherits(v, "predictor_dist") },
    population
  )

  sampled_data <- lapply(
    predictors,
    function(var) {
      args <- var$args
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

  sampled_data <- as_tibble(sampled_data)

  return(structure(
    sampled_data,
    population = population,
    class = c("population_sample", class(sampled_data))))
}

parent_population <- function(sample) {
  attr(sample, "population")
}

#' @param xs Data frame of predictor values drawn from the population, as
#'   obtained from `sample_x()`.
#' @importFrom cli cli_abort
#' @importFrom stats rbinom rpois rnorm
#' @examples
#' # A population with a simple linear relationship
#' pop <- population(
#'   x1 = predictor("rnorm", mean = 4, sd = 10),
#'   x2 = predictor("runif", min = 0, max = 10),
#'   y = response(0.7 + 2.2 * x1 - 0.2 * x2, error_scale = 1.0)
#' )
#'
#' xs <- pop |>
#'   sample_x(5)
#'
#' xs
#'
#' xs |>
#'   sample_y()
#' @export
#' @rdname sample_x
sample_y <- function(xs) {
  if (!inherits(xs, "population_sample")) {
    cli_abort(c("data passed to {.fn sample_y} must be from {.fn sample_x}",
                "i" = "other data frames do not have the necessary population attributes specifying the response distribution"))
  }

  n <- nrow(xs)
  population <- parent_population(xs)

  responses <- Filter(
    function(v) { inherits(v, "response_dist") },
    population
  )

  for (response_name in names(responses)) {
    response <- responses[[response_name]]

    # value on the response scale
    y_resp <- response$family$linkinv(eval(response$response_expr, envir = xs))

    family_name <- response$family$family

    if (family_name == "gaussian") {
      y_resp <- rnorm(n, mean = y_resp,
                      sd = eval(response$error_scale, envir = xs))
    } else if (family_name == "ols_with_error") {
      y_resp <- y_resp +
        response$family$simulate(NULL, 1, env = xs, ftd = rep(0, n)) *
        eval(response$error_scale, envir = xs)
    } else if (family_name == "binomial") {
      size <- eval(response$size, envir = xs)

      if (!isTRUE(all.equal(size, as.integer(size)))) {
        cli_abort("{.arg size} for {.fn binomial} families must be an integer or vector of integers")
      }

      if (!(length(size) == 1 || length(size) == length(y_resp))) {
        cli_abort(c("{.arg size} for {.fn binomial} families must be either length 1 or have one entry per observation",
                    "*" = "{.arg size} has length {.val {length(size)}}, but data has length {.val {length(y_resp)}}"))
      }

      y_resp <- rbinom(n, size = size, prob = y_resp)
    } else if (family_name == "poisson") {
      y_resp <- rpois(n, lambda = y_resp)
    } else if (family_name == "custom_family") {
      y_resp <- response$family$simulate(NULL, 1, env = xs, ftd = y_resp)
    } else {
      cli_abort(c("Unable to simulate from population family {.val {family_name}}",
                  "i" = "Supported families are gaussian, ols_with_error, binomial, custom_family, and poisson"))
    }

    xs[[response_name]] <- y_resp
  }

  return(xs)
}

#' @export
print.population_sample <- function(x, ...) {
  cat("Sample of ", nrow(x), " observations from\n", sep = "")

  print(parent_population(x))

  cat("\n")

  NextMethod("print")
}
