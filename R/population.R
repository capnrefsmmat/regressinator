
# TODO Allow creating a population from a (large) data frame, and sampling from
# it just by sampling with replacement

# TODO Make an example of using the framework for hypothesis testing.

# TODO Support mixture distributions of covariates and of errors, for easy study
# of outliers. Maybe "rmixture" that take multiple distributions as arguments?

#' Specify the distribution of a predictor variable
#'
#' Predictor variables can have any marginal distribution as long as a function
#' is provided to sample from the distribution. Multivariate distributions are
#' also supported: if the random generation function returns multiple columns,
#' multiple random variables will be created, successively numbered.
#'
#' The random generation function must take an argument named `n` specifying the
#' number of draws. For univariate distributions, it should return a vector of
#' length `n`; for multivariate distributions, it should return an array or
#' matrix with `n` rows and a column per variable.
#'
#' Multivariate predictors are successively numbered. For instance, if predictor
#' `X` is specified with
#'
#' ```
#' library(mvtnorm)
#' predictor(dist = "rmvnorm", mean = c(0, 1),
#'          sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2))
#' ```
#'
#' then the population predictors will be named `X1` and `X2`, and will have
#' covariance 0.5.
#'
#' @param dist Name (as character vector) of the function to generate draws from
#'   this predictor's distribution.
#' @param ... Additional arguments to pass to `dist` when generating draws.
#' @return A predictor distribution object
#' @examples
#' # Univariate normal distribution
#' predictor(dist = "rnorm", mean = 10, sd = 2.5)
#'
#' # Multivariate normal distribution
#' library(mvtnorm)
#' predictor(dist = "rmvnorm", mean = c(0, 1, 7))
#' @export
predictor <- function(dist, ...) {
  dist_args <- list(...)

  return(structure(
    list(dist = dist, args = dist_args),
    class = "predictor_dist"))
}

#' @export
print.predictor_dist <- function(x, ...) {
  cat(x$dist, "(", deparse(x$args), ")\n", sep = "")

  invisible(x)
}

#' Specify a response variable in terms of predictors
#'
#' Response variables are related to predictors (and other response variables)
#' through a link function and response distribution. First the expression
#' provided is evaluated using the predictors, to give this response variable's
#' value on the link scale; then the inverse link function and response
#' distribution are used to get the response value. See Details for more
#' information.
#'
#' Response variables are drawn based on a typical generalized linear model
#' setup. Let \eqn{Y} represent the response variable and \eqn{X} represent the
#' predictor variables. We specify that
#'
#' \deqn{Y \sim \text{SomeDistribution}(g^{-1}(\mu(X)))}{%
#'       Y ~ SomeDistribution(g^{-1}(\mu(X)))}
#'
#' where \eqn{\mu(X)} is the expression `expr`, and both the distribution and
#' link function \eqn{g} are specified by the `family` provided. For instance,
#' if the `family` is `gaussian()`, the distribution is Normal and the link is
#' the identity function; if the `family` is `binomial()`, the distribution is
#' binomial and the link is (by default) the logistic link.
#'
#' @param expr An expression, in terms of other predictor or response variables,
#'   giving this predictor's value on the link scale.
#' @param family The family of this response variable, e.g. `gaussian()` for an
#'   ordinary Gaussian linear relationship.
#' @param error_scale Scale factor for errors. Used only for linear families,
#'   such as `gaussian()` and `ols_with_error()`. Errors drawn while simulating
#'   the response variable will be multiplied by this scale factor. The scale
#'   factor can be a scalar value (such as a fixed standard deviation), or an
#'   expression in terms of the predictors, which will be evaluated when
#'   simulating response data. For generalized linear models, leave as `NULL`.
#' @importFrom stats gaussian
#' @export
response <- function(expr, family = gaussian(), error_scale = NULL) {
  response_expr <- substitute(expr)
  error_scale <- substitute(error_scale)

  if (!inherits(family, "family")) {
    family_class <- paste0(class(family), collapse = ", ")
    cli_abort(c("family provided must be a family object",
                "*" = "family provided has class {family_class}"))
  }

  if (!(family$family %in% c("gaussian", "ols_with_error")) &&
        !is.null(error_scale)) {
    cli_warn("error_scale was provided to population(), but family is not gaussian() or ols_with_error(), so it will be ignored")
  }

  if (family$family %in% c("gaussian", "ols_with_error") &&
        is.null(error_scale)) {
    cli_abort("error_scale must be provided for gaussian and ols_with_error families")
  }

  return(structure(
    list(response_expr = response_expr,
         family = family,
         error_scale = error_scale),
    class = "response_dist"))
}

#' @export
print.response_dist <- function(x, ...) {
  cat(x$family$family, "(", deparse(x$response_expr), sep = "")

  if (!is.null(x$error_scale)) {
    cat(", error_scale = ", deparse(x$error_scale), sep = "")
  }
  cat(")\n")

  invisible(x)
}

#' Define the population generalized regression relationship
#'
#' Specifies a hypothetical infinite population of cases. Each case has some
#' predictor variables and one or more response variables. The relationship
#' between the variables and response variables are defined, as well as the
#' population marginal distribution of each predictor variable.
#'
#' @param ... A sequence of named arguments defining predictor and response
#'   variables. These are evaluated in order, so later response variables may
#'   refer to earlier predictor and response variables. All predictors should be
#'   provided first, before any response variables.
#' @return A population object.
#' @seealso [predictor()] and [response()] to define the population;
#'   `sample_x()` and `sample_y()` to draw samples from it
#' @examples
#' # A population with a simple linear relationship
#' population(
#'   x1 = predictor("rnorm", mean = 4, sd = 10),
#'   x2 = predictor("runif", min = 0, max = 10),
#'   y = response(0.7 + 2.2 * x1 - 0.2 * x2, error_scale = 1.0)
#' )
#'
#' # A binary outcome Y, using a binomial family with logistic link
#' population(
#'   x1 = predictor("rnorm", mean = 4, sd = 10),
#'   x2 = predictor("runif", min = 0, max = 10),
#'   y = response(0.7 + 2.2 * x1 - 0.2 * x2,
#'                family = binomial(link = "logit"))
#' )
#'
#' # A population with a simple linear relationship and collinearity. Because X
#' # is bivariate, there will be two predictors, named x1 and x2.
#' library(mvtnorm)
#' population(
#'   x = predictor("rmvnorm", mean = c(0, 1),
#'                 sigma = matrix(c(1, 0.8, 0.8, 1), nrow = 2)),
#'   y = response(0.7 + 2.2 * x1 - 0.2 * x2, error_scale = 1.0)
#' )
#' @importFrom cli cli_abort cli_warn
#' @export
population <- function(...) {
  variables <- list(...)

  return(structure(variables, class = "population"))
}

#' @export
print.population <- function(x, ...) {
  cat("Population with variables:\n")

  for (name in names(x)) {
    cat(name, ": ", sep = "")
    print(x[[name]])
  }

  invisible(x)
}

#' Family representing a linear relationship with non-Gaussian errors
#'
#' This family can represent any non-Gaussian error, provided random variates
#' can be drawn by an R function. A family specified this way can be used to
#' specify a population (via `population()`), but can't be used to estimate a
#' model (such as with `glm()`).
#'
#' @param error Function that can draw random variables from the non-Gaussian
#'   distribution. For example, `rt` draws *t*-distributed random variates. The
#'   function must take an argument `n` indicating how many random variates to
#'   draw (as all random generation functions built into R do).
#' @param ... Further arguments passed to the `error` function to draw random
#'   variates, such as to specify degrees of freedom, shape parameters, or other
#'   parameters of the distribution. These arguments are evaluated with the
#'   model data in the environment, so they can be expressions referring to
#'   model data, such as values of the predictors.
#' @return A family object representing this family.
#' @seealso [custom_family()]
#' @examples
#' # t-distributed errors with 3 degrees of freedom
#' ols_with_error(rt, df = 3)
#'
#' # Cauchy-distributed errors
#' ols_with_error(rcauchy, scale = 3)
#' @importFrom rlang eval_tidy
#' @importFrom stats model.frame fitted gaussian
#' @export
ols_with_error <- function(error, ...) {
  fam <- gaussian()

  error_args <- substitute(alist(...))

  fam$family <- "ols_with_error"

  fam$initialize <- expression(
    cli_abort(c("ols_with_error() cannot be used to fit models, only to specify populations",
                "i" = "to fit models with custom error distribution assumptions,",
                "i" = "derive your own maximum likelihood estimator"))
  )

  fam$simulate <- function(object, nsim, env = model.frame(object), ftd = NULL) {
    if (is.null(ftd)) {
      ftd <- fitted(object)
    }

    n <- length(ftd) * nsim

    # Evaluate the error arguments in the model frame, so they can depend on the
    # covariates TODO What if the error arguments have the wrong length, because
    # the expression provided by the user is wrong? Will they be broadcast
    # correctly if nsim > 1?
    args <- eval_tidy(error_args, data = env)
    args$n <- n
    ftd + do.call(error, args)
  }

  return(fam)
}

#' Family representing a GLM with custom distribution and link function
#'
#' Allows specification of the random component and link function for a response
#' variable. In principle this could be used to specify any GLM family, but it
#' is usually easier to use the predefined families, such as `gaussian()` and
#' `binomial()`.
#'
#' A GLM is specified by a combination of:
#'
#' - Random component, i.e. the distribution that Y is drawn from
#' - Link function relating the mean of the random component to the linear predictor
#' - Linear predictor
#'
#' Using `custom_family()` we can specify the random component and link
#' function, while the linear predictor is set in `population()` when setting up
#' the population relationships. A family specified this way can be used to
#' specify a population (via `population()`), but can't be used to estimate a
#' model (such as with `glm()`).
#'
#' @param distribution The distribution of the random component. This should be
#'   in the form of a function taking one argument, the vector of values on the
#'   inverse link scale, and returning a vector of draws from the distribution.
#' @param inverse_link The inverse link function.
#' @return A family object representing this family
#' @seealso [ols_with_error()]
#' @examples
#' # A zero-inflated Poisson family
#' rzeroinfpois <- function(ys) {
#'   n <- length(ys)
#'   rpois(n, lambda = ys * rbinom(n, 1, prob = 0.4))
#' }
#'
#' custom_family(rzeroinfpois, exp)
#' @importFrom stats binomial
#' @export
custom_family <- function(distribution, inverse_link) {
  # sacrificial family
  fam <- binomial()

  fam$family <- "custom_family"

  fam$initialize <- expression(
    cli_abort(c("custom_family() cannot be used to fit models, only to specify populations",
                "i" = "see ?family for a list of families supported for model fits"))
  )
  fam$link <- paste0("inverse of ", deparse(substitute(inverse_link)))
  fam$simulate <- function(object, nsim, env = model.frame(object), ftd = NULL) {
    if (is.null(ftd)) {
      ftd <- fitted(object)
    }

    return(distribution(ftd))
  }

  fam$linkinv <- inverse_link

  return(fam)
}
