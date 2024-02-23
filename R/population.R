
# TODO Allow creating a population from a (large) data frame, and sampling from
# it just by sampling with replacement

#' Specify the distribution of a predictor variable
#'
#' Predictor variables can have any marginal distribution as long as a function
#' is provided to sample from the distribution. Multivariate distributions are
#' also supported: if the random generation function returns multiple columns,
#' multiple random variables will be created. If the columns are named, the
#' random variables will be named accordingly; otherwise, they will be
#' successively numbered.
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
#' predictor(dist = rmvnorm, mean = c(0, 1),
#'           sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2))
#' ```
#'
#' then the population predictors will be named `X1` and `X2`, and will have
#' covariance 0.5.
#'
#' If the multivariate predictor has named columns, the names will be used
#' instead. For instance, if predictor `X` generates a matrix with columns `A`
#' and `B`, the population predictors will be named `XA` and `XB`.
#'
#' @param dist Function to generate draws from this predictor's distribution,
#'   provided as a function or as a string naming the function.
#' @param ... Additional arguments to pass to `dist` when generating draws.
#' @return A `predictor_dist` object, to be used in `population()` to specify a
#'   population distribution
#' @examples
#' # Univariate normal distribution
#' predictor(dist = rnorm, mean = 10, sd = 2.5)
#'
#' # Multivariate normal distribution
#' library(mvtnorm)
#' predictor(dist = rmvnorm, mean = c(0, 1, 7))
#'
#' # Multivariate with named columns
#' rmulti <- function(n) {
#'   cbind(treatment = rbinom(n, size = 1, prob = 0.5),
#'         confounder = rnorm(n)
#'   )
#' }
#' predictor(dist = rmulti)
#' @export
predictor <- function(dist, ...) {
  dist_name <- deparse(substitute(dist))
  dist_args <- list(...)

  return(structure(
    list(dist = dist, dist_name = dist_name, args = dist_args),
    class = "predictor_dist"))
}

#' @export
print.predictor_dist <- function(x, ...) {
  if (length(x$args) > 0) {
    cat(x$dist_name, "(", deparse(x$args), ")\n", sep = "")
  } else {
    cat(x$dist_name, "()\n", sep = "")
  }

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
#' \deqn{Y \mid X \sim \text{SomeDistribution},}{%
#'       Y | X ~ SomeDistribution,}
#'
#' where
#'
#' \deqn{\mathbb{E}[Y \mid X = x] = g^{-1}(\mu(x)).}{%
#'       E[Y | X = x] = g^{-1}(\mu(x)).}
#'
#' Here \eqn{\mu(X)} is the expression `expr`, and both the distribution and
#' link function \eqn{g} are specified by the `family` provided. For instance,
#' if the `family` is `gaussian()`, the distribution is Normal and the link is
#' the identity function; if the `family` is `binomial()`, the distribution is
#' binomial and the link is (by default) the logistic link.
#'
#' ## Response families
#'
#' The following response families are supported.
#'
#' \describe{
#' \item{`gaussian()`}{
#' The default family is `gaussian()` with the identity link function,
#' specifying the relationship
#'
#' \deqn{Y \mid X \sim \text{Normal}(\mu(X), \sigma^2),}{%
#'       Y | X ~ Normal(mu(X), \sigma^2),}
#'
#' where \eqn{\sigma^2} is given by `error_scale`.
#' }
#'
#' \item{`ols_with_error()`}{Allows specification of custom non-Normal error
#' distributions, specifying the relationship
#'
#' \deqn{Y = \mu(X) + e,}
#'
#' where \eqn{e} is drawn from an arbitrary distribution, specified by the
#' `error` argument to `ols_with_error()`.
#' }
#'
#' \item{`binomial()`}{Binomial responses include binary responses (as in logistic
#' regression) and responses giving a total number of successes out of a number
#' of trials. The response has distribution
#'
#' \deqn{Y \mid X \sim \text{Binomial}(N, g^{-1}(\mu(X))),}{%
#'       Y | X ~ Binomial(N, g^{-1}(\mu(X))),
#' }
#'
#' where \eqn{N} is set by the `size` argument and \eqn{g} is the link function.
#' The default link is the logistic link, and others can be chosen with the
#' `link` argument to `binomial()`. The default \eqn{N} is 1, representing a
#' binary outcome.
#' }
#'
#' \item{`poisson()`}{Poisson-distributed responses with distribution
#'
#' \deqn{Y \mid X \sim \text{Poisson}(g^{-1}(\mu(X))),}{%
#'       Y | X ~ Poisson(g^{-1}(\mu(X))),
#' }
#'
#' where \eqn{g} is the link function. The default link is the log link, and
#' others can be chosen with the `link` argument to `poisson()`.
#' }
#'
#' \item{`custom_family()`}{Responses drawn from an arbitrary distribution with
#' arbitrary link function, i.e.
#'
#'  \deqn{Y \mid X \sim \text{SomeDistribution}(g^{-1}(\mu(X))),}{%
#'        Y | X ~ SomeDistribution(g^{-1}(\mu(X))),}
#'
#' where both \eqn{g} and SomeDistribution are specified by arguments to
#' `custom_family()`.
#' }
#' }
#'
#' ## Evaluation and scoping
#'
#' The `expr`, `error_scale`, and `size` arguments are evaluated only when
#' simulating data for this response variable. They are evaluated in an
#' environment with access to the predictor variables and the preceding response
#' variables, which they can refer to by name. Additionally, these arguments can
#' refer to variables in scope when the enclosing `population()` was defined.
#' See the Examples below.
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
#' @param size When the `family` is `binomial()`, this is the number of trials
#'   for each observation. Defaults to 1, as in logistic regression. May be
#'   specified either as a vector of the same length as the number of
#'   observations or as a scalar. May be written terms of other predictor or
#'   response variables. For other families, `size` is ignored.
#' @return A `response_dist` object, to be used in `population()` to specify a
#'   population distribution
#' @importFrom stats gaussian
#' @importFrom cli cli_abort cli_warn
#' @seealso [predictor()] and [population()] to define populations;
#'   [ols_with_error()] and [custom_family()] for custom response distributions
#' @examples
#' # Defining a binomial response. The expressions can refer to other predictors
#' # and to the environment where the `population()` is defined:
#' slope1 <- 2.5
#' slope2 <- -3
#' intercept <- -4.6
#' size <- 10
#' population(
#'   x1 = predictor(rnorm),
#'   x2 = predictor(rnorm),
#'   y = response(intercept + slope1 * x1 + slope2 * x2,
#'                family = binomial(), size = size)
#' )
#' @importFrom rlang enquo quo_is_null
#' @export
response <- function(expr, family = gaussian(), error_scale = NULL,
                     size = 1L) {
  response_expr <- enquo(expr)
  error_scale <- enquo(error_scale)
  size <- enquo(size)

  family <- normalize_family(family)

  if (!(family$family %in% c("gaussian", "ols_with_error")) &&
        !quo_is_null(error_scale)) {
    cli_warn("{.arg error_scale} was provided to {.fn population}, but family is not {.fn gaussian} or {.fn ols_with_error}, so it will be ignored")
  }

  if (family$family %in% c("gaussian", "ols_with_error") &&
        quo_is_null(error_scale)) {
    cli_abort("{.arg error_scale} must be provided for {.fn gaussian} and {.fn ols_with_error} families",
              class = "regressinator_error_scale")
  }

  return(structure(
    list(response_expr = response_expr,
         family = family,
         error_scale = error_scale,
         size = size),
    class = "response_dist"))
}

#' @export
#' @importFrom rlang quo_get_expr
print.response_dist <- function(x, ...) {
  cat(x$family$family, "(", deparse(quo_get_expr(x$response_expr)), sep = "")

  if (!quo_is_null(x$error_scale)) {
    cat(", error_scale = ", deparse(quo_get_expr(x$error_scale)), sep = "")
  }
  if (x$family$family == "binomial") {
    cat(", size = ", deparse(quo_get_expr(x$size)), sep = "")
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
#' linear_pop <- population(
#'   x1 = predictor(rnorm, mean = 4, sd = 10),
#'   x2 = predictor(runif, min = 0, max = 10),
#'   y = response(0.7 + 2.2 * x1 - 0.2 * x2, error_scale = 1.0)
#' )
#'
#' # A population whose response depends on local variables
#' slope <- 2.2
#' intercept <- 0.7
#' sigma <- 2.5
#' variable_pop <- population(
#'   x = predictor(rnorm),
#'   y = response(intercept + slope * x, error_scale = sigma)
#' )
#'
#' # Response error scale is heteroskedastic and depends on predictors
#' heteroskedastic_pop <- population(
#'   x1 = predictor(rnorm, mean = 4, sd = 10),
#'   x2 = predictor(runif, min = 0, max = 10),
#'   y = response(0.7 + 2.2 * x1 - 0.2 * x2,
#'                error_scale = 1 + x2 / 10)
#' )
#'
#' # A binary outcome Y, using a binomial family with logistic link
#' binary_pop <- population(
#'   x1 = predictor(rnorm, mean = 4, sd = 10),
#'   x2 = predictor(runif, min = 0, max = 10),
#'   y = response(0.7 + 2.2 * x1 - 0.2 * x2,
#'                family = binomial(link = "logit"))
#' )
#'
#' # A binomial outcome Y, with 10 trials per observation, using a logistic link
#' # to determine the probability of success for each trial
#' binomial_pop <- population(
#'   x1 = predictor(rnorm, mean = 4, sd = 10),
#'   x2 = predictor(runif, min = 0, max = 10),
#'   y = response(0.7 + 2.2 * x1 - 0.2 * x2,
#'                family = binomial(link = "logit"),
#'                size = 10)
#' )
#'
#' # Another binomial outcome, but the number of trials depends on another
#' # predictor
#' binom_size_pop <- population(
#'   x1 = predictor(rnorm, mean = 4, sd = 10),
#'   x2 = predictor(runif, min = 0, max = 10),
#'   trials = predictor(rpois, lambda = 20),
#'   y = response(0.7 + 2.2 * x1 - 0.2 * x2,
#'                family = binomial(link = "logit"),
#'                size = trials)
#' )
#'
#' # A population with a simple linear relationship and collinearity. Because X
#' # is bivariate, there will be two predictors, named x1 and x2.
#' library(mvtnorm)
#' collinear_pop <- population(
#'   x = predictor(rmvnorm, mean = c(0, 1),
#'                 sigma = matrix(c(1, 0.8, 0.8, 1), nrow = 2)),
#'   y = response(0.7 + 2.2 * x1 - 0.2 * x2, error_scale = 1.0)
#' )
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

#' Get the predictors of a population
#'
#' @param population Population object
#' @return Named list of predictors
#' @keywords internal
population_predictors <- function(population) {
  Filter(
    function(v) { inherits(v, "predictor_dist") },
    population
  )
}

#' Get the response variables of a population
#'
#' @param population Population object
#' @return Named list of response variables
#' @keywords internal
population_response <- function(population) {
  Filter(
    function(v) { inherits(v, "response_dist") },
    population
  )
}

#' Family representing a linear relationship with non-Gaussian errors
#'
#' The `ols_with_error()` family can represent any non-Gaussian error, provided
#' random variates can be drawn by an R function. A family specified this way
#' can be used to specify a population (via `population()`), but can't be used
#' to estimate a model (such as with `glm()`).
#'
#' @param error Function that can draw random variables from the non-Gaussian
#'   distribution, or a string giving the name of the function. For example,
#'   `rt` draws *t*-distributed random variates. The function must take an
#'   argument `n` indicating how many random variates to draw (as all random
#'   generation functions built into R do).
#' @param ... Further arguments passed to the `error` function to draw random
#'   variates, such as to specify degrees of freedom, shape parameters, or other
#'   parameters of the distribution. These arguments are evaluated with the
#'   model data in the environment, so they can be expressions referring to
#'   model data, such as values of the predictors.
#' @return A family object representing this family.
#' @seealso [custom_family()] for fully custom families, including for GLMs
#' @examples
#' # t-distributed errors with 3 degrees of freedom
#' ols_with_error(rt, df = 3)
#'
#' # A linear regression with t-distributed error, using error_scale to make
#' # errors large
#' population(
#'   x1 = predictor(rnorm, mean = 4, sd = 10),
#'   x2 = predictor(runif, min = 0, max = 10),
#'   y = response(0.7 + 2.2 * x1 - 0.2 * x2,
#'                family = ols_with_error(rt, df = 4),
#'                error_scale = 2.5)
#' )
#'
#' # Cauchy-distributed errors
#' ols_with_error(rcauchy, scale = 3)
#'
#' # A contaminated error distribution, where
#' # 95% of observations are Gaussian and 5% are Cauchy
#' rcontaminated <- function(n) {
#'   contaminant <- rbinom(n, 1, prob = 0.05)
#'
#'   return(ifelse(contaminant == 1,
#'                 rcauchy(n, scale = 20),
#'                 rnorm(n, sd = 1)))
#' }
#' ols_with_error(rcontaminated)
#' @importFrom rlang eval_tidy
#' @importFrom stats model.frame fitted gaussian
#' @export
ols_with_error <- function(error, ...) {
  fam <- gaussian()

  error_args <- substitute(alist(...))

  fam$family <- "ols_with_error"

  fam$initialize <- expression(
    cli_abort(c("{.fn ols_with_error} cannot be used to fit models, only to specify populations",
                "i" = "to fit models with custom error distribution assumptions, derive your own maximum likelihood estimator"))
  )

  fam$simulate <- function(object, nsim, env = model.frame(object), ftd = NULL) {
    if (is.null(ftd)) {
      ftd <- fitted(object)
    }

    if (nsim > 1) {
      cli_abort(c("{.fn ols_with_error} family simulation does not support {.arg nsim} > 1",
                  "*" = "This error should not happen; please report it as a bug"),
                class = "regressinator_error_nsim")
    }

    n <- length(ftd) * nsim

    # Evaluate the error arguments in the model frame, so they can depend on the
    # covariates
    args <- eval_tidy(error_args, data = env)
    args$n <- n
    err <- do.call(error, args)
    err_len <- length(err)

    if (err_len != n) {
      cli_abort(c("error function provided to {.fn ols_with_error} returned incorrect output length",
                  "*" = "data has {.val {n}} observations, but function only returned {.val {err_len}} value{?s}"),
                class = "regressinator_error_length")
    }

    return(ftd + err)
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
#' @seealso [ols_with_error()] for the special case of linear regression with
#'   custom error distribution
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
    cli_abort(c("{.fn custom_family} cannot be used to fit models, only to specify populations",
                "i" = "see {.topic stats::family} for a list of families supported for model fits"))
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
