
# TODO Allow creating a population from a (large) data frame, and sampling from
# it just by sampling with replacement

# TODO Allow multiple named response variables. Only question: is providing
# named arguments to `population()` misleading, because `y = x1 + x2` looks like
# exact quality instead of specifying a model that may include a link function?

# TODO Make an example of using the framework for hypothesis testing.

#' Define the population generalized regression relationship
#'
#' Specifies a hypothetical infinite population of cases. Each case has some
#' predictor variables, as well as a response variable. The relationship between
#' the variables and the response variable is defined, as well as the population
#' marginal distribution of each predictor variable.
#'
#' The relationship between predictors and response is inspired by generalized
#' linear models, but allows nonlinear relationships. Specifically,
#'
#' \deqn{g(E[Y \mid X]) = f(X),}{g(E[Y | X]) = f(X),}
#'
#' where `g` is a link function and `f` an arbitrary function of the predictors.
#' `f` is specified by the expression provided in `response`, while the link
#' function `g` is defined by the `family` chosen.
#'
#' The `predictors` are provided as a named list. The name defines the variable
#' name; the value is a list with entry `dist` specifying the function that
#' draws random variates from the predictor's distribution, and additional
#' entries giving arguments to be passed when drawing random variates (such as
#' mean, standard deviation, or shape parameters of the distribution).
#'
#' Each predictor is drawn independently of the other predictors. To specify
#' dependent or correlated predictors, use a function that generates
#' multivariate draws, such as `mvtnorm::rmnvorm()`. Each column of the
#' generated data will be assigned as a population predictor variable, numbered
#' successively. For example, if predictor `X` is specified with
#'
#' ```
#' predictors = list(
#'   X = list(dist = "rmvnorm", mean = c(0, 1),
#'            sigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2)))
#' ```
#'
#' then the population predictors will be named `X1` and `X2`, and will have
#' covariance 0.5.
#'
#' @param response An expression giving the response variable in terms of the
#'   predictor variables. (When the family is not Gaussian, gives the response
#'   on the link scale.)
#' @param predictors A named list of predictor variables in the population. See
#'   Details for specification of distributions.
#' @param family The family of this GLM, e.g. `gaussian()` for an ordinary
#'   Gaussian linear relationship.
#' @param error_scale Scale factor for errors. Used only for linear families,
#'   such as `gaussian()` and `ols_with_error()`. Errors drawn while simulating
#'   the response variable will be multiplied by this scale factor. The scale
#'   factor can be a scalar value (such as a fixed standard deviation), or an
#'   expression in terms of the predictors, which will be evaluated when
#'   simulating response data. For generalized linear models, leave as `NULL`.
#'
#' @examples
#' # A population with a simple linear relationship
#' population(0.7 + 2.2 * x1 - 0.2 * x2,
#'            predictors = list(
#'              x1 = list(dist = "rnorm", mean = 4, sd = 10),
#'              x2 = list(dist = "runif", min = 0, max = 10)
#'            ),
#'            error_scale = 1.0)
#'
#' # A binary outcome Y, using a binomial family with logistic link
#' population(0.7 + 2.2 * x1 - 0.2 * x2,
#'            predictors = list(
#'              x1 = list(dist = "rnorm", mean = 4, sd = 10),
#'              x2 = list(dist = "runif", min = 0, max = 10)
#'            ),
#'            family = binomial(link = "logit")
#' )
#'
#' # A population with a simple linear relationship and collinearity
#' library(mvtnorm)
#' population(0.7 + 2.2 * x1 - 0.2 * x2,
#'            predictors = list(
#'              x = list(dist = "rmvnorm", mean = c(0, 1),
#'                       sigma = matrix(c(1, 0.8, 0.8, 1), nrow = 2))
#'            ),
#'            error_scale = 1.0)
#' @importFrom cli cli_abort cli_warn
#' @importFrom stats gaussian
#' @export
population <- function(response, predictors, family = gaussian(),
                       error_scale = NULL) {
  response <- substitute(response)
  error_scale <- substitute(error_scale)

  if (!(family$family %in% c("gaussian", "ols_with_error")) &&
        !is.null(error_scale)) {
    cli_warn("error_scale was provided to population(), but family is not gaussian() or ols_with_error(), so it will be ignored")
  }

  if (family$family %in% c("gaussian", "ols_with_error") &&
        is.null(error_scale)) {
    cli_abort("error_scale must be provided for gaussian and ols_with_error families")
  }

  structure(
    list(response = response,
         predictors = predictors,
         family = family,
         error_scale = error_scale),
    class = "population")
}

print_variable <- function(name, desc) {
  cat(name, ": ", sep = "")

  dput(desc)
}

#' @export
print.population <- function(x, ...) {
  cat("\nPopulation with predictor variables:\n")

  for (name in names(x$predictors)) {
    print_variable(name, x$predictors[[name]])
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
  fam$simulate <- function(object, nsim, env = model.frame(object), ftd = NULL) {
    if (is.null(ftd)) {
      ftd <- fitted(object)
    }

    return(distribution(ftd))
  }

  fam$linkinv <- inverse_link

  return(fam)
}

## ## example 1: nonlinear OLS
pop <- population(
  0.7 + 4 * x1 - 2 * x2**2,
  predictors = list(
    x1 = list(dist = "rnorm", sd = 8),
    x2 = list(dist = "rnorm", sd = 4)
  ),
  family = gaussian(),
  error_scale = 1.0
)

## real_data <- sample_x(pop, 30) %>%
##   sample_y(pop)

## fit <- lm(y ~ x1 + x2, data = real_data)

## # 1.1: show residuals vs fitted values
## diagnose_model(fit) %>%
##   ggplot(aes(x = .fitted, y = .resid)) +
##   geom_point() +
##   facet_wrap(~ .sample)

## # 1.2: QQ plots
## diagnose_model(fit, n = 10) %>%
##   ggplot(aes(sample = .resid)) +
##   geom_qq() +
##   facet_wrap(~ .sample)


## # example 2: GLM
## pop <- population(
##   0.7 + 4 * x1 - 2 * x2,
##   predictors = list(
##     x1 = list(dist = "rnorm", sd = 0.8),
##     x2 = list(dist = "rnorm", sd = 0.4)
##   ),
##   family = binomial(link = "logit")
## )

## real_data <- sample_x(pop, 50) %>%
##   sample_y(pop)

## fit <- glm(y ~ x1 + x2, family = binomial, data = real_data)

## diagnose_model(fit)  %>%
##   ggplot(aes(x = .fitted, y = .resid)) +
##   geom_point() +
##   facet_wrap(~ .sample)

## ## example 3: OLS but the true error is non-normal
## pop <- population(
##   0.7 + 4 * x1 - 2 * x2,
##   predictors = list(
##     x1 = list(dist = "rnorm", sd = 8),
##     x2 = list(dist = "rnorm", sd = 4)
##   ),
##   family = ols_with_error(rt, df = 2),
##   error_scale = 1.0
## )

## real_data <- sample_x(pop, 30) %>%
##   sample_y(pop)

## fit <- lm(y ~ x1 + x2, data = real_data)
## cis <- confint(fit)

## coefs <- matrix(NA, nrow = 1000, ncol = 3)
## for (b in seq_len(1000)) {
##   new_sim <- real_data %>%
##     sample_y(pop)

##   new_fit <- lm(y ~ x1 + x2, data = new_sim)
##   coefs[b, ] <- coef(new_fit)
## }

## # Minimal bias
## colMeans(coefs)

## # But poor coverage TODO this isn't coverage!
## mean(coefs[, 2] > cis["x1", "2.5 %"] & coefs[, 2] < cis["x1", "97.5 %"])

## # Let's see QQ plots
## diagnose_model(fit, n = 20) %>%
##   ggplot(aes(sample = .resid)) +
##   geom_qq() +
##   geom_qq_line() +
##   facet_wrap(~ .sample)

## ## Example 4: Bad family
## pop <- population(
##   0.7 + 4 * x1 - 2 * x2,
##   predictors = list(
##     x1 = list(dist = "rnorm", sd = 0.8),
##     x2 = list(dist = "rnorm", sd = 0.4)
##   ),
##   family = quasibinomial(link = "logit")
## )

## real_data <- sample_x(pop, 50) %>%
##   sample_y(pop)

## ## Example 5: multivariate
## p <- population(0.7 + 2.2 * x1 - 0.2 * x2,
##                 predictors = list(
##                   x = list(dist = rmvnorm, mean = c(0, 1),
##                            sigma = matrix(c(1, 0.8, 0.8, 1), nrow = 2)),
##                   z = list(dist = rnorm, sd = 4)
##                 ),
##                 error_scale = 1.0)
