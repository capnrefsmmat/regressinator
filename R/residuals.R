#' Augment a model fit with partial residuals for all terms
#'
#' Construct a data frame containing the model data, partial residuals for all
#' predictors, and predictor effects, for use in residual diagnostic plots and
#' other analyses. The result is in tidy form (one row per predictor per
#' observation), allowing it to be easily manipulated for plots and simulations.
#'
#' # Predictors and regressors
#'
#' To define partial residuals, we must distinguish between the *predictors*,
#' the measured variables we are using to fit our model, and the *regressors*,
#' which are calculated from them. In a simple linear model, the regressors are
#' equal to the predictors. But in a model with polynomials, splines, or other
#' nonlinear terms, the regressors may be functions of the predictors.
#'
#' For example, in a regression with a single predictor \eqn{X}, the regression
#' model \eqn{Y = \beta_0 + \beta_1 X + e} has one regressor, \eqn{X}. But if we
#' choose a polynomial of degree 3, the model is \eqn{Y = \beta_0 + \beta_1 X +
#' \beta_2 X^2 + \beta_3 X^3}, and the regressors are \eqn{\{X, X^2, X^3\}}{{X,
#' X^2, X^3}}.
#'
#' Similarly, if we have predictors \eqn{X_1} and \eqn{X_2} and form a model
#' with main effects and an interaction, the regressors are \eqn{\{X_1, X_2, X_1
#' X_2\}}{{X_1, X_2, X_1 X_2}}.
#'
#' Partial residuals are defined in terms of the predictors, not the regressors,
#' and are intended to allow us to see the shape of the relationship between a
#' particular predictor and the response, and to compare it to how we have
#' chosen to model it with regressors. Partial residuals are not well-defined
#' for predictors that have an interaction with other predictors in the model,
#' as the shape of the modeled relationship varies depending on the other
#' predictors they interact with.
#'
#' # Linear models
#'
#' Consider a linear model where \eqn{\mathbb{E}[Y \mid X = x] = \mu(x)}{E[Y | X
#' = x] = \mu(x)}. The mean function \eqn{\mu(x)} is a linear combination of
#' regressors. Let \eqn{\hat \mu}{muhat} be the fitted model and \eqn{\hat
#' \beta_0}{beta0hat} be its intercept.
#'
#' Choose a predictor \eqn{X_f}, the *focal* predictor, to calculate partial
#' residuals for. Write the mean function as \eqn{\mu(X_f, X_o)}, where
#' \eqn{X_f} is the value of the focal predictor, and \eqn{X_o} represents all
#' other predictors.
#'
#' If \eqn{e_i} is the residual for observation \eqn{i}, the partial residual is
#'
#' \deqn{r_{if} = e_i + (\hat \mu(x_{if}, 0) - \hat \beta_0).}{
#' r_if = e_i + (muhat(x_if, 0) - beta0hat).}
#'
#' # Generalized linear models
#'
#' Consider a generalized linear model where \eqn{g(\mathbb{E}[Y \mid X = x]) =
#' \mu(x)}{g(E[Y | X = x]) = \mu(x)}, where \eqn{g} is a link function. Let
#' \eqn{\hat \mu}{muhat} be the fitted model and \eqn{\hat \beta_0}{beta0hat} be
#' its intercept.
#'
#' Let \eqn{e_i} be the *working residual* for observation \eqn{i}, defined to
#' be
#'
#' \deqn{e_i = (y_i - g^{-1}(x_i)) g'(x_i).}
#'
#' Choose a predictor \eqn{X_f}, the *focal* predictor, to calculate partial
#' residuals for. Write \eqn{\mu} as \eqn{\mu(X_f, X_o)}, where \eqn{X_f} is the
#' value of the focal predictor, and \eqn{X_o} represents all other predictors.
#' Hence \eqn{\mu(X_f, X_o)} gives the model's prediction on the link scale.
#'
#' The partial residual is again
#'
#' \deqn{r_{if} = e_i + (\hat \mu(x_{if}, 0) - \hat \beta_0).}{
#' r_if = e_i + (muhat(x_{if}, 0) - beta0hat).}
#'
#' # Interpretation
#'
#' Because the residuals \eqn{e_i} should have mean zero in a well-specified
#' model, plotting the partial residuals against \eqn{x_f} should produce a
#' shape matching the modeled relationship \eqn{\mu}. If the model is wrong, the
#' partial residuals will appear to deviate from the fitted relationship.
#'
#' Additionally, the function \eqn{\mu(X_f, 0)} can be used to show the
#' relationship between the focal predictor and the response. In a linear model,
#' the function is linear; with polynomial or spline regressors, it is
#' nonlinear. This function is the *predictor effect function*, and the
#' estimated predictor effects \eqn{\hat \mu(X_{fi}, 0)}{muhat(X_fi, 0)} are
#' included in this function's output.
#'
#' @param fit The model to obtain residuals for. This can be a model fit with
#'   `lm()` or `glm()`, or any model whose `residuals()` method supports a `type
#'   = "partial"` argument that returns a matrix or data frame of partial
#'   residuals.
#' @param predictors Predictors to calculate partial residuals for. Defaults to
#'   all predictors. Predictors can be specified using tidyselect syntax; see
#'   `help("language", package = "tidyselect")`.
#' @return Data frame (tibble) containing the model data and residuals in tidy
#'   form. There is one row *per predictor* per observation, with the following
#'   columns:
#'
#' \item{obs}{Row number of this observation in the original model data frame.}
#' \item{predictor_name}{Name of the predictor this row gives the partial
#' residual for.}
#' \item{predictor_value}{Value of the predictor this row gives the partial
#' residual for.}
#' \item{partial_resid}{Partial residual for this predictor for this
#' observation.}
#' \item{predictor_effect}{Predictor effect \eqn{\hat \mu(X_{fi},
#' 0)}{muhat(X_fi, 0)} for this observation.}
#'
#' @seealso [binned_residuals()] for the related binned residuals;
#'   `vignette("linear-regression-diagnostics")` for examples of plotting and
#'   interpreting partial residuals
#' @references R. Dennis Cook. "Exploring Partial Residual Plots",
#'   *Technometrics*, 35:4 (1993), 351-362.
#'   <https://doi.org/10.1080/00401706.1993.10485350>
#'
#' Cook, R. Dennis, and Croos-Dabrera, R.
#' "Partial Residual Plots in Generalized Linear Models." *Journal of the
#' American Statistical Association* 93, no. 442 (1998): 730–39.
#' <https://doi.org/10.2307/2670123>
#'
#' Fox, J., & Weisberg, S. (2018). "Visualizing Fit and Lack of Fit in Complex
#' Regression Models with Predictor Effect Plots and Partial Residuals." *Journal
#' of Statistical Software*, 87(9). <https://doi.org/10.18637/jss.v087.i09>
#' @importFrom stats predict
#' @importFrom insight get_predictors get_intercept
#' @importFrom tidyselect everything eval_select
#' @importFrom rlang enquo
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#' @examples
#' fit <- lm(mpg ~ cyl + disp + hp, data = mtcars)
#' partial_residuals(fit)
#'
#' # select predictors with tidyselect syntax
#' partial_residuals(fit, c(disp, hp))
#'
#' # predictors with multiple regressors
#' fit2 <- lm(mpg ~ poly(disp, 2), data = mtcars)
#' partial_residuals(fit2)
#' @export
partial_residuals <- function(fit, predictors = everything()) {
  # TODO Automatically omit predictors in interactions
  predictors <- enquo(predictors)

  pred_data <- get_predictors(fit)
  selection <- eval_select(predictors, pred_data)

  num_predictors <- length(selection)
  predictors <- pred_data[, selection, drop = FALSE]
  predictor_names <- names(predictors)
  nobs <- nrow(predictors)

  intercept <- get_intercept(fit)
  if (is.na(intercept)) {
    intercept <- 0
  }

  resids <- residuals(fit, type = "working")

  # FIXME
  ## interacting_predictors <- Filter(
  ##   function(p) in_interaction(formula(fit), p),
  ##   predictor_names
  ## )

  ## if (length(interacting_predictors) > 0) {
  ##   cli_warn(c("Partial residuals are not defined for predictors in interactions",
  ##              "*" = "Skipped predictors: {.var {interacting_predictors}}"))
  ## }

  ## if (length(interacting_predictors) == length(predictor_names)) {
  ##   return(tibble(
  ##     predictor_name = character(),
  ##     predictor_value = numeric(),
  ##     predictor_effect = numeric(),
  ##     partial_resid = numeric()
  ##   ))
  ## }

  out <- map_dfr(
    ## setdiff(predictor_names, interacting_predictors),
    predictor_names,
    function(predictor) {
      df <- pred_data
      for (p in names(df)) {
        if (p != predictor) {
          df[, p] <- 0
        }
      }

      effect <- predict(fit, newdata = df) - intercept

      tibble(
        predictor_name = predictor,
        predictor_value = predictors[, predictor],
        predictor_effect = effect,
        partial_resid = effect + resids
      )
    }
  )

  return(out)
}

#' Determine if a predictor is involved in an interaction
#'
#' This presently does not work for any term like `poly(x, 3)` or `I(x^2)`,
#' since the terms object contains those names, not the names of the underlying
#' predictors.
#'
#' @param formula Model formula
#' @param predictor Predictor to check for, as character vector
#' @return `TRUE` if in an interaction, `FALSE` otherwise
#' @keywords internal
in_interaction <- function(formula, predictor) {
  terms <- terms(formula)

  # indicates which predictors are in which model terms
  factors <- attr(terms, "factors")

  # drop columns where our predictor does not appear
  factors <- factors[, factors[predictor, ] > 0, drop = FALSE]

  # term with > 1 predictors is an interaction of some kind
  interactions <- colSums(factors) > 1

  return(any(interactions))
}

#' Obtain binned residuals for a model
#'
#' Construct a data frame by binning the fitted values or predictors of a model
#' into discrete bins of equal width, and calculating the average value of the
#' residuals within each bin.
#'
#' In many generalized linear models, the residual plots (Pearson or deviance)
#' are not useful because the response variable takes on very few possible
#' values, causing strange patterns in the residuals. For instance, in logistic
#' regression, plotting the residuals versus covariates usually produces two
#' curved lines.
#'
#' If we first bin the data, i.e. divide up the observations into `breaks` bins
#' based on their fitted values, we can calculate the average residual within
#' each bin. This can be more informative: if a region has 20 observations and
#' its average residual value is large, this suggests those observations are
#' collectively poorly fit. By default, the binning is with respect to the
#' fitted values.
#'
#' @param fit The model to obtain residuals for. This can be a model fit with
#'   `lm()` or `glm()`, or any model that has a `residuals()` method.
#' @param term If `NULL`, bin the residuals with respect to the fitted values.
#'   If non-`NULL`, this should be a string identifying a model term that will
#'   be binned.
#' @param breaks Number of bins to create, or a numeric vector of two or more
#'   unique cut points to use. Passed to `cut()` to create the cuts. If `NULL`,
#'   a default number of breaks is chosen based on the number of rows in the
#'   data.
#' @param ... Additional arguments passed on to `residuals()`. The most useful
#'   additional argument is typically `type`, to select the type of residuals to
#'   produce (such as standardized residuals or deviance residuals).
#' @return Data frame (tibble) with one row per bin, and the following columns:
#'
#' \item{bin}{Bin number.}
#' \item{n}{Number of observations in this bin.}
#' \item{min, max, mean, sd}{Minimum, maximum, mean, and standard deviation of
#' either the fitted values of observations in this bin, or of the term if
#' `term` is not `NULL`.}
#' \item{resid.mean}{Mean residual in this bin.}
#' \item{resid.sd}{Standard deviation of residuals in this bin.}
#'
#' @seealso [partial_residuals()] for the related partial residuals;
#'   `vignette("logistic-regression-diagnostics")` for examples of use and
#'   interpretation of binned residuals in logistic regression;
#'   [bin_by_interval()] and [bin_by_quantile()] to bin data and calculate other
#'   values in each bin
#' @references Gelman, A., Hill, J., and Vehtari, A. (2021). Regression and
#'   Other Stories. Section 14.5. Cambridge University Press.
#' @importFrom dplyr n summarize
#' @importFrom rlang sym .data
#' @importFrom stats residuals sd
#' @importFrom tibble as_tibble
#' @examples
#' fit <- lm(mpg ~ cyl + disp + hp, data = mtcars)
#'
#' binned_residuals(fit, breaks = 5)
#'
#' binned_residuals(fit, term = "cyl")
#' @export
binned_residuals <- function(fit, term = NULL, breaks = NULL, ...) {
  d <- model.frame(fit)

  if (is.null(term)) {
    # fit on the response scale
    d$.fitted <- fitted(fit)
    colname <- sym(".fitted")
  } else {
    if (!(term %in% colnames(d))) {
      cli_abort("term {.var {term}} is not in the model frame provided")
    }

    colname <- sym(term)
  }

  d$.resid <- residuals(fit, ...)

  out <- d |>
    bin_by_interval(!!colname, breaks = breaks) |>
    summarize(
      n = n(),
      min = min(!!colname),
      max = max(!!colname),
      mean = mean(!!colname),
      sd = sd(!!colname),
      resid.mean = mean(.data$.resid),
      resid.sd = sd(.data$.resid)
    )

  return(as_tibble(out))
}

# adapted from https://stackoverflow.com/a/13217607, available CC-BY-SA
#' @importFrom stats terms
response_var <- function(formula) {
    tt <- terms(formula)
    vars <- as.character(attr(tt, "variables"))[-1] # [1] is the list call
    response <- attr(tt, "response") # index of response var

    return(vars[response])
}

#' Augment a model fit with residuals, in "long" format
#'
#' Use `broom::augment()` to augment a model fit with residual and fit
#' information, then reformat the resulting data frame into a "long" format with
#' one row per predictor per observation, to facilitate plotting of the result.
#'
#' The name comes by analogy to `tidyr::pivot_longer()`, and the concept of long
#' versus wide data formats.
#'
#' @param x A model fit object, such as those returned by `lm()` or `glm()`. See
#'   the broom documentation for the full list of model types supported.
#' @param ... Additional arguments passed to `broom::augment()`.
#' @return A data frame (tibble) in similar form to those produced by
#'   `broom::augment()`, but expanded to have one row per predictor per
#'   observation. Columns `.predictor_name` and `.predictor_value` identify the
#'   predictor and its value. An additional column `.obs` records the original
#'   observation numbers so results can be matched to observations in the
#'   original model data.
#' @importFrom broom augment
#' @importFrom tidyr pivot_longer starts_with any_of
#' @examples
#' fit <- lm(mpg ~ cyl + disp + hp, data = mtcars)
#'
#' # each observation appears 3 times, once per predictor:
#' augment_longer(fit)
#' @export
augment_longer <- function(x, ...) {
  out <- augment(x, ...)
  out$.obs <- rownames(out)
  response <- response_var(x)

  return(pivot_longer(out, cols = !starts_with(".") & !any_of(response),
                      names_to = ".predictor_name",
                      values_to = ".predictor_value"))
}

#' Group a data frame into bins
#'
#' Groups a data frame (similarly to `dplyr::group_by()`) based on the values of
#' a column, either by dividing up the range into equal pieces or by quantiles.
#'
#' `bin_by_interval()` breaks the numerical range of that column into
#' equal-sized intervals, or into intervals specified by `breaks`.
#' `bin_by_quantile()` splits the range into pieces based on quantiles of the
#' data, so each interval contains roughly an equal number of observations.
#'
#' @param .data Data frame to bin
#' @param col Column to bin by
#' @param breaks Number of bins to create. `bin_by_interval()` also accepts a
#'   numeric vector of two or more unique cut points to use. If `NULL`, a
#'   default number of breaks is chosen based on the number of rows in the data.
#' @return Grouped data frame, similar to those returned by `dplyr::group_by()`.
#'   An additional column `.bin` indicates the bin number for each group. Use
#'   `dplyr::summarize()` to calculate values within each group, or other dplyr
#'   operations that work on groups.
#' @export
#' @examples
#' suppressMessages(library(dplyr))
#' cars |>
#'   bin_by_interval(speed, breaks = 5) |>
#'   summarize(mean_speed = mean(speed),
#'             mean_dist = mean(dist))
#'
#' cars |>
#'   bin_by_quantile(speed, breaks = 5) |>
#'   summarize(mean_speed = mean(speed),
#'             mean_dist = mean(dist))
#' @importFrom dplyr mutate group_by
#' @importFrom rlang .data
bin_by_interval <- function(.data, col, breaks = NULL) {
  n <- nrow(.data)

  if (is.null(breaks)) {
    breaks <- size_heuristic(n)
  }

  mutate(.data,
         .bin = cut({{ col }}, breaks = breaks, labels = FALSE)) |>
    group_by(.data$.bin)
}

#' @rdname bin_by_interval
#' @importFrom ggplot2 cut_number
#' @export
bin_by_quantile <- function(.data, col, breaks = NULL) {
  n <- nrow(.data)

  if (is.null(breaks)) {
    breaks <- size_heuristic(n)
  }

  mutate(.data,
         .bin = cut_number({{ col }}, n = breaks, labels = FALSE)) |>
    group_by(.data$.bin)
}

size_heuristic <- function(n) {
  # size heuristic taken from the arm package's binned residuals function
  if (n <= 10) {
    floor(n / 2)
  } else if (n < 100) {
    10
  } else {
    floor(sqrt(n))
  }
}

#' Empirically estimate response values on the link scale
#'
#' Calculates the average value of the response variable, and places this on the
#' link scale. Plotting these against a predictor (by dividing the dataset into
#' bins) can help assess the choice of link function.
#'
#' @param response Vector of response variable values.
#' @param family Family object representing the response distribution and link
#'   function. Only the link function will be used.
#' @param na.rm Should `NA` values of the response be stripped? Passed to
#'   `mean()` when calculating the mean of the response.
#' @return Mean response value, on the link scale.
#' @export
#' @examples
#' suppressMessages(library(dplyr))
#' suppressMessages(library(ggplot2))
#'
#' mtcars |>
#'   bin_by_interval(disp, breaks = 5) |>
#'   summarize(
#'     mean_disp = mean(disp),
#'     link = empirical_link(am, binomial())
#'   ) |>
#'   ggplot(aes(x = mean_disp, y = link)) +
#'   geom_point()
empirical_link <- function(response, family, na.rm = FALSE) {
  family <- normalize_family(family)

  ybar <- mean(response, na.rm = na.rm)

  return(family$linkfun(ybar))
}
