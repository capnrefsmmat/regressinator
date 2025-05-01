#' Augment a model fit with partial residuals for all terms
#'
#' Construct a data frame containing the model data, partial residuals for all
#' quantitative predictors, and predictor effects, for use in residual
#' diagnostic plots and other analyses. The result is in tidy form (one row per
#' predictor per observation), allowing it to be easily manipulated for plots
#' and simulations.
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
#' chosen to model it with regressors. Partial residuals are not useful for
#' categorical (factor) predictors, and so these are omitted.
#'
#' Besides regressors, a model may have offset terms, which enter the model with
#' a fixed coefficient of 1. These are fixed to their mean value for partial
#' residual calculations.
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
#' Setting \eqn{X_o = 0} means setting all other numeric predictors to 0; factor
#' predictors are set to their first (baseline) level.
#'
#' # Generalized linear models
#'
#' Consider a generalized linear model where \eqn{g(\mathbb{E}[Y \mid X = x]) =
#' \eta(x)}{g(E[Y | X = x]) = \eta(x)}, where \eqn{g} is a link function and
#' \eqn{\eta(x)} is the linear predictor. Let \eqn{\hat \eta}{etahat} be the
#' fitted linear predictor and \eqn{\hat \beta_0}{beta0hat} be its intercept.
#' Let \eqn{\hat \mu(x) = g^{-1}(\hat \eta(x))}{muhat(x) = g^{-1}(etahat(x))} be
#' the fitted mean function.
#'
#' The *working residual* \eqn{e_i} for observation \eqn{i} is
#'
#' \deqn{e_i = (y_i - \hat \mu(x_i)) g'(\hat \mu(x_i)).}{e_i = (y_i -
#' muhat(x_i)) g'(muhat(x_i)).}
#'
#' Choose a predictor \eqn{X_f}, the *focal* predictor, to calculate partial
#' residuals for. Write \eqn{\eta} as \eqn{\eta(X_f, X_o)}, where \eqn{X_f} is the
#' value of the focal predictor, and \eqn{X_o} represents all other predictors.
#' Hence \eqn{\eta(X_f, X_o)} gives the model's prediction on the link scale.
#'
#' The partial residual is again
#'
#' \deqn{r_{if} = e_i + (\hat \eta(x_{if}, 0) - \hat \beta_0).}{
#' r_if = e_i + (etahat(x_{if}, 0) - beta0hat).}
#'
#' # Interpretation
#'
#' In linear regression, because the residuals \eqn{e_i} should have mean zero
#' in a well-specified model, plotting the partial residuals against \eqn{x_f}
#' should produce a shape matching the modeled relationship \eqn{\mu}. If the
#' model is wrong, the partial residuals will appear to deviate from the fitted
#' relationship. Provided the regressors are uncorrelated or approximately
#' linearly related to each other, the plotted trend should approximate the true
#' relationship between \eqn{x_f} and the response.
#'
#' In generalized linear models, this is approximately true if the link function
#' \eqn{g} is approximately linear over the range of observed \eqn{x} values.
#'
#' Additionally, the function \eqn{\mu(X_f, 0)} (in linear models) or
#' \eqn{\eta(X_f, 0)} (in generalized linear models) can be used to show the
#' relationship between the focal predictor and the response. In a linear model,
#' the function is linear; with polynomial or spline regressors, it is
#' nonlinear. This function is the *predictor effect function*, and the
#' estimated predictor effects are included in this function's output.
#'
#' # Limitations
#'
#' Factor predictors (as factors, logical, or character vectors) are detected
#' automatically and omitted. However, if a numeric variable is converted to
#' factor in the model formula, such as with `y ~ factor(x)`, the function
#' cannot determine the appropriate type and will raise an error. Create factors
#' as needed in the source data frame *before* fitting the model to avoid this
#' issue.
#'
#' @param fit The model to obtain residuals for. This can be a model fit with
#'   `lm()` or `glm()`, or any model with a `predict()` method that accepts a
#'   `newdata` argument.
#' @param predictors Predictors to calculate partial residuals for. Defaults to
#'   all predictors, skipping factors. Predictors can be specified using
#'   tidyselect syntax; see `help("language", package = "tidyselect")` and the
#'   examples below.
#' @return Data frame (tibble) containing the model data and residuals in tidy
#'   form. There is one row *per selected predictor* per observation. All
#'   predictors are included as columns, plus the following additional columns:
#'
#' \item{.obs}{Row number of this observation in the original model data frame.}
#' \item{.predictor_name}{Name of the predictor this row gives the partial
#' residual for.}
#' \item{.predictor_value}{Value of the predictor this row gives the partial
#' residual for.}
#' \item{.partial_resid}{Partial residual for this predictor for this
#' observation.}
#' \item{.predictor_effect}{Predictor effect \eqn{\hat \mu(X_{if},
#' 0)}{muhat(X_if, 0)} for this observation.}
#'
#' @seealso [binned_residuals()] for the related binned residuals;
#'   [augment_longer()] for a similarly formatted data frame of ordinary
#'   residuals; `vignette("linear-regression-diagnostics")`,
#'   `vignette("logistic-regression-diagnostics")`, and
#'   `vignette("other-glm-diagnostics")` for examples of plotting and
#'   interpreting partial residuals
#' @references R. Dennis Cook (1993). "Exploring Partial Residual Plots",
#'   *Technometrics*, 35:4, 351-362. \doi{10.1080/00401706.1993.10485350}
#'
#' Cook, R. Dennis, and Croos-Dabrera, R. (1998).
#' "Partial Residual Plots in Generalized Linear Models." *Journal of the
#' American Statistical Association* 93, no. 442: 730–39. \doi{10.2307/2670123}
#'
#' Fox, J., & Weisberg, S. (2018).
#' "Visualizing Fit and Lack of Fit in Complex Regression Models with Predictor
#' Effect Plots and Partial Residuals." *Journal of Statistical Software*,
#' 87(9). \doi{10.18637/jss.v087.i09}
#' @importFrom stats predict formula
#' @importFrom insight find_offset get_predictors get_intercept
#' @importFrom tidyselect everything eval_select
#' @importFrom rlang enquo
#' @importFrom purrr map_dfr
#' @importFrom tibble as_tibble
#' @examples
#' fit <- lm(mpg ~ cyl + disp + hp, data = mtcars)
#' partial_residuals(fit)
#'
#' # You can select predictors with tidyselect syntax:
#' partial_residuals(fit, c(disp, hp))
#'
#' # Predictors with multiple regressors are supported:
#' fit2 <- lm(mpg ~ poly(disp, 2), data = mtcars)
#' partial_residuals(fit2)
#'
#' # Allowing an interaction by number of cylinders is fine, but partial
#' # residuals are not generated for the factor. Notice the factor must be
#' # created first, not in the model formula:
#' mtcars$cylinders <- factor(mtcars$cyl)
#' fit3 <- lm(mpg ~ cylinders * disp + hp, data = mtcars)
#' partial_residuals(fit3)
#' @export
partial_residuals <- function(fit, predictors = everything()) {
  # Detect and reject factor() in formulas
  detect_transmutation(formula(fit))

  predictors <- enquo(predictors)

  pred_data <- get_predictors(fit)
  selection <- eval_select(predictors, pred_data,
                           allow_empty = FALSE, allow_rename = FALSE)

  predictors <- drop_factors(pred_data[, selection, drop = FALSE])
  predictor_names <- names(predictors)

  # If there is an offset, we must provide this in the prototype. Fix it at its
  # mean value
  offset <- find_offset(fit)
  if (!is.null(offset)) {
    offset_value <- mean(get_data(fit)[[offset]])
  }

  intercept <- get_intercept(fit)
  if (is.na(intercept)) {
    intercept <- 0
  }

  resids <- residuals(fit, type = "working")

  out <- map_dfr(
    predictor_names,
    function(predictor) {
      df <- prototype_for(pred_data, predictor)

      if (!is.null(offset)) {
        df[[offset]] <- offset_value
      }

      effect <- predict(fit, newdata = df) - intercept

      out <- pred_data
      out$.predictor_name <- predictor
      out$.predictor_value <- predictors[[predictor]]
      out$.predictor_effect <- effect
      out$.partial_resid <- effect + resids
      return(out)
    }
  )

  return(as_tibble(out))
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
#' collectively poorly fit. We can also bin each predictor and calculate
#' averages within those bins, allowing the detection of misspecification for
#' specific model terms.
#'
#' # Limitations
#'
#' Factor predictors (as factors, logical, or character vectors) are detected
#' automatically and omitted. However, if a numeric variable is converted to
#' factor in the model formula, such as with `y ~ factor(x)`, the function
#' cannot determine the appropriate type and will raise an error. Create factors
#' as needed in the source data frame *before* fitting the model to avoid this
#' issue.
#'
#' @param fit The model to obtain residuals for. This can be a model fit with
#'   `lm()` or `glm()`, or any model that has `residuals()` and `fitted()`
#'   methods.
#' @param predictors Predictors to calculate binned residuals for. Defaults to
#'   all predictors, skipping factors. Predictors can be specified using
#'   tidyselect syntax; see `help("language", package = "tidyselect")` and the
#'   examples below. Specify `predictors = .fitted` to obtain binned residuals
#'   versus fitted values.
#' @param breaks Number of bins to create. If `NULL`, a default number of breaks
#'   is chosen based on the number of rows in the data.
#' @param ... Additional arguments passed on to `residuals()`. The most useful
#'   additional argument is typically `type`, to select the type of residuals to
#'   produce (such as standardized residuals or deviance residuals).
#' @return Data frame (tibble) with one row per bin *per selected predictor*,
#'   and the following columns:
#'
#' \item{.bin}{Bin number.}
#' \item{n}{Number of observations in this bin.}
#' \item{predictor_name}{Name of the predictor that has been binned.}
#' \item{predictor_min, predictor_max, predictor_mean, predictor_sd}{Minimum,
#' maximum, mean, and standard deviation of the predictor (or fitted values).}
#' \item{resid_mean}{Mean residual in this bin.}
#' \item{resid_sd}{Standard deviation of residuals in this bin.}
#'
#' @seealso [partial_residuals()] for the related partial residuals;
#'   `vignette("logistic-regression-diagnostics")` and
#'   `vignette("other-glm-diagnostics")` for examples of use and interpretation
#'   of binned residuals in logistic regression and GLMs; [bin_by_interval()]
#'   and [bin_by_quantile()] to bin data and calculate other values in each bin
#' @references Gelman, A., Hill, J., and Vehtari, A. (2021). *Regression and
#'   Other Stories*. Section 14.5. Cambridge University Press.
#' @importFrom dplyr n summarize
#' @importFrom insight get_predictors
#' @importFrom purrr map_dfr
#' @importFrom rlang sym .data
#' @importFrom stats residuals sd formula
#' @importFrom tibble as_tibble
#' @examples
#' fit <- lm(mpg ~ disp + hp, data = mtcars)
#'
#' # Automatically bins both predictors:
#' binned_residuals(fit, breaks = 5)
#'
#' # Just bin one predictor, selected with tidyselect syntax. Multiple could be
#' # selected with c().
#' binned_residuals(fit, disp, breaks = 5)
#'
#' # Bin the fitted values:
#' binned_residuals(fit, predictors = .fitted)
#'
#' # Bins are made using the predictor, not regressors derived from it, so here
#' # disp is binned, not its polynomial
#' fit2 <- lm(mpg ~ poly(disp, 2), data = mtcars)
#' binned_residuals(fit2)
#' @export
binned_residuals <- function(fit, predictors = !".fitted", breaks = NULL,
                             ...) {
  # Detect and reject factor() in formulas
  detect_transmutation(formula(fit))

  predictors <- enquo(predictors)

  pred_data <- get_predictors(fit)
  pred_data$.fitted <- fitted(fit)

  selection <- eval_select(predictors, pred_data,
                           allow_empty = FALSE, allow_rename = FALSE)

  predictors <- drop_factors(pred_data[, selection, drop = FALSE])
  predictor_names <- names(predictors)

  predictors$.resid <- residuals(fit, ...)

  out <- map_dfr(
    predictor_names,
    function(predictor) {
      pred <- sym(predictor)
      predictors |>
        bin_by_quantile(!!pred, breaks = breaks) |>
        summarize(
          n = n(),
          predictor_name = predictor,
          predictor_min = min(!!pred),
          predictor_max = max(!!pred),
          predictor_mean = mean(!!pred),
          predictor_sd = sd(!!pred),
          resid_mean = mean(.data$.resid),
          resid_sd = sd(.data$.resid)
        )
    }
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
#' # Limitations
#'
#' Factor predictors (as factors, logical, or character vectors) can't coexist
#' with numeric variables in the `.predictor_value` column. If there are some
#' numeric and some factor predictors, the factor predictors will automatically
#' be omitted. If all predictors are factors, they will be combined into one
#' factor with all levels. However, if a numeric variable is converted to factor
#' in the model formula, such as with `y ~ factor(x)`, the function cannot
#' determine the appropriate types and will raise an error. Create factors as
#' needed in the source data frame *before* fitting the model to avoid this
#' issue.
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
#' @importFrom dplyr relocate select
#' @importFrom tidyr pivot_longer starts_with any_of
#' @seealso [partial_residuals()], [binned_residuals()]
#' @examples
#' fit <- lm(mpg ~ cyl + disp + hp, data = mtcars)
#'
#' # each observation appears 3 times, once per predictor:
#' augment_longer(fit)
#' @export
augment_longer <- function(x, ...) {
  # Detect and reject factor() in formulas
  detect_transmutation(formula(x))

  out <- augment(x, ...)
  out$.obs <- rownames(out)
  response <- response_var(x)

  predictor_cols <- select(out, !starts_with(".") & !any_of(response))

  factor_cols <- factor_columns(predictor_cols)
  if (all(factor_cols)) {
    factor_names <- NULL
  } else {
    # Omit factors if they're going to be combined with numerics
    factor_names <- names(predictor_cols)[factor_cols]
  }

  pivot_longer(out,
               cols = !starts_with(".") & !any_of(response) &
                 !any_of(factor_names),
               names_to = ".predictor_name",
               values_to = ".predictor_value") |>
    relocate(".predictor_name", ".predictor_value")
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
#'   In `bin_by_quantile()`, if the number of unique values of the column is
#'   smaller than `breaks`, fewer bins will be produced.
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
#' @importFrom rlang .data enquo as_name
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
  col <- enquo(col)

  # drop = TRUE is necessary because in data frames, this indexing produces a
  # vector for which length() is appropriate; but tibble indexing always
  # produces a tibble, where length() is the number of columns, unless we force
  # drop behavior.
  n_unique <- length(unique(.data[, as_name(col), drop = TRUE]))

  if (is.null(breaks)) {
    breaks <- size_heuristic(n)
  }

  breaks <- min(breaks, n_unique - 1)

  mutate(.data,
         .bin = cut_number(!!col, n = breaks, labels = FALSE)) |>
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

#' Augment data with randomized quantile residuals
#'
#' Generates a data frame containing a model's predictors, the residuals, and
#' the randomized quantile residuals as additional columns.
#'
#' Randomized quantile residuals provide more interpretable residuals for
#' generalized linear models (GLMs), such as logistic regression. See Dunn and
#' Smyth (1996) for details, or review the examples provided in
#' `vignette("DHARMa", package="DHARMa")`.
#'
#' Let \eqn{F_Y(y; x, \beta)} be the predicted cumulative distribution function
#' for \eqn{Y} when \eqn{X = x}, using the fitted GLM. When the response is
#' continuous, the randomized quantile residual for observation \eqn{i} is
#'
#' \deqn{r_{q,i} = F_Y(y_i; x_i, \hat \beta).}
#'
#' When the response is discrete, let
#'
#' \deqn{a_i = \lim_{y \uparrow y_i} F_Y(y; x_i, \hat \beta)}
#'
#' and
#'
#' \deqn{b_i = F_Y(y_i; x_i, \hat \beta),}
#'
#' then draw the randomized quantile residual as
#'
#' \deqn{r_{q,i} \sim \text{Uniform}(a_i, b_i).}{r_{q,i} ~ Uniform(a_i, b_i).}
#'
#' As cumulative distributions are left-continuous, this "jitters" the values
#' between the discrete steps, resulting in a residual that is uniformly
#' distributed when the model is correct.
#'
#' Some definitions of randomized quantile residuals transform the resulting
#' values using the standard normal inverse cdf, so they are normally
#' distributed. That step is omitted here, as uniform residuals are easy to work
#' with.
#'
#' # Implementation details
#'
#' Uses `broom::augment()` to generate the data frame, then uses the [DHARMa
#' package](https://cran.r-project.org/package=DHARMa) to generate randomized
#' quantile residuals for the model.
#'
#' @param x Fitted model to obtain randomized quantile residuals from
#' @param ... Additional arguments to pass to `broom::augment()`
#' @return Data frame with one row per observation used to fit `x`, including a
#'   `.quantile.resid` column containing the quantile residuals. See
#'   `broom::augment()` and its methods for details of other columns.
#'
#' For `augment_quantile_longer()`, the output is in "long" format with one row
#' per predictor per observation. Columns `.predictor_name` and
#' `.predictor_value` identify the predictor and its value. An additional column
#' `.obs` records the original observation numbers so results can be matched to
#' observations in the original model data. See Limitations in
#' `augment_longer()` for limitations on factor predictors.
#' @importFrom DHARMa simulateResiduals
#' @importFrom broom augment
#' @references Dunn, Peter K., and Gordon K. Smyth (1996).
#'   "Randomized Quantile Residuals." *Journal of Computational and Graphical
#'   Statistics* 5 (3): 236–44. \doi{10.2307/1390802}
#' @seealso `vignette("logistic-regression-diagnostics")` and
#'   `vignette("other-glm-diagnostics")` for examples of plotting and
#'   interpreting randomized quantile residuals; [augment_longer()];
#'   [broom::augment()]
#' @export
augment_quantile <- function(x, ...) {
  out <- augment(x, ...)

  dh <- simulateResiduals(x)

  out$.quantile.resid <- residuals(dh)

  return(out)
}

#' @export
#' @rdname augment_quantile
augment_quantile_longer <- function(x, ...) {
  # Detect and reject factor() in formulas
  detect_transmutation(formula(x))

  out <- augment_quantile(x, ...)
  out$.obs <- rownames(out)
  response <- response_var(x)

  predictor_cols <- select(out, !starts_with(".") & !any_of(response))

  factor_cols <- factor_columns(predictor_cols)
  if (all(factor_cols)) {
    factor_names <- NULL
  } else {
    # Omit factors if they're going to be combined with numerics
    factor_names <- names(predictor_cols)[factor_cols]
  }

  pivot_longer(out,
               cols = !starts_with(".") & !any_of(response) &
                 !any_of(factor_names),
               names_to = ".predictor_name",
               values_to = ".predictor_value") |>
    relocate(".predictor_name", ".predictor_value")
}
