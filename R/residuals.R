#' Augment a model fit with partial residuals for all terms
#'
#' Construct a data frame containing the model data, partial residuals, fitted
#' values, and linear predictor, for use in residual diagnostic plots and other
#' analyses. The result is in tidy form (one row per predictor per observation),
#' allowing it to be easily manipulated for plots and simulations.
#'
#' Consider a generalized linear model setting. We define \eqn{\mathbb{E}[Y \mid
#' X] = \mu(X)}{E[Y | X] = \mu(X)}. Let \eqn{g} be the link function relating
#' \eqn{\mu(X)} to the linear predictors, so \eqn{g(\mu(X)) = \beta^T X}. We
#' define the partial residual for covariate \eqn{j} to be
#'
#' \deqn{(Y - \hat \mu) g'(\hat \mu) + \hat \beta^T X_j}
#'
#' where \eqn{g'(\hat \mu)} is the first derivative of the link function with
#' respect to \eqn{\mu}, and \eqn{\hat \mu} gives the predictions from our
#' fitted model.
#'
#' For example, in a linear model, \eqn{g'(\hat \mu) = 1}, so the partial
#' residuals for covariate \eqn{j} are the ordinary residuals plus the
#' contribution of \eqn{X_j} to the model fit. Since \eqn{\hat \mu} is also a
#' linear combination of the predictors, this gives us the partial residual of
#' \eqn{Y} on every predictor *except* \eqn{X_j}:
#'
#' \deqn{Y - \sum_{k \neq j} \hat \beta_k X_k}
#'
#' When plotted against \eqn{X_j}, we expect this to produce a line with slope
#' \eqn{\hat \beta_j}. If instead the scatterplot has a curve, this suggests the
#' relationship between \eqn{X_j} and \eqn{Y} is nonlinear.
#'
#' @param fit The model to obtain residuals for. This can be a model fit with
#'   `lm()` or `glm()`, or any model whose `residuals()` method supports a `type
#'   = "partial"` argument that returns a matrix or data frame of partial
#'   residuals.
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
#' \item{fitted}{Fitted value (on the response scale) for this observation,
#' obtained using `fitted()`.}
#' \item{linear_predictor}{Linear predictor for this observation, obtained using
#' `predict()`. For linear models, the fitted values and linear predictor are
#' identical.}
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
#' @importFrom stats predict fitted
#' @importFrom tibble as_tibble
#' @examples
#' fit <- lm(mpg ~ cyl + disp + hp, data = mtcars)
#' partial_residuals(fit)
#' @export
partial_residuals <- function(fit) {
  partials <- as.data.frame(residuals(fit, type = "partial"))

  d <- model.frame(fit)
  predictor_names <- names(partials)
  num_predictors <- ncol(partials)
  nobs <- nrow(partials)

  yhat <- fitted(fit)
  lp <- predict(fit) # linear predictor, since glms default to type="link"

  out <- data.frame(obs = rep(seq_len(num_predictors), times = nobs))

  for (obs in seq_len(nobs)) {
    for (pred in seq_len(num_predictors)) {
      row <- (obs - 1) * num_predictors + pred

      out$fitted[row] <- yhat[obs]
      out$linear_predictor[row] <- lp[obs]

      out$predictor_name[row] <- predictor_names[pred]
      out$predictor_value[row] <- d[obs, predictor_names[pred]]

      out$partial_resid[row] <- partials[obs, predictor_names[pred]]
    }
  }

  return(as_tibble(out))
}

#' Obtain binned residuals for a model
#'
#' Construct a data frame by binning the fitted values or predictors of a model
#' into discrete bins, and calculating the average value of the residuals within
#' each bin.
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
#'   [bin_by_interval()] to bin data and calculate other values in each bin
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
#' a column, by dividing the numerical range of that column into equal-sized
#' intervals and collecting all rows in each interval.
#'
#' @param .data Data frame to bin
#' @param col Column to bin by
#' @param breaks Number of bins to create, or a numeric vector of two or more
#'   unique cut points to use. Passed to `cut()` to create the cuts. If `NULL`,
#'   a default number of breaks is chosen based on the number of rows in the
#'   data.
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
#' @importFrom ggplot2 cut_interval
#' @importFrom dplyr mutate group_by
#' @importFrom rlang .data
bin_by_interval <- function(.data, col, breaks = NULL) {
  n <- nrow(.data)

  if (is.null(breaks)) {
    # size heuristic taken from the arm package's binned residuals function
    breaks <- if (n <= 10) {
      floor(n / 2)
    } else if (n < 100) {
      10
    } else {
      floor(sqrt(n))
    }
  }

  mutate(.data,
         .bin = cut({{ col }}, breaks = breaks, labels = FALSE)) |>
    group_by(.data$.bin)
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
