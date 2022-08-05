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
#' \eqn{Y} on every predictor *except* \eqn{X_j}.
#'
#' @param fit The model to obtain residuals for. This can be a model fit with
#'   `lm()` or `glm()`, or any model whose `residuals()` method supports a `type
#'   = "partial"` argument that returns a matrix or data frame of partial
#'   residuals.
#' @return Data frame containing the model data and residuals in tidy form.
#'   There is one row *per predictor* per observation, with the following
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
#' @seealso [binned_residuals()]
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

  return(out)
}

#' Obtained binned residuals for a model
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
#' If we first bin the data, i.e. divide up the observations into `n_bins` bins
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
#' @param n_bins The number of discrete bins to use. If `NULL`, the number of
#'   bins is chosen heuristically to scale with the number of observations, to
#'   ensure each bin has a reasonable number of observations.
#' @param ... Additional arguments passed on to `residuals()`. The most useful
#'   additional argument is typically `type`, to select the type of residuals to
#'   produce (such as standardized residuals or deviance residuals).
#' @return Data frame with one row per bin.
#'
#' \item{.bin}{Bin number.}
#' \item{.n}{Number of observations in this bin.}
#' \item{.fitted.min, .fitted.max, .fitted.mean, .fitted.sd}{Minimum, maximum,
#' mean, and standard deviation of the fitted values of observations in this
#' bin. If `term` is not `NULL`, these are instead named after the term, e.g. if
#' `term = "x1"` these are `x1.min` and so on.}
#' \item{.resid.mean}{Mean residual in this bin.}
#' \item{.resid.sd}{Standard deviation of residuals in this bin.}
#'
#' @seealso [partial_residuals()]
#' @references Gelman, A. and Hill, J. (2006). Data Analysis Using Regression
#'   and Multilevel/Hierarchical Models. Cambridge University Press.
#' @importFrom stats residuals sd
#' @export
binned_residuals <- function(fit, term = NULL, n_bins = NULL, ...) {
  if (is.null(term)) {
    # fit on the response scale
    x_vals <- fitted(fit)
    colname <- ".fitted"
  } else {
    d <- model.frame(fit)

    if (!(term %in% colnames(d))) {
      cli_abort("term '{term}' is not in the model frame provided")
    }

    x_vals <- d[, term]
    colname <- term
  }

  resids <- residuals(fit, ...)

  n <- length(x_vals)

  if (is.null(n_bins)) {
    # size heuristic taken from the arm package
    n_bins <- if (n <= 10) {
      floor(n / 2)
    } else if (n < 100) {
      100
    } else {
      floor(sqrt(n))
    }
  }

  binned <- cut(x_vals, n_bins, labels = FALSE)

  uniq_bins <- unique(binned)

  out <- data.frame(.bin = uniq_bins)

  for (row in seq_len(nrow(out))) {
    bin <- out$.bin[row]

    out$.n[row] <- sum(binned == bin)

    bin_x <- x_vals[binned == bin]

    out[row, paste0(colname, c(".min", ".max", ".mean", ".sd"))] <-
      c(min(bin_x), max(bin_x), mean(bin_x), sd(bin_x))

    bin_resids <- resids[binned == bin]
    out$.resid.mean[row] <- mean(bin_resids)
    out$.resid.sd[row] <- sd(bin_resids)
  }

  return(out)
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
#' @return A data frame in similar form to those produced by `broom::augment()`,
#'   but expanded to have one row per predictor per observation. Columns
#'   `.predictor_name` and `.predictor_value` identify the predictor and its
#'   value. An additional column `.obs` records the original observation numbers
#'   so results can be matched to observations in the original model data.
#' @importFrom broom augment
#' @importFrom tidyr pivot_longer starts_with any_of
#' @export
augment_longer <- function(x, ...) {
  out <- augment(x, ...)
  out$.obs <- rownames(out)
  response <- response_var(x)

  return(pivot_longer(out, cols = !starts_with(".") & !any_of(response),
                      names_to = ".predictor_name",
                      values_to = ".predictor_value"))
}
