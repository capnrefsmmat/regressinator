#' Calculate diagnostics for a model, and produce a lineup of those diagnostics.
#'
#' A lineup hides the diagnostics among "null" diagnostics, i.e. the same
#' diagnostics calculated using models fit to data where all model assumptions
#' are correct. For each null diagnostic, `diagnose_model()` simulates new
#' responses from the model using the fitted covariate values and the model's
#' error distribution, link function, and so on. Hence the new response values
#' are generated under ideal conditions: the fitted model is true and all
#' assumptions hold. `decrypt()` reveals which diagnostics are the true
#' diagnostics.
#'
#' To generate different kinds of diagnostics, the user can provide a custom
#' `fn`. The `fn` should take a model fit as its argument and return a data
#' frame. For instance, the data frame might contain one row per observation and
#' include the residuals and fitted values for each observation; or it might
#' be a single row containing a summary statistic or test statistic.
#'
#' `fn` will be called on the original `fit` provided. Then `simulate()` will be
#' used to simulate data from the model fit `nsim - 1` times, the model will be
#' refit to each of these datasets, and `fn` will be run on each refit model.
#' The null distribution is conditional on X, i.e. the covariates used will be
#' identical, and only the response values will be simulated. The data frames
#' are concatenated with an additional `.sample` column identifying which fit
#' each row came from.
#'
#' When called, this function will print a message such as
#' `decrypt("sD0f gCdC En JP2EdEPn ZY")`. This is how to get the location of the
#' true diagnostics among the null diagnostics: evaluating this in the R console
#' will produce a string such as `"True data in position 5"`.
#'
#' Because `diagnose_model()` uses the S3 generic methods `simulate()` and
#' `update()`, it can be used with any model fit for which methods are provided.
#' In base R, this includes `lm()` and `glm()`.
#'
#' @param fit A model fit to data, such as by `lm()` or `glm()`
#' @param fn A diagnostic function. The function should take one argument, a
#'   fitted model, and return a data frame. The default is `broom::augment()`,
#'   which produces a data frame containing the original data and additional
#'   columns `.fitted`, `.resid`, and so on. To see a list of model types
#'   supported by `broom::augment()`, and to find documentation on the columns
#'   reported for each type of model, load the `broom` package and use
#'   `methods(augment)`.
#' @param nsim Number of total diagnostics. For example, if `nsim = 20`, the
#'   diagnostics for `fit` are hidden among 19 null diagnostics.
#' @return For `diagnose_model()`, a data frame with columns corresponding to
#'   the columns returned by `fn`. The additional column `.sample` indicates
#'   which set of diagnostics each row is from. For instance, if the true data
#'   is in position 5, selecting rows with `.sample == 5` will retrieve the
#'   diagnostics from the original model fit.
#'
#' For `decrypt()`, the decrypted message.
#' @importFrom broom augment
#' @importFrom nullabor lineup
#' @importFrom stats simulate update
#' @seealso [sampling_distribution()] to simulate draws from the population
#'   distribution, rather than the null
#' @examples
#' fit <- lm(dist ~ speed, data = cars)
#' diagnose_model(fit, nsim = 5)
#'
#' resids_vs_speed <- function(f) {
#'   data.frame(resid = residuals(f),
#'              speed = model.frame(f)$speed)
#' }
#' head(diagnose_model(fit, fn = resids_vs_speed, nsim = 5))
#' @export
diagnose_model <- function(fit, fn = augment, nsim = 20) {
  # TODO Wrap fn to ensure it always returns a data frame, regardless of its
  # return type
  true <- fn(fit)

  orig_data <- fit$model
  # get an empty data frame with the same column names
  simulated_diagnostics <- true[c(), ]
  simulated_ys <- simulate(fit, nsim = nsim - 1)
  for (ii in seq_len(ncol(simulated_ys))) {
    sim_data <- orig_data
    sim_data$y <- simulated_ys[, ii]

    sim_fit <- update(fit, data = sim_data)

    diagnostics <- fn(sim_fit)
    diagnostics$.n <- ii

    simulated_diagnostics <- rbind(simulated_diagnostics, diagnostics)
  }

  return(lineup(true = true, samples = simulated_diagnostics, n = nsim))
}

#' @rdname diagnose_model
#' @importFrom nullabor decrypt
#' @usage decrypt(...)
#' @param ... Message to decrypt, specifying the location of the true
#'   diagnostics
#' @export decrypt
#' @name decrypt
NULL
