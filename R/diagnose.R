#' Calculate diagnostics for a model, and produce a lineup of those diagnostics.
#'
#' A lineup hides the diagnostics among "null" diagnostics, i.e. the same
#' diagnostics calculated using models fit to data where all model assumptions
#' are correct.
#'
#' Specifically, the null case conditions on X, meaning the covariates are
#' identical. For each null diagnostic, `diagnose_model()` simulates new
#' responses from the model using the fitted covariate values and the model's
#' error distribution, link function, and so on. Hence the new response values
#' are generated under ideal conditions: the fitted model is true and all
#' assumptions hold.
#'
#' When called, this function will print a message such as
#' `decrypt("sD0f gCdC En JP2EdEPn ZY")`. This is how to get the location of the
#' true diagnostics among the null diagnostics: evaluating this in the R console
#' will produce a string such as `"True data in position 5"`.
#'
#' @param fit A model fit to data, such as by `lm()` or `glm()`
#' @param fn A diagnostic function. The function should take one argument, a
#'   fitted model, and return a data frame. For example, to evaluate residuals
#'   versus a particular covariate, the function should return a data frame
#'   where one column is that covariate, and the other column is the
#'   corresponding residuals. The default is `broom::augment()`, which produces
#'   a data frame containing the original data and additional columns `.fitted`,
#'   `.resid`, and so on. To see a list of model types supported by
#'   `broom::augment()`, and to find documentation on the columns reported for
#'   each type of model, load the `broom` package and use `methods(augment)`.
#' @param n Number of total diagnostics. For example, if `n = 20`, the
#'   diagnostics for `fit` are hidden among 19 null diagnostics.
#' @return A data frame with columns corresponding to the columns returned by
#'   `fn`. The additional column `.sample` indicates which set of diagnostics
#'   each row is from. For instance, if the true data is in position 5,
#'   selecting rows with `.sample == 5` will retrieve the diagnostics from the
#'   original model fit.
#' @importFrom broom augment
#' @importFrom nullabor lineup
#' @importFrom stats simulate update
#' @export
diagnose_model <- function(fit, fn = augment, n = 20) {
  # TODO Wrap fn to ensure it always returns a data frame, regardless of its
  # return type
  true <- fn(fit)

  orig_data <- fit$model
  # get an empty data frame with the same column names
  simulated_diagnostics <- true[c(), ]
  simulated_ys <- simulate(fit, nsim = n - 1)
  for (ii in seq_len(ncol(simulated_ys))) {
    sim_data <- orig_data
    sim_data$y <- simulated_ys[, ii]

    sim_fit <- update(fit, data = sim_data)

    diagnostics <- fn(sim_fit)
    diagnostics$.n <- ii

    simulated_diagnostics <- rbind(simulated_diagnostics, diagnostics)
  }

  return(lineup(true = true, samples = simulated_diagnostics, n = n))
}

#' Reveal the position of the real data in lineups
#'
#' `diagnose_model()` will print an encrypted message stating where the real
#' data is in the lineup it creates. For instance, it might print
#' `decrypt("sD0f gCdC En JP2EdEPn ZY")`, which might decrypt to a string such
#' as `"True data in position 5"`.
#'
#' See the nullabor package documentation for more details on `decrypt()` and
#' its lineup system.
#'
#' @importFrom nullabor decrypt
#' @export decrypt
#' @name decrypt
NULL
