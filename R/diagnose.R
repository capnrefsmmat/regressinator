#' Produce a lineup for a fitted model
#'
#' A lineup hides diagnostics among "null" diagnostics, i.e. the same
#' diagnostics calculated using models fit to data where all model assumptions
#' are correct. For each null diagnostic, `model_lineup()` simulates new
#' responses from the model using the fitted covariate values and the model's
#' error distribution, link function, and so on. Hence the new response values
#' are generated under ideal conditions: the fitted model is true and all
#' assumptions hold. `decrypt()` reveals which diagnostics are the true
#' diagnostics.
#'
#' To generate different kinds of diagnostics, the user can provide a custom
#' `fn`. The `fn` should take a model fit as its argument and return a data
#' frame. For instance, the data frame might contain one row per observation and
#' include the residuals and fitted values for each observation; or it might be
#' a single row containing a summary statistic or test statistic.
#'
#' `fn` will be called on the original `fit` provided. Then
#' `parametric_boot_distribution()` will be used to simulate data from the model
#' fit `nsim - 1` times, refit the model to each simulated dataset, and run `fn`
#' on each refit model. The null distribution is conditional on X, i.e. the
#' covariates used will be identical, and only the response values will be
#' simulated. The data frames are concatenated with an additional `.sample`
#' column identifying which fit each row came from.
#'
#' When called, this function will print a message such as
#' `decrypt("sD0f gCdC En JP2EdEPn ZY")`. This is how to get the location of the
#' true diagnostics among the null diagnostics: evaluating this in the R console
#' will produce a string such as `"True data in position 5"`.
#'
#' # Model limitations
#'
#' Because this function uses S3 generic methods such as `simulate()` and
#' `update()`, it can be used with any model fit for which methods are provided.
#' In base R, this includes `lm()` and `glm()`.
#'
#' The model provided as `fit` must be fit using the `data` argument to provide
#' a data frame. For example:
#'
#' ```
#' fit <- lm(dist ~ speed, data = cars)
#' ```
#'
#' When simulating new data, this function provides the simulated data as the
#' `data` argument and re-fits the model. If you instead refer directly to local
#' variables in the model formula, this will not work. For example, if you fit a
#' model this way:
#'
#' ```
#' # will not work
#' fit <- lm(cars$dist ~ cars$speed)
#' ```
#'
#' It will not be possible to refit the model using simulated datasets, as that
#' would require modifying your environment to edit `cars`.
#'
#' @param fit A model fit to data, such as by `lm()` or `glm()`
#' @param fn A diagnostic function. The function's first argument should be the
#'   fitted model, and it must return a data frame. Defaults to
#'   `broom::augment()`, which produces a data frame containing the original
#'   data and additional columns `.fitted`, `.resid`, and so on. To see a list
#'   of model types supported by `broom::augment()`, and to find documentation
#'   on the columns reported for each type of model, load the `broom` package
#'   and use `methods(augment)`.
#' @param nsim Number of total diagnostics. For example, if `nsim = 20`, the
#'   diagnostics for `fit` are hidden among 19 null diagnostics.
#' @param ... Additional arguments passed to `fn` each time it is called.
#' @return A data frame (tibble) with columns corresponding to the columns
#'   returned by `fn`. The additional column `.sample` indicates which set of
#'   diagnostics each row is from. For instance, if the true data is in position
#'   5, selecting rows with `.sample == 5` will retrieve the diagnostics from
#'   the original model fit.
#' @importFrom broom augment
#' @importFrom nullabor lineup
#' @importFrom stats simulate update
#' @importFrom tibble as_tibble
#' @seealso [parametric_boot_distribution()] to simulate draws by using the
#'   fitted model to draw new response values; [sampling_distribution()] to
#'   simulate draws from the population distribution, rather than from the model
#' @examples
#' fit <- lm(dist ~ speed, data = cars)
#' model_lineup(fit, nsim = 5)
#'
#' resids_vs_speed <- function(f) {
#'   data.frame(resid = residuals(f),
#'              speed = model.frame(f)$speed)
#' }
#' model_lineup(fit, fn = resids_vs_speed, nsim = 5)
#'
#' @references Buja et al. (2009). Statistical inference for exploratory data
#'   analysis and model diagnostics. *Philosophical Transactions of the Royal
#'   Society A*, 367 (1906), pp. 4361-4383. \doi{10.1098/rsta.2009.0120}
#'
#'
#' Wickham et al. (2010). Graphical inference for infovis. *IEEE Transactions on
#' Visualization and Computer Graphics*, 16 (6), pp. 973-979.
#' \doi{10.1109/TVCG.2010.161}
#' @export
model_lineup <- function(fit, fn = augment, nsim = 20, ...) {
  true <- fn(fit, ...)
  check_fn_output(true)

  simulated_diagnostics <- parametric_boot_distribution(fit, fn = fn,
                                                        nsim = nsim - 1, ...)

  return(as_tibble(lineup(true = true, samples = simulated_diagnostics,
                          n = nsim)))
}

#' @importFrom rlang caller_env
#' @importFrom cli cli_abort
check_fn_output <- function(x) {
  if (!inherits(x, "data.frame")) {
    cli_abort(
      c("Diagnostic function {.arg fn} must return a data frame or tibble",
        "x" = "{.arg fn} returned a result of class {.cls {class(x)}}"),
      class = "regressinator_diagnostic_class",
      call = caller_env()
    )
  }
}

#' Decrypt message giving the location of the true plot in a lineup
#'
#' Decrypts the message printed by `model_lineup()` indicating the location of
#' the true diagnostics in the lineup.
#'
#' @importFrom nullabor decrypt
#' @usage decrypt(...)
#' @param ... Message to decrypt, specifying the location of the true
#'   diagnostics
#' @return The decrypted message.
#' @export decrypt
#' @name decrypt
NULL

#' Simulate the distribution of estimates by parametric bootstrap
#'
#' Repeatedly simulates new response values by using the fitted model, holding
#' the covariates fixed. By default, refits the same model to each simulated
#' dataset, but an alternative model can be provided. Estimates, confidence
#' intervals, or other quantities are extracted from each fitted model and
#' returned as a tidy data frame.
#'
#' The default behavior samples from a model and refits the same model to the
#' sampled data; this is useful when, for example, exploring how model
#' diagnostics look when the model is well-specified. Another common use of the
#' parametric bootstrap is hypothesis testing, where we might simulate from a
#' null model and fit an alternative model to the data, to obtain the null
#' distribution of a particular estimate or statistic. Provide `alternative_fit`
#' to have a specific model fit to each simulated dataset, rather than the model
#' they are simulated from.
#'
#' Only the response variable from the `fit` (or `alternative_fit`, if given) is
#' redrawn; other response variables in the population are left unchanged from
#' their values in `data`.
#'
#' @inheritSection model_lineup Model limitations
#'
#' @param fit A model fit to data, such as by `lm()` or `glm()`, to simulate new
#'   response values from.
#' @param alternative_fit A model fit to data, to refit to the data sampled from
#'   `fit`. Defaults to `fit`, but an alternative model can be provided to
#'   examine its behavior when `fit` is the true model.
#' @param data Data frame to be used in the simulation. Must contain the
#'   predictors needed for both `fit` and `alternative_fit`. Defaults to the
#'   predictors used in `fit`.
#' @param fn Function to call on each new model fit to produce a data frame of
#'   estimates. Defaults to `broom::tidy()`, which produces a tidy data frame of
#'   coefficients, estimates, standard errors, and hypothesis tests.
#' @param nsim Number of total simulations to run.
#' @param ... Additional arguments passed to `fn` each time it is called.
#' @return A data frame (tibble) with columns corresponding to the columns
#'   returned by `fn`. The additional column `.sample` indicates which fit each
#'   row is from.
#' @seealso [model_lineup()] to use resampling to aid in regression diagnostics;
#'   [sampling_distribution()] to simulate draws from the population
#'   distribution, rather than the null
#' @importFrom broom tidy
#' @importFrom insight get_data
#' @importFrom purrr map_dfr
#' @importFrom tibble as_tibble
#' @examples
#' # Bootstrap distribution of estimates:
#' fit <- lm(mpg ~ hp, data = mtcars)
#' parametric_boot_distribution(fit, nsim = 5)
#'
#' # Bootstrap distribution of estimates for a quadratic model, when true
#' # relationship is linear:
#' quad_fit <- lm(mpg ~ poly(hp, 2), data = mtcars)
#' parametric_boot_distribution(fit, quad_fit, nsim = 5)
#'
#' # Bootstrap distribution of estimates for a model with an additional
#' # predictor, when it's truly zero. data argument must be provided so
#' # alternative fit has all predictors available, not just hp:
#' alt_fit <- lm(mpg ~ hp + wt, data = mtcars)
#' parametric_boot_distribution(fit, alt_fit, data = mtcars, nsim = 5)
#' @export
parametric_boot_distribution <- function(fit, alternative_fit = fit,
                                         data = get_data(fit),
                                         fn = tidy, nsim = 100, ...) {
  if (length(nsim) > 1 || nsim <= 0 || nsim %% 1 != 0) {
    cli_abort(c("Number of simulations must be a positive integer",
                "x" = "Received {.arg nsim} = {.val {nsim}}"))
  }

  check_data_arg(fit)

  # ensure fit uses the same predictors that alternative_fit will get
  orig_data <- data

  env <- current_env()

  # Check that update() works at all. If it works this time, it'll likely work
  # in future iterations, so we only prettify this error.
  tryCatch(
    fit <- update(fit, data = data),
    error = function(e) {
      cli_abort(
        c("Failed to re-fit the provided model {.arg fit}",
          x = conditionMessage(e)),
        class = "regressinator_update",
        call = env
      )
    }
  )


  simulated_ys <- simulate(fit, nsim = nsim)

  response <- response_var(alternative_fit)

  out <- map_dfr(
    seq_len(ncol(simulated_ys)),
    function(ii) {
      sim_data <- orig_data
      sim_data[, response] <- simulated_ys[, ii]

      sim_fit <- update(alternative_fit, data = sim_data)

      diagnostics <- fn(sim_fit, ...)

      check_fn_output(diagnostics)

      diagnostics$.n <- ii
      return(diagnostics)
    }
  )

  return(as_tibble(out))
}

#' Simulate the sampling distribution of estimates from a population
#'
#' Repeatedly refits the model to new samples from the population, calculates
#' estimates for each fit, and compiles a data frame of the results.
#'
#' To generate sampling distributions of different quantities, the user can
#' provide a custom `fn`. The `fn` should take a model fit as its argument and
#' return a data frame. For instance, the data frame might contain one row per
#' estimated coefficient and include the coefficient and its standard error; or
#' it might contain only one row of model summary statistics.
#'
#' @inheritSection model_lineup Model limitations
#'
#' @param fit A model fit to data, such as by `lm()` or `glm()`, to refit to
#'   each sample from the population.
#' @param data Data drawn from a `population()`, using `sample_x()` and possibly
#'   `sample_y()`. The `population()` specification is used to draw the samples.
#' @param fn Function to call on each new model fit to produce a data frame of
#'   estimates. Defaults to `broom::tidy()`, which produces a tidy data frame of
#'   coefficients, estimates, standard errors, and hypothesis tests.
#' @param nsim Number of simulations to run.
#' @param fixed_x If `TRUE`, the default, the predictor variables are held fixed
#'   and only the response variables are redrawn from the population. If
#'   `FALSE`, the predictor and response variables are drawn jointly.
#' @param ... Additional arguments passed to `fn` each time it is called.
#' @return Data frame (tibble) of `nsim + 1` simulation results, formed by
#'   concatenating together the data frames returned by `fn`. The `.sample`
#'   column identifies which simulated sample each row came from. Rows with
#'   `.sample == 0` come from the original `fit`.
#' @seealso [parametric_boot_distribution()] to simulate draws from a fitted
#'   model, rather than from the population
#' @importFrom broom tidy
#' @examples
#' pop <- population(
#'   x1 = predictor(rnorm, mean = 4, sd = 10),
#'   x2 = predictor(runif, min = 0, max = 10),
#'   y = response(0.7 + 2.2 * x1 - 0.2 * x2, error_scale = 4.0)
#' )
#'
#' d <- sample_x(pop, n = 20) |>
#'   sample_y()
#'
#' fit <- lm(y ~ x1 + x2, data = d)
#' # using the default fn = broom::tidy(). conf.int argument is passed to
#' # broom::tidy()
#' samples <- sampling_distribution(fit, d, conf.int = TRUE)
#' samples
#'
#' suppressMessages(library(dplyr))
#' # the model is correctly specified, so the estimates are unbiased:
#' samples |>
#'   group_by(term) |>
#'   summarize(mean = mean(estimate),
#'             sd = sd(estimate))
#'
#' # instead of coefficients, get the sampling distribution of R^2
#' rsquared <- function(fit) {
#'   data.frame(r2 = summary(fit)$r.squared)
#' }
#' sampling_distribution(fit, d, rsquared, nsim = 10)
#' @importFrom tibble as_tibble
#' @export
sampling_distribution <- function(fit, data, fn = tidy, nsim = 100,
                                  fixed_x = TRUE, ...) {
  if (length(nsim) > 1 || nsim <= 0 || nsim %% 1 != 0) {
    cli_abort(c("Number of simulations must be a positive integer",
                "x" = "Received {.arg nsim} = {.val {nsim}}"))
  }

  check_data_arg(fit)

  out <- fn(fit, ...)
  check_fn_output(out)

  out <- as_tibble(out)
  out$.sample <- 0

  for (b in seq_len(nsim)) {
    new_data <- if (fixed_x) {
      sample_y(data)
    } else {
      sample_y(sample_x(parent_population(data), nrow(data)))
    }
    new_fit <- fn(update(fit, data = new_data), ...)
    check_fn_output(new_fit)
    new_fit$.sample <- b

    out <- rbind(out, new_fit)
  }

  return(out)
}
