# Utilities for examining models, formulas, and predictors

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

#' Check whether each column in a data frame is a factor
#'
#' @param df Data frame to process
#' @return Logical vector with same number of entries as columns in `df`. Each
#'   entry is `TRUE` if the corresponding column is a factor.
#' @keywords internal
factor_columns <- function(df) {
  is_factor <- function(obj) {
    inherits(obj, c("factor", "logical", "character"))
  }
  factors <- unlist(lapply(df, is_factor))

  return(factors)
}

#' Drop factor columns from a data frame
#'
#' Issues messages for the columns dropped.
#'
#' @param df Data frame to process
#' @return Data frame without columns that are factors
#' @keywords internal
#' @importFrom cli cli_inform
drop_factors <- function(df) {
  factors <- factor_columns(df)

  if (!any(factors)) {
    return(df)
  }

  factor_names <- names(df)[factors]
  cli_inform(c("Factor predictors were dropped:",
               "*" = "{.var {factor_names}}"))

  return(df[, !factors, drop = FALSE])
}

#' Get a prototype data frame for partial residuals
#'
#' All predictors, except the one we are calculating partial residuals for, are
#' set to 0 (or their baseline level, for factors).
#'
#' @param df data frame of predictors
#' @param predictor character vector identifying one predictor
#' @return prototype data frame
#' @keywords internal
prototype_for <- function(df, predictor) {
  for (p in names(df)) {
    if (p != predictor) {
      if (is.factor(df[[p]])) {
        lev <- levels(df[[p]])
        df[[p]] <- factor(lev[1], levels = lev)
      } else if (is.logical(df[[p]])) {
        df[[p]] <- FALSE
      } else if (is.character(df[[p]])) {
        # first factor level is the first string, in order
        df[[p]] <- sort(unique(df[[p]]))[1]
      } else {
        df[[p]] <- 0
      }
    }
  }
  return(df)
}

#' Detect transmutation in formulas, such as factor(), and raise an error
#'
#' We rely on predictors to occur in models with only one type (such as numeric
#' or factor), but the use of factor() would make it possible for a predictor to
#' appear both as a factor or as numeric. The use of factor() also makes it
#' harder to correctly detect the types of predictors, since the methods for
#' obtaining model predictors provide them before they are converted to factor,
#' not after. So we reject formulas that transmute types inside the formula,
#' such as with factor().
#'
#' Presently only factor() calls are rejected, but if other transmutations (such
#' as conversions to logical or numeric) prove to be problems, they can be
#' rejected as well.
#'
#' @param formula Model formula
#' @param call Environment in which to raise the error, defaulting to the
#'   calling environment. As this function is recursive, this reduces the
#'   complexity of backtraces.
#' @return No value. Raises an error if transmutation is present.
#' @importFrom cli cli_abort
#' @importFrom rlang is_formula is_call
#' @keywords internal
detect_transmutation <- function(formula, call = parent.frame()) {
  if (is_formula(formula)) {
    # look on the RHS of the formula
    detect_transmutation(formula[[3]], call)
  }
  if (is_call(formula, "factor")) {
    cli_abort(c("Model formula contains a call to {.fun factor}",
                "*" = "In term {.code {deparse(formula)}}",
                "i" = "Convert variables to factors before fitting the model"),
              class = "regressinator_transmutation_factor",
              call = call)
  }
  if (is.call(formula)) {
    sapply(as.list(formula),
           function(el) {
             detect_transmutation(el, call)
           })
  }
}

#' Check that the model fit uses the data argument to provide data
#'
#' We simulate by passing simulated data arguments to update(). If the original
#' fit does not use the data argument, and instead refers directly to variables
#' in the environment, the simulations will not behave as expected. This may
#' result in the "simulated" fits all using the original data, for instance.
#'
#' For example, in `lm(mtcars$mpg ~ mtcars$drat)`, simulating new data and
#' providing it in `data =` will not change the data used for fitting.
#'
#' Detect a missing `data` argument and abort. It is still possible to provide
#' `data` but also refer directly to the calling environment, but this is harder
#' to detect.
#'
#' @param fit A fitted model object, whose call is to be examined
#' @return No value. Raises an error if no `data` argument was used in `fit`.
#'
#' @importFrom cli cli_abort
#' @importFrom stats getCall
#' @keywords internal
check_data_arg <- function(fit) {
  call <- getCall(fit)

  if (!("data" %in% names(call))) {
    cli_abort(c("Model fit does not contain a {.arg data} argument; cannot refit with simulated data",
                "*" = "Model call was: {.code {deparse(call)}}",
                "i" = "Simulations work by updating {.arg data} argument and refitting",
                ">" = "Refit model with a formula referring to columns in {.arg data}"),
              class = "regressinator_data_arg")
  }
}
