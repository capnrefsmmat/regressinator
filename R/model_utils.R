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

#' Drop factor columns from a data frame
#'
#' Issues messages for the columns dropped.
#'
#' @param df Data frame to process
#' @return Data frame without columns that are factors
#' @keywords internal
#' @importFrom cli cli_inform
drop_factors <- function(df) {
  is_factor <- function(obj) {
    inherits(obj, c("factor", "logical", "character"))
  }
  factors <- unlist(lapply(df, is_factor))

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
      if (is.factor(df[, p])) {
        lev <- levels(df[, p])
        df[, p] <- factor(lev[1], levels = lev)
      } else if (is.logical(df[, p])) {
        df[, p] <- FALSE
      } else if (is.character(df[, p])) {
        # first factor level is the first string, in order
        df[, p] <- sort(unique(df[, p]))[1]
      } else {
        df[, p] <- 0
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
