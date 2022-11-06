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
  classes <- unlist(lapply(df, class))

  factors <- (classes == "factor" | classes == "logical")

  if (!any(factors)) {
    return(df)
  }

  factor_names <- names(classes)[factors]
  cli_inform(c("Factor predictors were dropped:",
               "*" = "{.var {factor_names}}"))

  return(df[, !factors])
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
