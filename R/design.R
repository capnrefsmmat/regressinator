#' Define a finite population with potential outcomes
#' TODO
potential_outcomes <- function() { }


#' Intervene to set factors according to an experimental design
#'
#' Specifies values of predictor variables according to a chosen experimental
#' design, so that simulated experiments can be conducted.
#'
#' The regressinator supports two paradigms of experimental design:
#'
#' - **Super-population inference.** In this paradigm, experimental units are
#' considered to be a random sample from a larger population. The experiment
#' fixes values of certain predictors (the experimental factors), and the
#' response variables are random because the units are sampled from the
#' population.
#'
#' - **Finite sample inference.** In this paradigm, the experimental units in
#' the study are the population; there is no super-population from which they
#' are drawn. Their potential outcomes (under each level of treatment) are fixed
#' in advance of the experiment. Randomness comes from the random assignment
#' process.
#'
#' Most experimental design textbooks consider super-population inference, since
#' it aligns with how inference is usually presented in linear regression.
#' Causal inference texts often use finite sample inference, since it fits the
#' potential outcomes framework well. See References below for further detail.
#'
#' # Super-population inference
#'
#' When given a population object from `population()`, `design_x()` supports
#' super-population inference, as used in most experimental design textbooks.
#' The `design` argument is a data frame giving the experimental design: each
#' row is an experimental unit, and each column is a treatment factor. This sets
#' the number of observations and their treatment assignments. `sample_y()` then
#' samples response values from the population.
#'
#' # Finite sample inference
#'
#' TODO
#'
#' @references Imbens and Rubin (2015). *Causal Inference for Statistics,
#' Social, and Biomedical Sciences: An Introduction*. Cambridge University
#' Press.
#'
#' Ding, Li, and Miratrix (2017). "Bridging finite and super population causal
#' inference", *Journal of Causal Inference*, 5:2. \doi{10.1515/jci-2016-0027}
#'
#' @param population Population object to draw response values from.
#' @param ... TODO
#' @export
design_x <- function(population, ...) {
  UseMethod("design_x")
}

#' @param design Data frame specifying the experimental design. Each column is
#'   one experimental factor; each row is one unit in the experiment.
#' @rdname design_x
#' @examples
#' # Super-population inference
#' pop <- population(
#'   treatment = predictor(rfactor, levels = c("treatment", "control")),
#'   response = response(by_level(treatment, treatment = 10, control = 0),
#'                       error_scale = 1)
#' )
#'
#' des <- data.frame(treatment = c(rep("treatment", 5), rep("control", 5)))
#'
#' # Repeated samples will draw new responses from the population
#' design_x(pop, des) |> sample_y()
#' @export
design_x.population <- function(population, design, ...) {
  n <- nrow(design)

  output <- sample_x(population, n)

  factor_names <- names(design)
  response_names <- names(population_response(population))

  for (fact in factor_names) {
    if (!(fact %in% names(output))) {
      cli_abort(c("{.arg design} argument must be a data frame containing assignments for predictors in the population",
                  "i" = "Population does not contain this predictor: {.val {fact}}"),
                class = "regressinator_design_preds")
    }

    output[[fact]] <- design[[fact]]
  }

  class(output) <- c("designed_sample", class(output))

  return(output)
}
