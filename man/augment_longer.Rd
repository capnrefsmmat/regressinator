% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/residuals.R
\name{augment_longer}
\alias{augment_longer}
\title{Augment a model fit with residuals, in "long" format}
\usage{
augment_longer(x, ...)
}
\arguments{
\item{x}{A model fit object, such as those returned by \code{lm()} or \code{glm()}. See
the broom documentation for the full list of model types supported.}

\item{...}{Additional arguments passed to \code{broom::augment()}.}
}
\value{
A data frame (tibble) in similar form to those produced by
\code{broom::augment()}, but expanded to have one row per predictor per
observation. Columns \code{.predictor_name} and \code{.predictor_value} identify the
predictor and its value. An additional column \code{.obs} records the original
observation numbers so results can be matched to observations in the
original model data.
}
\description{
Use \code{broom::augment()} to augment a model fit with residual and fit
information, then reformat the resulting data frame into a "long" format with
one row per predictor per observation, to facilitate plotting of the result.
}
\details{
The name comes by analogy to \code{tidyr::pivot_longer()}, and the concept of long
versus wide data formats.
}
\section{Limitations}{
Factor predictors (as factors, logical, or character vectors) can't coexist
with numeric variables in the \code{.predictor_value} column. If there are some
numeric and some factor predictors, the factor predictors will automatically
be omitted. If all predictors are factors, they will be combined into one
factor with all levels. However, if a numeric variable is converted to factor
in the model formula, such as with \code{y ~ factor(x)}, the function cannot
determine the appropriate types and will raise an error. Create factors as
needed in the source data frame \emph{before} fitting the model to avoid this
issue.
}

\examples{
fit <- lm(mpg ~ cyl + disp + hp, data = mtcars)

# each observation appears 3 times, once per predictor:
augment_longer(fit)

fit <- lm(mpg ~ cyl + disp + I(disp^2) + hp, data = mtcars)

# each observation still appears 3 times, as disp and disp^2 are one
# predictor:
augment_longer(fit)
}
\seealso{
\code{\link[=partial_residuals]{partial_residuals()}}, \code{\link[=binned_residuals]{binned_residuals()}}
}
