# regressinator 0.2.1

- The documentation for `partial_residuals()` now gives the correct definition
  of partial residuals for GLMs. (The code gave the correct output already, but
  the stated definition was incorrect.)

- `parametric_boot_distribution()`, and by extension `model_lineup()`, now
  supports models fit with transformations. Previously, using transformations
  like `I(x^2)` in model formulas (or using `ns()` or other functions that
  transform predictors into matrix columns) would cause the model refitting to
  fail.

# regressinator 0.2.0

- New `augment_quantile()` and `augment_quantile_longer()` functions augment
  model training data with randomized quantile residuals.

- If a `predictor()` is multivariate (it generates a matrix of multiple columns
  when sampled from), `sample_x()` will now use the column names to set the
  predictor names appropriately, rather than simply numbering them.

- It is now recommended to provide a distribution function to `predictor()`,
  rather than a function name as a string. (For example, use `predictor(rnorm,
  mean = 1)` rather than `predictor("rnorm", mean = 1)`.) Using a string
  introduces problems when the distribution function is defined in an
  environment not on the search path when `sample_x()` needs to find it; by
  passing a function directly, `sample_x()` will always be able to find it.

- `partial_residuals()` and `binned_residuals()` now reject tidyselect syntax
  that tries to rename predictors, or that results in no predictors being
  selected. This syntax already caused the functions to fail with strange error
  messages, so it is now more explicitly rejected.

# regressinator 0.1.3

This version fixes several bugs that arose during classroom use.

- Simulation functions (`model_lineup()`, `parametric_boot_distribution()`, and
  `sampling_distribution()`) now check to determine if the model being simulated
  from was fit using the `data =` argument, and issue an error if it was not.
  The simulations work by calling `update(fit, data = ...)` with newly simulated
  data, and `update()` uses this to call the model fit function again with the
  specified `data =` argument. But if the model was fit without one, the
  argument is unused, and the simulations just reuse the original data.

  For example, if you fit this model:

  ```
  bad_fit <- lm(cars$dist ~ cars$speed)
  ```

  the simulation functions cannot work correctly because even with a different
  `data =` argument, the model fit will still refer to `cars`. The model should
  be fit like this:

  ```
  good_fit <- lm(dist ~ speed, data = cars)
  ```

  To prevent simulation problems, a suitable error is issued, so the user can
  refit the model correctly.

- `response()` now correctly detects when the `error_scale` argument was missing
  and issues the appropriate error.

- `augment_longer()` now supports models with factor predictors. If there are
  some factors and some continuous predictors, the factors are omitted from the
  result; if the predictors are all factors, they are kept.

- `parametric_boot_distribution()` now supports simulations when
  `alternative_fit` uses predictors that were not used in `fit`. Previously,
  these would fail because the simulated data frame only contained the
  predictors used in `fit`. Supply the new `data` argument to specify the data
  frame used in simulations.

# regressinator 0.1.2

First version released to CRAN.
