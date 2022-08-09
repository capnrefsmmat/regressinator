
<!-- WARNING: README.md is generated from README.Rmd; edit that
     instead, then use rmarkdown::render("README.Rmd") to regenerate -->

# the regressinator

The regressinator is a pedagogical tool for conducting simulations of
regression analyses and diagnostics. It can:

-   Simulate populations with predictor variables from arbitrary
    distributions
-   Simulate response variables that are function of the predictor
    variables plus error, or are drawn from a distribution related to
    the predictors
-   Given a model fit to simulated data, generate new simulated data
    based on the model fit
-   Facilitate lineup plots comparing diagnostics on the fitted model to
    diagnostics where all model assumptions are met.

Or, in other words: We can specify a population with a known
relationship between predictors in response, simulate data from that
population, fit a model to that data, and simulate data from the model.
If the fitted model is misspecified, we can simulate data from the model
(where it is *correctly* specified), and make a lineup to compare the
fitted model’s diagnostics to null diagnostics.

Here’s a simple regression example:

``` r
library(regressinator)
library(ggplot2)

# true relationship is nonlinear
pop <- population(
  x = predictor("runif", min = -5, max = 5),
  y = response(10 + 0.7 * x**2, family = gaussian(), error_scale = 2)
)

nonlin_sample <- pop |>
  sample_x(n = 20) |>
  sample_y()

fit <- lm(y ~ x, data = nonlin_sample)

diagnose_model(fit) |>
  ggplot(aes(x = x, y = .resid)) +
  geom_point() +
  facet_wrap(~ .sample) +
  labs(x = "x", y = "Residual")
```

    ## decrypt("23eg MuPu NE KwWNPNwE Ft")

<img src="man/figures/README-example-regression-lineup-1.png" width="672" />

Each of the 20 residual plots represents a simple regression fit. One of
the regressions was fit to the simulated data with a nonlinear
relationship; the others were fit to data simulated *from the linear
model*, using the same model. (That is, they were simulated from `fit`,
and then the same model was refit to them.) So one of these residual
plots shows the residuals of a linear model fit to nonlinear data, and
the others show the residuals of a linear model fit to *linear* data.

The premise of the regressinator is that this kind of simulation can be
a valuable teaching tool when students are learning to interpret
regression diagnostics and determine how different kinds of
misspecification might affect model fits.

Check the Get Started guide for more detail.

## Inspirations

Visual inference – hiding a plot of real data among “null plots” – has
been around for quite a while, and I certainly did not invent this idea
myself. For example:

-   Buja, A., Cook, D., Hofmann, H., Lawrence, M., Lee, E. K.,
    Swayne, D. F., & Wickham, H. (2009). [Statistical inference for
    exploratory data analysis and model
    diagnostics](https://doi.org/10.1098/rsta.2009.0120). *Philosophical
    Transactions of the Royal Society A*, 367(1906), 4361–4383.

-   Wickham, H., Cook, D., Hofmann, H., & Buja, A. (2010). [Graphical
    inference for infovis](https://doi.org/10.1109/TVCG.2010.161). *IEEE
    Transactions on Visualization and Computer Graphics*, 16(6),
    973–979.

The idea of using visual inference to teach regression diagnostics is
also not new:

-   Loy, A. (2021). [Bringing Visual Inference to the
    Classroom](https://doi.org/10.1080/26939169.2021.1920866). *Journal
    of Statistics and Data Science Education*, 29(2), 171-182.

The regressinator simply adds an easy-to-use framework to allow all
kinds of teaching activities to be constructed quickly. Instructors can
design example to display in lecture, or students with R experience can
run interactive examples and explore different situations. Ideally, if a
student asks “But what if *\[some problem with the model\]* happens?”,
you should be able to reply with a quick simulation.

## Compared to other packages

There have been several past efforts to support pedagogical simulation
and lineup plots in R:

-   [nullabor](https://cran.r-project.org/package=nullabor), which
    supports lineup plots. The regressinator uses nullabor underneath
    when building lineups.
-   [mosaic](https://cran.r-project.org/package=mosaic), part of
    [Project MOSAIC](http://www.mosaic-web.org/), provides a simple set
    of functions for doing EDA and basic inferential statistics,
    including the `do()` function for easy simulation without loops.
-   [infer](https://infer.tidymodels.org/), part of the
    [tidymodels](https://www.tidymodels.org/) framework, can conduct
    simulation-based inference for proportions, means, regression
    slopes, and other estimates commonly used in statistics courses,
    using only a few simple functions.

Unlike these packages, the regressinator provides a simple tool for
specifying a *population* and sampling from it, rather than conducting
bootstrapping or permutation on an observed dataset. This makes the
regressinator suitable for, say, exploring the properties of regression
estimates and diagnostics in known populations, but less suitable for
simulation-based hypothesis testing.

Unlike infer, the regressinator does not wrap R statistical methods or
provide its own inference functions. Users must use `lm()`, `glm()`, or
whatever other methods they need for their modeling. This makes the
regressinator less suitable for introductory courses where extra
complexity should be hidden away from students, but more suitable for
more advanced work: as students advance to more complex models provided
by other packages, they can use those models in the regressinator,
without any special support being required.

A useful counterpart to the regressinator might be
[rsample](https://rsample.tidymodels.org/), a general framework for
methods that resample from the observed data, such as bootstrapping. In
the same way that the regressinator supports general-purpose simulation
from the population without hard-coding specific use cases, rsample
supports resampling and cross-validation in a general way that could be
used for any kind of modeling, not just models built into the package.
