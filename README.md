# The `buildmer` package

`buildmer` is an `R` package written to simplify the process of testing whether the terms in your `lmer` (or equivalent) models make a significant contribution to the log likelihood, AIC, BIC, or explained deviance (the latter is not a formal statistical test, but informally tests if the model fit improved by at least an ounce of a percent). The aim of the package is to **fully automate model selection**, as is already possible for non-mixed regression models and lmerTest models using the `step` function.

In addition, the package (optionally) **determines the order of your predictors** by their contribution to the model fit. This is intended to mimic the results you might get from the stepwise functions available in, e.g., *SPSS* (although *SPSS* support stepwise elimination only for fixed-effect models...). In sum, the `buildmer` package aims to take the complex and time-consuming parts of the model fitting procedure out of your hands -- all you need to do is specify your intended maximal model and your dataset, and the package will take care of the rest.

**Nonconvergence** of models is handled properly by removing a random slope if the maximal model you have specified turns out to not converge. **p-values** are calculated for you using Wald z-scores. Better alternatives (Kenward-Roger, Satterthwaite) are available at the user's option (if the `lmerTest` package is installed).

The package supports the fitting of a wide variety of models, if the relevant packages are available. The following `buildmer` functions make it possible to fit the following types of models:
 * *buildmer*: `lm`, `glm`, `lmer`, `glmer`
 * *buildgamm4*: `gamm4` (package `gamm4`)
 * *buildbam*: `bam` (package `mgcv`)
 * *buildgam*: `gam` (package `mgcv`)
 * *buildgamm*: `gamm` (package `mgcv`)
 * *buildglmmTMB*: `glmmTMB` (package `glmmTMB`)
 * *buildmultinom*: `multinom` (package `nnet`)
 * *buildmertree*: `lmertree`, `glmertree` (package `glmertree`)
 * *buildlme*: `lme` (package `nlme`)
 * *buildclmm*: `clmm` (package `ordinal`)
 * *buildmer.nb*: `glm.nb` (package `MASS`) and `glmer.nb` (package `lme4`)
 * *buildcustom*: enables you to write your own wrapper function to use the buildmer term-reordering and elimination features with any type of model you want.

Automatic elimination of fixed, random, and/or smooth terms, is possible and enabled by default using the backward (default) or forward stepwise method. Bi-directional elimination is also possible, by passing e.g. `direction=c('forward','backward','forward')`, although I would not want to recommend doing this.

The intention of the `buildmer` package is to make your life as simple as:

```
library(buildmer)
model <- buildmer(Reaction ~ Days + (Days|Subject),lme4::sleepstudy)
```

In other words, the intention of the package is for you to specify the maximal model that you *would like* to fit, and let the package worry about whether or not your model converges, and whether or not your effect structure before, during, or after non-significant term removal needs to be passed to `gamm4`, `lmer`, or `lm` (or `glm`(`er`) or `{g,b}am`). More options are available; please see the documentation for details.

If certain terms need to be added together, you can construct your own buildmer formula object by calling `tab <- tabulate.formula(formula)`, modifying the `tab` to suit your needs, and passing it to the call to `buildmer()` (or another buildmer fitter) in the formula part, with an extra argument `dep` passing the dependent variable as a character string. This can be used to force certain parameters to only be evaluated together (see the `block` column in the `tab` object); this makes it possible to use `buildmer` to build models with diagonal covariance structures if using a specially-prepared data frame. See the vignette for details.

## Installation
```
remotes::install_github('cvoeten/buildmer')
```

## Known issues

Bugs will be fixed as they are uncovered. If you notice a bug, please file an issue for it and I will look into it if I have time.
