## [Unreleased]

### Breaking Changes

- **No underscores in public names.** The keywords and arguments `plot_area`, `total_area`,
  `strata_area` and `point_area` are now `plotarea`, `totalarea`, `strataarea` and `pointarea`.
- **Arguments are never uppercase.** The population sizes `N`, `M`, `N1` and `N2` of `twostagesampling`,
  `independentoccasionssampling`, `completereplacementsampling`, `doublesampling` and
  `partialreplacementsampling` are now `npop`, `mpop`, `npop1` and `npop2` (they are positional, so the
  calls do not change). The result columns `N` and `M` keep their names.
- **The columns of every table are in the symbols of Péllico Netto & Brena (1997)**, in ASCII, in
  the form of the other packages of the ecosystem (lowercase for the estimates of the sample,
  uppercase for totals and values of the population, the index next to the symbol, `l` and `u`
  for the limits). They replace the English camelCase names:

  | Before | Now | Quantity |
  |:-------|:----|:---------|
  | `mean` | `xm` | mean ``\bar{x}`` per sampling unit |
  | `meanPerHectare` | `xha` | mean per hectare (acre) |
  | `coefficientOfVariation` | `cv` | coefficient of variation (%) |
  | `varianceOfMean`, `standardError` | `s2xm`, `sxm` | variance of the mean, standard error |
  | `absoluteError`, `relativeError` | `Ea`, `Er` | sampling errors |
  | `total`, `lowerLimit`, `upperLimit` | `X`, `Xl`, `Xu` | total and the limits of its interval |
  | `samplingFraction` | `f` | **now ``f = n/N``**; it used to hold the correction factor ``1 - f`` |
  | `requiredSampleSize`, `additionalSampleSize` | `nreq`, `nadd` | sampling intensity |
  | `area` | `A` | area of the forest |
  | `admissibleError`, `samplingInterval` | `e`, `m` | admissible error, number of lines |
  | `withinVariance`, `betweenVariance`, `totalVariance` | `s2d`, `s2e`, `s2x` | variances within, between clusters, total |
  | `intraclassCorrelation` | `r` | intracluster correlation |
  | `variance`, `n`, `mean` of the cluster tables | `s2i`, `n`, `xmi` | per cluster/start/primary unit |
  | `n`, `mean`, `variance`, `standardDeviation` of `auxiliarytable` | `nh`, `xmh`, `s2h`, `sh` | per stratum |
  | `areaShare`, `weightedStandardDeviation`, `weightedVariance` | `Wh`, `Whsh`, `Whs2h` | weight ``W_h`` and its products |
  | `sampleSizeByStratum` (and `required`/`additional`), `areaByStratum` | `nh`, `nhreq`, `nhadd`, `Ah` | tuples of the result table |
  | `Source`, `DOF`, `SS`, `MS`, `F`, `p(>F)` of `anova` | `source`, `gl`, `SQ`, `QM`, `F`, `p` | analysis of variance |
  | `treeCount`, `treesPerHectare`, `basalAreaPerHectare` | `n`, `N`, `G` | per point or line |
  | `basalAreaTotal` | `Gt` | total basal area |
  | `distance` of the Prodan point table | `L` | distance of the k-th tree |
  | `...SD` of `prodansampling` | `sN`, `sG`, `sxha` | standard deviation among points |
  | `averageRatio`, `sumPredicted` | `Rm`, `K` | `threepsampling` |
  | `regressionSlope`, `residualVariance`, `temporaryCount` | `b`, `s2yx`, `u` | `doublesampling` |
  | `regressionEstimate`, `temporaryMean`, `regressionWeight` | (gone, see below) | `partialreplacementsampling` |
  | `meanGrowth`, `totalGrowth`, and the limits of `change` | `dm`, `D`, `Dl`, `Du` | change between occasions (`sdm` is its standard error) |

### Documentation

- **The docstring of every design has the complete formulas** of Péllico Netto & Brena (1997),
  with the number of the equation, the symbols of the columns and the notation of the book: the
  mean, the variance, the coefficient of variation, the variance of the mean with the correction of
  the finite population, the standard error, the sampling errors, the total and its confidence
  interval, and the sampling intensity (solved iteratively), plus what is specific of each design
  (the Satterthwaite degrees of freedom and the Neyman allocation of the stratified sampling, the
  analysis of variance between and within clusters, the two stages, the regression of the double
  sampling, the weights of the partial replacement, the factors of the horizontal line sampling).
- The module docstring has the table of the notation, the README and the tutorial use the new
  columns, and the tutorial now has the horizontal line (Strand), Prodan and 3P designs.
- The examples of the docstrings pass the areas with units, and their outputs were run.

### Aligned with the book

Every design now applies the formulas of Péllico Netto & Brena (1997), and the test suite reproduces
the worked example of each chapter (5 to 11) with the data of the book, which were checked against the
sums and statistics that the book prints:

- `clustersampling`, `multistartsystematicsampling`: the finite-population correction multiplies only
  the variance between clusters, ``s_{ar{x}}^2 = (1 - f) s_e^2/n + s_d^2/(nM)`` (Eq. 9.17). It used to
  multiply the whole variance of the mean. The tables have the new column `s2xm`.
- `twostagesampling`: the variance of the mean is Eq. 8.14, ``(1 - f_1) s_e^2/n + (1 - f_2) s_d^2/(nm)``,
  each correction dropped when negligible (Eqs. 8.15 and 8.16); the factor ``f_1`` of the second term
  of Cochran's formula is gone. The table has the new column `s2xm`.
- `independentoccasionssampling`, `completereplacementsampling`, `doublesampling`,
  `partialreplacementsampling`: the degrees of freedom of ``t`` for the change are
  ``(n_1 - 1) + (n_2 - 1)``, one more than before.
- `doublesampling`: the variance of the change is Eq. 11.75, which uses the correlation of the
  permanent plots, instead of the variance of the estimate of occasion 2 alone.
- `partialreplacementsampling` **was rewritten with the estimators of Ware & Cunia (1962)** (Eqs.
  11.85 to 11.109): the constants ``a`` and ``c``, the current mean ``ar{y}``, its variance (11.86),
  and the direct estimate of the change with ``A`` and ``B`` (11.103 to 11.107). The previous combination
  of a regression estimate and the mean of the new plots (Cochran) is gone. The columns of `occasion2`
  are now `xm`, `sxm`, `Ea`, `Er`, `xha`, `X`, `Xl`, `Xu`, `a`, `c`, `r`, `m`, `u`, `n`, `N` (here `n` is the
  number of new plots, as in the book, and `v`, `yr`, `yn`, `b` are gone). It needs at least two plots of
  each of the groups `u` and `n`.
- `stratifiedsampling`: `cv` is the coefficient of variation of the stratified population,
  ``s_{st}/ar{x}_{st}`` (it was the relative standard error); the second term of the variance of the
  mean is dropped when the sampling fraction is negligible (Eq. 6.13); the degrees of freedom of
  Satterthwaite are rounded up as in the book; the required plots are distributed proportionally to the
  weights (Eqs. 6.27 to 6.29), as in the first stratification of the book, and the new keyword
  `allocation=:neyman` gives Neyman's optimum allocation (Eqs. 6.33 to 6.35). Every stratum needs two
  plots.
- `systematicsampling`: `cv` is the coefficient of variation of the plots, from their sample variance
  (it was the relative standard error of the mean).

### Differences from the book that remain

- **The tabulated t.** The book reads Student's t from a table with the degrees of freedom rounded
  (``t = 1.99`` for 69 degrees of freedom, ``2.0`` for 55, ...), the package computes it exactly, so the
  errors, `Ea` and `Er`, and the limits agree with the book in the first two or three digits.
- **The degrees of freedom of the cluster designs.** The examples of the book use ``t = 1.99`` (cluster,
  25 clusters), ``1.996`` (random starts, 9 starts) and ``1.96`` (two-stage, 36 primary units), which are
  the ``t`` of a large number of degrees of freedom; the package uses ``n - 1``, the clusters or primary
  units being the sampling units, so the errors are larger (2.064, 2.306 and 2.030). The variances, the
  means, the totals and the sampling intensities are those of the book.
- **The sampling intensity** is solved iteratively until it stops changing and rounded up, while the book
  stops after two or three steps and rounds to the nearest integer in some of them (18.05 plots of the
  stratified example become 18, the package gives 19). The other examples give the number of the book
  (70, 17, 13 and 36).
- **Errors of the book.** The book prints MQe = 2 218.968 in the two-stage example, but its own means give
  2 193.009, so the package gets ``s_e^2 = 276.228`` and ``s_{ar{x}}^2 = 7.390`` where the book prints
  280.555 and 7.467; the cell 32 x III of that table (125,77) and the unit 34 of Quadro 11.3 (181,8) are
  misprints, since the means and sums of the book need 35,92 and 180,8.
- **Not computed:** the sampling intensity of the successive occasions (Eqs. 11.82, 11.83 and the one of
  Sec. 11.4), the ratio estimator of Sec. 5.5 and the optimum number of subunits per cluster or primary
  unit.
- The data of some examples of the previous documentation had transcription errors (a plot of the cluster
  example, a cell of the systematic example, another of the random-start example); they were replaced by
  the data of the book. The systematic example used the 5 columns of Quadro 7.1 as lines; the lines of the
  book are its 15 rows.
