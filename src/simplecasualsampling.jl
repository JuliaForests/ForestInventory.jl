# A plot occupying a full reference area (1 ha in the metric system, 1 ac in the
# imperial one) is treated as effectively infinite regardless of the sampling fraction,
# matching the convention used by the original tool this design was ported from.
_isreferencearea(a::Area) = a == 1 * referencearea(a)

"""
    simpleCasualSampling(value::AbstractVector, plot_area::Area, total_area::Area; e::Real=10, α::Real=0.95)
    simpleCasualSampling(value::AbstractVector, plot_area::Real, total_area::Real; kwargs...)

Performs simple random sampling for forest inventory analysis with specified plot area and total area.

# Description

Calculates various statistical parameters for a forest inventory using simple random
sampling methodology. It is designed to estimate the total of some measured quantity
(volume, basal area, tree count, ...) within a forest area by analyzing measurements from
sample plots.

Simple random sampling is a statistical method where each unit (in this case, a plot) has an equal chance of being selected. This method is commonly used in forest inventories to estimate parameters like the mean per hectare, the total, and confidence intervals.

# Arguments

- `value`: vector of measurements from each sampled plot (volume, basal area, tree count,
  or anything else being sampled). If given as a `Unitful` quantity, its unit is preserved
  through every reported quantity; if given as plain numbers, **no unit is assumed** and
  every reported quantity stays a plain number too — unlike other functions in this
  package, this design has no way to know what a bare number here represents.
- `plot_area`: the area of each sample plot, as an `Area` quantity. All plots are assumed to have the same area. A plain number is taken to be hectares.
- `total_area`: the total area of the forest or stand being inventoried, as an `Area` quantity. A plain number is taken to be hectares.
- `e::Real=10`: the desired relative error margin as a percentage (default is 10%). This represents the maximum acceptable error in the estimate relative to the true mean.
- `α::Real=0.95`: the confidence level for the statistical estimates (default is 95%). This determines the confidence interval for the total estimate.

# Returns

- `DataFrame`: a single-row report, one column per statistic. Every column that carries a physical unit is a real `Unitful` quantity (exactly like `dmetrics`/`standmetrics` elsewhere in this package), so `removeunits`/`restoreunits` work on it directly. Columns:
  - `mean`: mean value per plot.
  - `coefficientOfVariation`: coefficient of variation (%).
  - `varianceOfMean`: variance of the mean.
  - `standardError`: standard error of the mean.
  - `absoluteError`: absolute sampling error.
  - `relativeError`: relative sampling error (%).
  - `meanPerHectare`: estimated value per hectare.
  - `total`: estimated total.
  - `lowerLimit`/`upperLimit`: confidence interval for the total.
  - `population`: `"finite"` or `"infinite"`, the population classification used.
  - `samplingFraction`: finite-population correction factor.
  - `n`: number of measured plots.
  - `requiredSampleSize`: required number of plots for the target error `e`.
  - `additionalSampleSize`: additional plots still needed (`0` if `n ≥ requiredSampleSize`).
  - `N`: number of possible plots in the population.
  - `area`: total area.

# Mathematical basis

The mean plot value and its coefficient of variation:
```math
\\bar{x} = \\frac{\\sum_{i=1}^{n} x_i}{n} \\qquad CV = \\frac{s}{\\bar{x}} \\times 100
```

For an infinite population, the variance of the mean is `s²/n`; for a finite population of
size `N = \\text{total\\_area}/\\text{plot\\_area}`, a finite-population correction `f = 1 - n/N` is applied:
```math
s^2_{\\bar{x}} = \\frac{s^2}{n} \\qquad \\text{or} \\qquad s^2_{\\bar{x}} = \\frac{s^2}{n}\\left(1 - \\frac{n}{N}\\right)
```

The required sample size for the desired error `e` is solved iteratively (see `_requiredsamplesize`), since the t quantile depends on the sample size through its own degrees of freedom.

# Technical description

Simple random sampling is the reference design against which every other sampling method
in this module is compared: stratified, cluster, and systematic sampling all reduce to it
in the degenerate case of a single stratum/cluster/full census. The population is
classified as infinite whenever the sampling fraction is negligible (`f ≥ 0.98`) or the
plot occupies a full reference area (1 ha or 1 ac), in which case the finite-population
correction is dropped from the variance of the mean.

# Examples

```julia-repl
julia> v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1]u"m^3";

julia> report = simpleCasualSampling(v, 0.05, 10; e=10, α=0.95)
1×17 DataFrame
 Row │ mean        coefficientOfVariation  varianceOfMean  standardError  absoluteError  relativeError  meanPerHectare  total       lowerLimit  upperLimit  population  samplingFraction  n      requiredSampleSize  additionalSampleSize  N      area
     │ Quantity…   Float64                 Quantity…       Quantity…      Quantity…      Float64        Quantity…       Quantity…   Quantity…   Quantity…   String      Float64           Int64  Int64                Int64                  Int64  Quantity…
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 441.691 m^3                   10.03     168.498 m^6    12.9807 m^3    28.9227 m^3           6.55  8833.82 m^3 ha^-1  88338.2 m^3  82553.6 m^3  94122.7 m^3  finite                 0.945     11                   7                      0    200  10 ha

julia> report.mean
441.6909090909091 m^3
```
"""
function simpleCasualSampling(value::AbstractVector, plot_area::Area, total_area::Area;
  e::Real=10, α::Real=0.95)

  n = length(value)
  n > 1 || throw(ArgumentError("At least two sampled plots are required."))

  N = round(Int, ustrip(uconvert(NoUnits, total_area / plot_area)))
  f = 1 - n / N
  x̅ = mean(value)
  cv = std(value) / x̅ * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  infinite = _isinfinitepopulation(n, N) || _isreferencearea(plot_area)
  if infinite
    population = "infinite"
    s²x̅ = var(value) / n
    requiredplots = _requiredsamplesize(cv^2, n, e, α)
  else
    population = "finite"
    s²x̅ = (var(value) / n) * f
    requiredplots = _requiredsamplesize(cv^2, n, e, α, N)
  end
  sx̅ = sqrt(s²x̅)
  missingplots = n > requiredplots ? 0 : requiredplots - n

  absoluteerror = Ttab * sx̅
  relativeerror = absoluteerror / x̅ * 100
  hectarevalue = x̅ / uconvert(referencearea(plot_area), plot_area)
  totalvalue = x̅ * N
  ciupper = totalvalue + N * absoluteerror
  cilower = totalvalue - N * absoluteerror

  return DataFrame(
    mean=x̅, coefficientOfVariation=round(cv, digits=2), varianceOfMean=s²x̅, standardError=sx̅,
    absoluteError=absoluteerror, relativeError=round(relativeerror, digits=2), meanPerHectare=hectarevalue,
    total=totalvalue, lowerLimit=cilower, upperLimit=ciupper, population=population,
    samplingFraction=round(f, digits=3), n=n, requiredSampleSize=requiredplots,
    additionalSampleSize=missingplots, N=N, area=total_area,
  )
end

simpleCasualSampling(value::AbstractVector, plot_area::Real, total_area::Real; kwargs...) =
  simpleCasualSampling(value, plot_area * AUNIT, total_area * AUNIT; kwargs...)
