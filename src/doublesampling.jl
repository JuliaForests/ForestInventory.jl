"""
    doubleSampling(value1::AbstractVector, value2::AbstractVector, plot_area::Area, N::Integer; α::Real=0.95)
    doubleSampling(value1::AbstractVector, value2::AbstractVector, plot_area::Real, N::Integer; kwargs...)

Estimates growth between two inventory occasions using a regression estimator over a large first-occasion sample and a smaller remeasured (permanent) subset, generic to any sample size.

# Description

`n` plots are measured at occasion 1; only `m ≤ n` of them ("permanent" plots) are
remeasured at occasion 2 — the rest are "temporary", occasion-1-only plots that exist
purely to pin down `x̄₁` more precisely at low cost. The occasion-2 mean is then estimated
by regressing the permanent plots' occasion-2 value on their occasion-1 value and
applying that regression to the full, larger occasion-1 sample (Cochran, 1977, §12.9).

# Arguments

- `value1`: occasion-1 measurement for every plot (`n` entries). If given as a `Unitful` quantity, its unit is preserved; plain numbers are kept unitless throughout.
- `value2`: occasion-2 measurement, one entry per plot in the same order as `value1`; `missing` marks a temporary plot that was not remeasured.
- `plot_area`: the area of each plot, as an `Area` quantity, used only to report the per-hectare value. A plain number is taken to be hectares.
- `N::Integer`: total number of possible plots in the population.
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- [`SamplingReport`](@ref) with `occasion1` (simple-random-sampling report over all `n`
  plots, same columns as [`independentOccasionsSampling`](@ref)'s occasion tables),
  `occasion2` (the regression estimate) and `change` (the estimated growth). `occasion2`
  columns:
  - `mean`: the regression-estimated occasion-2 mean.
  - `regressionSlope`: the regression slope.
  - `residualVariance`: residual variance of the regression.
  - `standardError`, `absoluteError`, `relativeError`, `meanPerHectare`, `total`, `lowerLimit`, `upperLimit`: as in [`simpleCasualSampling`](@ref).
  - `m`: number of permanent (remeasured) plots.
  - `temporaryCount`: number of temporary plots (`n - m`).
  - `N`: number of possible plots in the population.

# Mathematical basis

The regression slope `b` is solved via the normal equations `(X'X)β = X'y` over the
permanent plots (see `_olsslope`), the same transposed-design-matrix approach used for
the ANOVA in [`stratifiedSampling`](@ref):
```math
\\hat{y}_{reg} = \\bar{y}_m + b(\\bar{x}_n - \\bar{x}_m)
```
where `x̄ₙ` is the occasion-1 mean over all `n` plots and `x̄ₘ`/`ȳₘ` are the occasion-1/2
means over only the `m` permanent plots. Its variance combines the regression's residual
variance `S²yx` (from the `m` permanent plots) with the extra precision gained from the
larger first-phase sample:
```math
s^2_{\\hat{y}_{reg}} = \\frac{S^2_{yx}}{m} + \\frac{S^2_y - S^2_{yx}}{n}
```

# Examples

```julia-repl
julia> value1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0, 23.1, 17.6, 20.8, 21.9]u"m^3";
julia> value2 = [22.1, 25.8, 23.4, missing, 26.6, missing, 27.5, missing, missing, 26.0]u"m^3";

julia> report = doubleSampling(value1, value2, 0.05, 200);
julia> occasion2(report).mean
24.441857923497274 m^3
```
"""
function doubleSampling(value1::AbstractVector, value2::AbstractVector, plot_area::Area, N::Integer; α::Real=0.95)

  length(value1) == length(value2) || throw(DimensionMismatch("value1 and value2 must have the same length."))
  permanentidx = findall(!ismissing, value2)
  m = length(permanentidx)
  n = length(value1)
  m >= 3 || throw(ArgumentError("At least three permanent (remeasured) plots are required."))

  o1 = _occasionmean(value1, N)
  occasion1 = _occasiontable(o1, N, plot_area, α)
  vunit = unit(o1.x̅)

  x1perm = ustrip.(value1[permanentidx])
  y2perm = ustrip.(value2[permanentidx])
  x̄1perm, ȳ2perm = mean(x1perm), mean(y2perm)
  x̄1all = ustrip(o1.x̅)

  b = _olsslope(y2perm, x1perm)
  sstY = sum(abs2, y2perm .- ȳ2perm)
  ssxy = sum((x1perm .- x̄1perm) .* (y2perm .- ȳ2perm))
  syx² = (sstY - b * ssxy) / (m - 2)
  sy² = sstY / (m - 1)

  yreg = ȳ2perm + b * (x̄1all - x̄1perm)
  meanvar = syx² / m + (sy² - syx²) / n
  yregq, sx̅ = yreg * vunit, sqrt(meanvar) * vunit

  t = -quantile(TDist(m - 1), (1 - α) / 2)
  absoluteerror = t * sx̅
  relativeerror = absoluteerror / yregq * 100
  hectarevalue = yregq / plot_area
  totalvalue = yregq * N

  occasion2 = DataFrame(
    mean=yregq, regressionSlope=round(b, digits=4), residualVariance=syx² * vunit^2, standardError=sx̅,
    absoluteError=absoluteerror, relativeError=round(relativeerror, digits=2), meanPerHectare=hectarevalue,
    total=totalvalue, lowerLimit=totalvalue - N * absoluteerror, upperLimit=totalvalue + N * absoluteerror,
    m=m, temporaryCount=n - m, N=N,
  )

  growth = yregq - o1.x̅
  change = _changetable(growth, meanvar, m - 1, N, α)

  return SamplingReport((; occasion1, occasion2, change))
end

function doubleSampling(value1::AbstractVector, value2::AbstractVector, plot_area::Real, N::Integer; kwargs...)
  doubleSampling(value1, value2, plot_area * AUNIT, N; kwargs...)
end
