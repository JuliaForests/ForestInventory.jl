"""
    multistartSystematicSampling(start::AbstractVector, value::AbstractVector, plot_area::Area, total_area::Area;
                                  e::Real=10, α::Real=0.95)
    multistartSystematicSampling(start::AbstractVector, value::AbstractVector, plot_area::Real, total_area::Real; kwargs...)

Performs systematic sampling with multiple random starts for forest inventory analysis, generic to any number of starts.

# Description

A single systematic line (see [`systematicSampling`](@ref)) gives no internal replication
to estimate its own precision from first principles — the successive-differences method
is an approximation. Laying out several independent systematic lines instead, each at its
own random starting point along the sampling interval, turns the design into a genuine
probability sample: every line is itself a systematic subsample of `M` plots, and the `n`
lines play the role of primary sampling units drawn from the `N` possible starting points.

# Arguments

- `start::AbstractVector`: each plot's random-start (systematic line) label, same order as `value`.
- `value::AbstractVector`: measurement for each plot, same order as `start`. If given as a `Unitful` quantity, its unit is preserved; plain numbers are kept unitless throughout.
- `plot_area`: the area of each individual plot, as an `Area` quantity. A plain number is taken to be hectares.
- `total_area`: the total area of the forest or stand, as an `Area` quantity. A plain number is taken to be hectares.
- `e::Real=10`: desired relative error margin as a percentage (default 10%).
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- [`SamplingReport`](@ref) with `startTable` (per-start descriptive statistics) and
  `resultTable` (one row, one column per statistic — same shape and column names as
  [`clusterSampling`](@ref), since both share the same estimator).

# Mathematical basis

Statistically identical to [`clusterSampling`](@ref)'s between/within decomposition —
each start is a primary unit of `M` plots — with the population size `N` given the
systematic-sampling interpretation of "number of possible random starts" rather than
"number of possible spatial clusters": `N = total_area / (plot_area × M)`, the number of
non-overlapping systematic subsamples that tile the population.

# Technical description

The only reason this is a separate design from [`clusterSampling`](@ref) rather than a
mode flag is field procedure, not statistics: a "start" is a systematically-spread
subsample covering the whole area, while a "cluster" is a spatially contiguous group of
plots. Their estimators share the same variance decomposition, so both call the same
internal `_clustervariance` engine.

# Examples

```julia-repl
julia> start = repeat(1:4, inner=5);

julia> value = [18.2, 19.1, 17.8, 18.9, 19.4, 22.4, 23.1, 21.9, 22.8, 23.4,
                 15.1, 16.0, 14.8, 15.6, 15.9, 20.5, 21.2, 19.8, 20.6, 21.0]u"m^3";

julia> report = multistartSystematicSampling(start, value, 0.02, 15);
julia> resultTable(report).mean
19.375 m^3
```
"""
function multistartSystematicSampling(start::AbstractVector, value::AbstractVector, plot_area::Area, total_area::Area;
  e::Real=10, α::Real=0.95)

  length(start) == length(value) || throw(DimensionMismatch("start and value must have the same length."))

  table = _clustertable(start, value; groupname=:start)
  M = _equalclustersize(table.n)
  n = nrow(table)

  N = round(Int, ustrip(uconvert(NoUnits, total_area / (plot_area * M))))
  x̅ = mean(table.mean)
  vunit = unit(x̅)
  v2unit = vunit^2

  variance = _clustervariance(ustrip.(table.mean), ustrip.(table.variance), M, N)
  cv = sqrt(variance.totalvar) / ustrip(x̅) * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  infinite = _isinfinitepopulation(n, N)
  E = (e / 100) * x̅
  requiredplots = _requiredsamplesize(variance.varterm, n, ustrip(E), α, infinite ? nothing : N)
  missingplots = n > requiredplots ? 0 : requiredplots - n

  population = infinite ? "infinite" : "finite"
  sx̅ = sqrt(variance.meanvar) * vunit
  absoluteerror = Ttab * sx̅
  relativeerror = absoluteerror / x̅ * 100
  hectarevalue = x̅ / uconvert(referencearea(plot_area), plot_area)
  totalvalue = x̅ * M * N
  ciupper = totalvalue + M * N * absoluteerror
  cilower = totalvalue - M * N * absoluteerror

  resulttable = DataFrame(
    mean=x̅, coefficientOfVariation=round(cv, digits=2), withinVariance=variance.withinvar * v2unit,
    betweenVariance=variance.betweenvar * v2unit, totalVariance=variance.totalvar * v2unit,
    intraclassCorrelation=round(variance.icc, digits=4), standardError=sx̅, absoluteError=absoluteerror,
    relativeError=round(relativeerror, digits=2), meanPerHectare=hectarevalue, total=totalvalue,
    lowerLimit=cilower, upperLimit=ciupper, population=population, samplingFraction=round(1 - n / N, digits=3),
    M=M, n=n, requiredSampleSize=requiredplots, additionalSampleSize=missingplots, N=N, area=total_area,
  )

  return SamplingReport((; startTable=table, resultTable=resulttable))
end

function multistartSystematicSampling(start::AbstractVector, value::AbstractVector, plot_area::Real, total_area::Real; kwargs...)
  multistartSystematicSampling(start, value, plot_area * AUNIT, total_area * AUNIT; kwargs...)
end
