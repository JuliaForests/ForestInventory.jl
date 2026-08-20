# Per-primary-unit sample size, mean and internal variance, one row per primary unit
# (cluster/start/primary), generic to any number of them via `groupby`. `groupname` is
# the column name to expose in the returned table (`:cluster`, `:start`, `:primary`) --
# shared by `clusterSampling`, `multistartSystematicSampling`, and `twoStageSampling`,
# which differ only in what the primary unit is called and how the population size is
# derived from it, not in this statistics step.
function _clustertable(group::AbstractVector, value::AbstractVector; groupname::Symbol=:cluster)
  data = DataFrame(; groupname => group, value=value)
  return combine(groupby(data, groupname, sort=true)) do df
    (n=nrow(df), mean=mean(df.value), variance=var(df.value))
  end
end

# Validates that every primary unit has the same number of secondary units and returns
# it; `_clustervariance`'s decomposition assumes equal-size primary units.
function _equalclustersize(n::AbstractVector{<:Integer})
  sizes = unique(n)
  length(sizes) == 1 ||
    throw(ArgumentError("Every cluster must contain the same number of secondary units; found sizes $(sort(sizes))."))
  return only(sizes)
end

"""
    clusterSampling(cluster::AbstractVector, value::AbstractVector, plot_area::Area, total_area::Area;
                     e::Real=10, α::Real=0.95)
    clusterSampling(cluster::AbstractVector, value::AbstractVector, plot_area::Real, total_area::Real; kwargs...)

Performs one-stage cluster sampling for forest inventory analysis, generic to any number of clusters.

# Description

In cluster sampling the sampling unit is a cluster of `M` neighboring plots (a
"conglomerate") rather than an individual plot — cheaper to lay out in the field, at the
cost of the within-cluster homogeneity typically inflating the variance of the mean
relative to an equivalent simple random sample. Every cluster is assumed to contain the
same number of plots `M` (Cochran, 1977, §9.3, "clusters of equal size").

# Arguments

- `cluster::AbstractVector`: each plot's cluster label, same order as `value`.
- `value::AbstractVector`: measurement for each plot, same order as `cluster`. If given as a `Unitful` quantity, its unit is preserved; plain numbers are kept unitless throughout.
- `plot_area`: the area of each individual plot (secondary unit), as an `Area` quantity. A plain number is taken to be hectares.
- `total_area`: the total area of the forest or stand, as an `Area` quantity. A plain number is taken to be hectares.
- `e::Real=10`: desired relative error margin as a percentage (default 10%).
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- [`SamplingReport`](@ref) with `clusterTable` (per-cluster descriptive statistics) and
  `resultTable` (one row, one column per statistic). `resultTable` columns:
  - `mean`, `coefficientOfVariation`, `standardError`, `absoluteError`, `relativeError`, `meanPerHectare`, `total`, `lowerLimit`, `upperLimit`, `population`, `samplingFraction`: as in [`simpleCasualSampling`](@ref), using the cluster-sampling mean/variance.
  - `withinVariance`, `betweenVariance`, `totalVariance`: within-cluster, between-cluster, and total variance (per plot).
  - `intraclassCorrelation`: intraclass correlation coefficient `ρ`.
  - `M`: number of secondary units (plots) per cluster.
  - `n`, `requiredSampleSize`, `additionalSampleSize`: measured, required, and missing number of clusters.
  - `N`: number of possible clusters in the population.
  - `area`: total area.

# Mathematical basis

With `n` sampled clusters of `M` plots each out of `N` possible clusters in the
population, the within- and between-cluster variance (per plot) are:
```math
s^2_w = \\overline{s^2_j} \\qquad
s^2_b = \\frac{M\\sum_j (\\bar{x}_j-\\bar{x})^2/(n-1) - s^2_w}{M}
```
and the variance of the overall mean, with finite-population correction `f = 1-n/N`:
```math
s^2_{\\bar{x}} = f\\frac{s^2_b}{n} + \\frac{s^2_w}{nM}
```
The intraclass correlation coefficient `ρ = s²_b/(s²_b+s²_w)` measures how similar plots
within the same cluster are: `ρ = 0` recovers the simple-random-sampling variance exactly.

# Technical description

`clusterSampling` with `M = 1` (one plot per cluster) is numerically identical to
[`simpleCasualSampling`](@ref) — clustering stops mattering once there is nothing left
to cluster. The population size `N` is the number of possible *clusters*, not plots:
`total_area / (plot_area × M)`.

# Examples

```julia-repl
julia> cluster = repeat(1:6, inner=4);

julia> value = [18.2, 19.1, 17.8, 18.9, 22.4, 23.1, 21.9, 22.8, 15.1, 16.0, 14.8, 15.6,
                 27.3, 28.1, 26.9, 27.8, 19.8, 20.5, 19.1, 20.0, 24.5, 25.2, 23.9, 24.8]u"m^3";

julia> report = clusterSampling(cluster, value, 0.02, 15);
julia> resultTable(report).mean
21.4 m^3
julia> clusterTable(report)
```
"""
function clusterSampling(cluster::AbstractVector, value::AbstractVector, plot_area::Area, total_area::Area;
  e::Real=10, α::Real=0.95)

  length(cluster) == length(value) || throw(DimensionMismatch("cluster and value must have the same length."))

  table = _clustertable(cluster, value; groupname=:cluster)
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

  return SamplingReport((; clusterTable=table, resultTable=resulttable))
end

clusterSampling(cluster::AbstractVector, value::AbstractVector, plot_area::Real, total_area::Real; kwargs...) =
  clusterSampling(cluster, value, plot_area * AUNIT, total_area * AUNIT; kwargs...)
