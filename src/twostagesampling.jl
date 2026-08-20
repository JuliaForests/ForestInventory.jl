# Two-stage sampling's variance decomposition, unlike one-stage cluster sampling's
# `_clustervariance`, has two independent sources of sampling error: which `n` (of `N`)
# primary units were drawn, and which `m` (of `M`) secondary units were measured within
# each of them. Cochran (1977), eq. 10.4, equal-size primaries, SRS at both stages.
function _twostagevariance(means::AbstractVector{<:Real}, vars::AbstractVector{<:Real}, m::Real, M::Real, N::Real)
  n = length(means)
  ȳ = mean(means)
  withinvar = mean(vars)
  msb = m * sum(abs2, means .- ȳ) / (n - 1)
  betweenvar = (msb - withinvar) / m
  f1, f2 = n / N, m / M
  meanvar = (1 - f1) / n * betweenvar + f1 * (1 - f2) / (n * m) * withinvar
  varterm = betweenvar + withinvar / m
  correctionvarterm = betweenvar + withinvar / M
  return (; withinvar, betweenvar, meanvar, varterm, correctionvarterm)
end

"""
    twoStageSampling(primary::AbstractVector, value::AbstractVector, plot_area::Area, N::Integer, M::Integer;
                      e::Real=10, α::Real=0.95)
    twoStageSampling(primary::AbstractVector, value::AbstractVector, plot_area::Real, N::Integer, M::Integer; kwargs...)

Performs two-stage sampling for forest inventory analysis, generic to any number of primary units.

# Description

Two-stage sampling draws `n` primary units from a population of `N` (e.g. stands or
blocks), then sub-samples `m` secondary units (plots) from within each of the `M` possible
secondary units of every drawn primary — unlike [`clusterSampling`](@ref), which
measures every secondary unit inside a drawn cluster. Sub-sampling trades some precision
for a further reduction in fieldwork, since not every drawn primary unit needs a full
census.

# Arguments

- `primary::AbstractVector`: each plot's primary-unit label, same order as `value`.
- `value::AbstractVector`: measurement for each plot, same order as `primary`. If given as a `Unitful` quantity, its unit is preserved; plain numbers are kept unitless throughout.
- `plot_area`: the area of each plot (secondary unit), as an `Area` quantity, used only to report the per-hectare value. A plain number is taken to be hectares.
- `N::Integer`: total number of primary units in the population.
- `M::Integer`: total number of secondary units per primary unit in the population (assumed equal across primary units).
- `e::Real=10`: desired relative error margin as a percentage (default 10%).
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- [`SamplingReport`](@ref) with `primaryTable` (per-primary-unit descriptive statistics)
  and `resultTable` (one row, one column per statistic). `resultTable` columns:
  - `mean`, `coefficientOfVariation`, `standardError`, `absoluteError`, `relativeError`, `meanPerHectare`, `total`, `lowerLimit`, `upperLimit`, `population`, `samplingFraction`: as in [`simpleCasualSampling`](@ref), using the two-stage mean/variance.
  - `withinVariance`, `betweenVariance`: within-primary and between-primary variance.
  - `m`: number of secondary units measured per primary unit.
  - `M`: number of secondary units per primary unit in the population.
  - `n`, `requiredSampleSize`, `additionalSampleSize`: measured, required, and missing number of primary units.
  - `N`: number of possible primary units in the population.

# Mathematical basis

With sampling fractions `f₁ = n/N` (primary stage) and `f₂ = m/M` (secondary stage):
```math
s^2_{\\bar{x}} = \\frac{1-f_1}{n}s^2_b + \\frac{f_1(1-f_2)}{nm}s^2_w
```
Setting `m = M` (every drawn primary fully censused) removes the secondary-stage
contribution entirely and this reduces to [`clusterSampling`](@ref)'s formula.

# Technical description

The required number of primary units is solved the same way as every other design in
this module (see `_requiredsamplesize`), except the finite-population correction term
uses the *population* secondary-unit count `M` while the point estimate itself uses the
*sampled* count `m` — the only design in this module where those two differ.

# Examples

```julia-repl
julia> primary = repeat(1:5, inner=3);

julia> value = [18.2, 19.1, 17.8, 22.4, 23.1, 21.9, 15.1, 16.0, 14.8, 27.3, 28.1, 26.9, 19.8, 20.5, 19.1]u"m^3";

julia> report = twoStageSampling(primary, value, 0.02, 40, 6);
julia> resultTable(report).mean
20.673333333333336 m^3
```
"""
function twoStageSampling(primary::AbstractVector, value::AbstractVector, plot_area::Area, N::Integer, M::Integer;
  e::Real=10, α::Real=0.95)

  length(primary) == length(value) || throw(DimensionMismatch("primary and value must have the same length."))

  table = _clustertable(primary, value; groupname=:primary)
  m = _equalclustersize(table.n)
  n = nrow(table)
  m <= M || throw(ArgumentError("The number of measured secondary units per primary unit (m=$m) cannot exceed the population size M=$M."))

  x̅ = mean(table.mean)
  vunit = unit(x̅)
  v2unit = vunit^2

  variance = _twostagevariance(ustrip.(table.mean), ustrip.(table.variance), m, M, N)
  cv = sqrt(variance.withinvar + variance.betweenvar) / ustrip(x̅) * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  infinite = _isinfinitepopulation(n, N)
  E = (e / 100) * x̅
  requiredprimaries = _requiredsamplesize(variance.varterm, n, ustrip(E), α, infinite ? nothing : N;
    correctionvarterm=variance.correctionvarterm)
  missingprimaries = n > requiredprimaries ? 0 : requiredprimaries - n

  population = infinite ? "infinite" : "finite"
  sx̅ = sqrt(variance.meanvar) * vunit
  absoluteerror = Ttab * sx̅
  relativeerror = absoluteerror / x̅ * 100
  hectarevalue = x̅ / uconvert(referencearea(plot_area), plot_area)
  totalvalue = x̅ * N * M
  ciupper = totalvalue + N * M * absoluteerror
  cilower = totalvalue - N * M * absoluteerror

  resulttable = DataFrame(
    mean=x̅, coefficientOfVariation=round(cv, digits=2), withinVariance=variance.withinvar * v2unit,
    betweenVariance=variance.betweenvar * v2unit, standardError=sx̅, absoluteError=absoluteerror,
    relativeError=round(relativeerror, digits=2), meanPerHectare=hectarevalue, total=totalvalue,
    lowerLimit=cilower, upperLimit=ciupper, population=population, samplingFraction=round(1 - n / N, digits=3),
    m=m, M=M, n=n, requiredSampleSize=requiredprimaries, additionalSampleSize=missingprimaries, N=N,
  )

  return SamplingReport((; primaryTable=table, resultTable=resulttable))
end

function twoStageSampling(primary::AbstractVector, value::AbstractVector, plot_area::Real, N::Integer, M::Integer; kwargs...)
  twoStageSampling(primary, value, plot_area * AUNIT, N, M; kwargs...)
end
