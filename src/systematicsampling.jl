# Successive differences between consecutive plots along a systematic line, skipping the
# pair that straddles two different lines when a `line` grouping is given -- generic to
# any number of lines (transects), reducing to a single running difference when there is
# only one.
function _successivedifferences(value::AbstractVector, line::Union{AbstractVector,Nothing})
  n = length(value)
  keep = isnothing(line) ? trues(n - 1) : [line[i] == line[i+1] for i in 1:n-1]
  return [value[i] - value[i+1] for i in 1:n-1 if keep[i]]
end

"""
    systematicSampling(value::AbstractVector, plot_area::Area, total_area::Area;
                        line::Union{AbstractVector,Nothing}=nothing, e::Real=10, α::Real=0.95)
    systematicSampling(value::AbstractVector, plot_area::Real, total_area::Real; kwargs...)

Performs systematic sampling for forest inventory analysis using the method of successive differences.

# Description

In systematic sampling, plots are laid out at a fixed interval along one or more lines
(transects) instead of being chosen at random, which is often far cheaper to execute in
the field. Because consecutive plots tend to be more alike than two random plots, its
precision cannot be estimated with the simple-random-sampling variance formula; the
method of successive differences below is the standard substitute (Cochran, 1977, §8.6).

# Arguments

- `value`: vector of plot measurements, **in the order the plots were measured along the line(s)**. If given as a `Unitful` quantity, its unit is preserved; plain numbers are kept unitless throughout.
- `plot_area`: the area of each sample plot, as an `Area` quantity. A plain number is taken to be hectares.
- `total_area`: the total area of the forest or stand, as an `Area` quantity. A plain number is taken to be hectares.
- `line::Union{AbstractVector,Nothing}=nothing`: an optional grouping vector, one entry per plot, identifying which line/transect each plot belongs to (same order as `value`). When given, the difference between the last plot of one line and the first plot of the next is excluded, since they are not actually adjacent on the ground. `nothing` (the default) treats every plot as a single line.
- `e::Real=10`: desired relative error margin as a percentage, only used to report whether it was met (default 10%).
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- `DataFrame`: a single-row report, one column per statistic (same shape as [`simpleCasualSampling`](@ref)). Columns:
  - `mean`, `coefficientOfVariation`, `varianceOfMean`, `standardError`, `absoluteError`, `relativeError`, `meanPerHectare`, `total`, `lowerLimit`, `upperLimit`, `population`, `samplingFraction`: as in [`simpleCasualSampling`](@ref).
  - `admissibleError`: the target relative error `e`, echoed back for comparison against `relativeError` (this design reports no required sample size — see Technical description below).
  - `samplingInterval`: number of lines (transects), `k`.
  - `n`: number of measured plots.
  - `N`: number of possible plots in the population.
  - `area`: total area.

# Mathematical basis

For `k` lines totalling `n` plots, the variance of the mean is estimated from the sum of
squared successive differences within each line:
```math
s^2_{\\bar{x}} = \\frac{\\sum (x_i - x_{i+1})^2}{2n(n-k)}\\left(1 - \\frac{n}{N}\\right)
```
which reduces to the classic single-line formula `Σ(xᵢ-xᵢ₊₁)² / (2n(n-1))` when `k = 1`.

# Technical description

Systematic sampling has no general closed-form sample-size formula the way simple or
stratified sampling do — its precision depends on how the underlying variable happens to
be arranged along the line, not just on `n` — so, matching the source this design was
ported from, this function reports whether the achieved relative error (`relativeError`)
meets the target (`admissibleError`) instead of a required plot count.

# Examples

```julia-repl
julia> v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1]u"m^3";

julia> report = systematicSampling(v, 0.05, 10; e=10, α=0.95)
julia> report.mean
441.6909090909091 m^3
julia> report.relativeError <= report.admissibleError
true
```
"""
function systematicSampling(value::AbstractVector, plot_area::Area, total_area::Area;
  line::Union{AbstractVector,Nothing}=nothing, e::Real=10, α::Real=0.95)

  n = length(value)
  n > 1 || throw(ArgumentError("At least two sampled plots are required."))
  isnothing(line) || length(line) == n || throw(DimensionMismatch("line must have the same length as value."))

  N = round(Int, ustrip(uconvert(NoUnits, total_area / plot_area)))
  k = isnothing(line) ? 1 : length(unique(line))
  x̅ = mean(value)
  diffs = _successivedifferences(value, line)

  infinite = _isinfinitepopulation(n, N)
  f = infinite ? 1.0 : 1 - n / N
  population = infinite ? "infinite" : "finite"
  s²x̅ = sum(abs2, diffs) / (2n * (n - k)) * f
  sx̅ = sqrt(s²x̅)
  cv = sx̅ / x̅ * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  absoluteerror = Ttab * sx̅
  relativeerror = absoluteerror / x̅ * 100
  hectarevalue = x̅ / uconvert(referencearea(plot_area), plot_area)
  totalvalue = x̅ * N
  ciupper = totalvalue + N * absoluteerror
  cilower = totalvalue - N * absoluteerror

  return DataFrame(
    mean=x̅, coefficientOfVariation=round(cv, digits=2), varianceOfMean=s²x̅, standardError=sx̅,
    absoluteError=absoluteerror, relativeError=round(relativeerror, digits=2), admissibleError=e,
    meanPerHectare=hectarevalue, total=totalvalue, lowerLimit=cilower, upperLimit=ciupper,
    population=population, samplingFraction=round(f, digits=3), samplingInterval=k, n=n, N=N, area=total_area,
  )
end

systematicSampling(value::AbstractVector, plot_area::Real, total_area::Real; kwargs...) =
  systematicSampling(value, plot_area * AUNIT, total_area * AUNIT; kwargs...)
