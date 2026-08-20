"""
    independentOccasionsSampling(value1::AbstractVector, value2::AbstractVector, plot_area::Area,
                                  N1::Integer, N2::Integer; α::Real=0.95)
    independentOccasionsSampling(value1::AbstractVector, value2::AbstractVector, plot_area::Real,
                                  N1::Integer, N2::Integer; kwargs...)

Estimates growth between two inventory occasions using two independent (unmatched) samples, generic to any sample sizes.

# Description

The simplest of the successive-occasions designs: occasion 1 and occasion 2 are sampled
independently, with no attempt to remeasure the same plots. This is easy to execute but
the least efficient at detecting growth, since none of the natural plot-to-plot
correlation between the two occasions is exploited — see [`completeReplacementSampling`](@ref)
and [`partialReplacementSampling`](@ref) for designs that remeasure some or all plots.

# Arguments

- `value1`, `value2`: independent vectors of plot measurements from occasion 1 and occasion 2 respectively (need not be the same length). If given as `Unitful` quantities, their unit is preserved; plain numbers are kept unitless throughout.
- `plot_area`: the area of each plot, as an `Area` quantity, used only to report the per-hectare value. A plain number is taken to be hectares.
- `N1`, `N2`: total number of possible plots in the population at occasion 1 and occasion 2 respectively.
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- [`SamplingReport`](@ref) with `occasion1`, `occasion2` and `change`, each a single-row,
  one-column-per-statistic `DataFrame`:
  - `occasion1`/`occasion2` columns: `mean`, `standardError`, `absoluteError`, `relativeError`, `meanPerHectare`, `total`, `lowerLimit`, `upperLimit`, `n`, `N`.
  - `change` columns: `meanGrowth`, `standardError`, `absoluteError`, `relativeError`, `totalGrowth`, `lowerLimit`, `upperLimit`.

# Mathematical basis

Since the two samples are independent, the variance of the estimated growth is simply the
sum of each occasion's variance of the mean, with no covariance term:
```math
\\bar{g} = \\bar{x}_2 - \\bar{x}_1 \\qquad s^2_{\\bar{g}} = s^2_{\\bar{x}_1} + s^2_{\\bar{x}_2}
```
using `df = (n_1-1)+(n_2-1)-1` degrees of freedom for the growth confidence interval.

# Examples

```julia-repl
julia> v1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0]u"m^3";
julia> v2 = [24.1, 23.8, 22.9, 25.6, 24.0]u"m^3";

julia> report = independentOccasionsSampling(v1, v2, 0.05, 200, 200);
julia> change(report).meanGrowth
3.913333333333334 m^3
```
"""
function independentOccasionsSampling(value1::AbstractVector, value2::AbstractVector, plot_area::Area,
  N1::Integer, N2::Integer; α::Real=0.95)

  o1 = _occasionmean(value1, N1)
  o2 = _occasionmean(value2, N2)
  occasion1 = _occasiontable(o1, N1, plot_area, α)
  occasion2 = _occasiontable(o2, N2, plot_area, α)

  growth = o2.x̅ - o1.x̅
  growthvar = ustrip(o1.meanvar) + ustrip(o2.meanvar)
  df = (o1.n - 1) + (o2.n - 1) - 1
  change = _changetable(growth, growthvar, df, N2, α)

  return SamplingReport((; occasion1, occasion2, change))
end

function independentOccasionsSampling(value1::AbstractVector, value2::AbstractVector, plot_area::Real,
  N1::Integer, N2::Integer; kwargs...)
  independentOccasionsSampling(value1, value2, plot_area * AUNIT, N1, N2; kwargs...)
end
