"""
    completeReplacementSampling(value1::AbstractVector, value2::AbstractVector, plot_area::Area,
                                 N1::Integer, N2::Integer; α::Real=0.95)
    completeReplacementSampling(value1::AbstractVector, value2::AbstractVector, plot_area::Real,
                                 N1::Integer, N2::Integer; kwargs...)

Estimates growth between two inventory occasions by remeasuring the exact same plots on both occasions, generic to any sample size.

# Description

The same `n` permanent plots are measured at occasion 1 and remeasured at occasion 2
("complete/total replacement" in the sampling-on-successive-occasions terminology). Since
`value1[i]` and `value2[i]` refer to the same plot, the two occasions are correlated
rather than independent — usually positively, since a plot with above-average value at
occasion 1 tends to remain above average at occasion 2 — and exploiting that correlation
tightens the growth estimate compared to [`independentOccasionsSampling`](@ref).

# Arguments

- `value1`, `value2`: matched vectors of plot measurements at occasion 1 and occasion 2 (`value1[i]`/`value2[i]` are the same plot), of equal length. If given as `Unitful` quantities, their unit is preserved; plain numbers are kept unitless throughout.
- `plot_area`: the area of each plot, as an `Area` quantity, used only to report the per-hectare value. A plain number is taken to be hectares.
- `N1`, `N2`: total number of possible plots in the population at occasion 1 and occasion 2 respectively.
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- [`SamplingReport`](@ref) with `occasion1`, `occasion2` and `change` — same column
  layout as [`independentOccasionsSampling`](@ref).

# Mathematical basis

The matched design subtracts twice the sample covariance between the two occasions from
the variance of the growth estimate:
```math
\\bar{g} = \\bar{x}_2 - \\bar{x}_1 \\qquad
s^2_{\\bar{g}} = s^2_{\\bar{x}_1} + s^2_{\\bar{x}_2} - \\frac{2\\,\\mathrm{Cov}(x_1,x_2)}{n}
```
using `df = 2(n-1)-1` degrees of freedom for the growth confidence interval. Positive
correlation between occasions makes `s²_ḡ` smaller than the independent-samples case.

# Examples

```julia-repl
julia> v1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0, 24.6, 17.3]u"m^3";
julia> v2 = [22.1, 25.8, 23.4, 21.0, 26.6, 20.9, 26.0, 22.5]u"m^3";

julia> report = completeReplacementSampling(v1, v2, 0.05, 200, 200);
julia> change(report).meanGrowth
3.1750000000000007 m^3
```
"""
function completeReplacementSampling(value1::AbstractVector, value2::AbstractVector, plot_area::Area,
  N1::Integer, N2::Integer; α::Real=0.95)

  length(value1) == length(value2) ||
    throw(DimensionMismatch("value1 and value2 must have the same length: the same plots are remeasured on both occasions."))

  o1 = _occasionmean(value1, N1)
  o2 = _occasionmean(value2, N2)
  occasion1 = _occasiontable(o1, N1, plot_area, α)
  occasion2 = _occasiontable(o2, N2, plot_area, α)

  n = o1.n
  growth = o2.x̅ - o1.x̅
  growthcovariance = ustrip(cov(value1, value2))
  # with a small sample and strong occasion-to-occasion correlation, the covariance
  # estimate can occasionally overshoot enough to make this formula dip below zero;
  # clamped at zero, since a variance estimate cannot legitimately be negative.
  growthvar = max(ustrip(o1.meanvar) + ustrip(o2.meanvar) - 2 * growthcovariance / n, 0.0)
  df = 2 * (n - 1) - 1
  change = _changetable(growth, growthvar, df, N2, α)

  return SamplingReport((; occasion1, occasion2, change))
end

function completeReplacementSampling(value1::AbstractVector, value2::AbstractVector, plot_area::Real,
  N1::Integer, N2::Integer; kwargs...)
  completeReplacementSampling(value1, value2, plot_area * AUNIT, N1, N2; kwargs...)
end
