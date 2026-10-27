"""
    completereplacementsampling(value1::AbstractVector, value2::AbstractVector, plotarea::Area,
                                 npop1::Integer, npop2::Integer; α::Real=0.95)

Estimates the growth (change) of a forest between two inventory occasions by remeasuring the same plots, generic to any sample size.

# Description

The same ``n`` permanent plots are measured on occasion 1 and remeasured on occasion 2 (sampling with complete replacement, or total repetition). Since `value1[i]` and `value2[i]` are the same plot, the two occasions are correlated, usually positively, and subtracting that correlation tightens the estimate of the growth with respect to [`independentoccasionssampling`](@ref). The notation and the formulas follow Péllico Netto & Brena (1997), chapter 11, Sec. 11.2.

# Arguments

- `value1`, `value2`: the matched vectors of plot measurements on occasion 1 and occasion 2 (`value1[i]` and `value2[i]` are the same plot), of equal length. A `Unitful` quantity keeps its unit; a plain number stays unitless.
- `plotarea`: the area of each plot, an `Area`, used only for the mean per hectare.
- `npop1`, `npop2`: the possible plots of the population on occasion 1 and on occasion 2.
- `α::Real=0.95`: the confidence level.

# Returns

A [`SamplingReport`](@ref) with `occasion1`, `occasion2` and `change`, with the same columns as [`independentoccasionssampling`](@ref):

- `occasion1`, `occasion2`: one row each, with `xm` (mean per plot of the occasion, ``\\bar{x}`` and ``\\bar{y}``), `sxm` (standard error), `Ea`, `Er`, `xha` (mean per hectare), `X` (total, ``N \\bar{x}``), `Xl`, `Xu` (limits of its confidence interval), `n` (plots measured) and `N` (possible plots).
- `change`: one row with `dm` (mean change per plot, ``\\bar{d}``), `sdm` (its standard error), `Ea`, `Er`, `D` (total change, ``\\hat{D} = N \\bar{d}``, with the possible plots of the second occasion) and `Dl`, `Du` (limits of its confidence interval).

# Mathematical basis

On each occasion the mean, the variance of the mean with the correction of the finite population, the standard error and the sampling errors are those of [`simplecasualsampling`](@ref) (Eqs. 11.1 to 11.9), with ``n`` plots of ``N`` possible and Student's ``t`` of ``n - 1`` degrees of freedom:

```math
\\bar{x} = \\frac{\\sum_{i=1}^{n} x_i}{n} \\qquad
s_{\\bar{x}}^2 = \\frac{s_x^2}{n} \\left(1 - \\frac{n}{N}\\right) \\qquad
E_a = t \\, s_{\\bar{x}} \\qquad E_r = \\frac{t \\, s_{\\bar{x}}}{\\bar{x}} \\times 100 \\qquad
\\hat{X} = N \\bar{x}
```

The change is the difference between the means of the ``m = n`` remeasured plots (Eq. 11.33). The variance of the change subtracts the covariance between the occasions (Eqs. 11.34 and 11.35), which is why the precision is larger than that of independent samples:

```math
\\bar{d} = \\bar{y} - \\bar{x} \\qquad
s_{\\bar{d}}^2 = s_{\\bar{x}}^2 + s_{\\bar{y}}^2 - \\frac{2 s_{xy}}{m} \\qquad
s_{xy} = \\frac{\\sum_{j=1}^{m} (x_j - \\bar{x})(y_j - \\bar{y})}{m - 1}
```

With few plots and a strong correlation the estimate of the covariance may overshoot and make ``s_{\\bar{d}}^2`` negative, and it is set to zero in that case, as a variance cannot be negative. The errors and the total change are ``E_a = t \\, s_{\\bar{d}}``, ``E_r = t \\, s_{\\bar{d}} / \\bar{d} \\times 100`` and ``\\hat{D} = N_2 \\bar{d}`` (Eqs. 11.37 to 11.41). The degrees of freedom of ``t`` are ``(n_1 - 1) + (n_2 - 1) = 2(n - 1)`` (text of Eq. 11.39).

# Examples

The example of the book (chapter 11, Quadro 11.2): *Eucalyptus* sp. in Lençóis Paulista, SP, with 45 permanent plots measured in 1976 and 1977, in a population of 1 500 units. The volumes of the book are per hectare, so they are the volumes of units of 1 ha.

```julia-repl
julia> v1 = [141.89, 147.74, 126.26, 137.08, 39.30, 133.19, 110.21, 126.79, 132.09, 73.79, 76.22, 63.42,
        77.47, 37.24, 98.34, 87.04, 116.52, 83.75, 110.99, 77.87, 139.73, 182.11, 94.68, 181.69,
        161.44, 172.58, 191.29, 137.93, 167.08, 110.10, 123.59, 107.45, 137.35, 137.77, 79.56,
        80.79, 118.52, 133.61, 174.91, 140.90, 176.79, 185.32, 135.84, 109.97, 73.65]u"m^3";
julia> v2 = [174.99, 187.58, 153.81, 176.28, 61.76, 160.87, 150.23, 164.15, 170.99, 101.36, 114.85, 87.89,
        102.06, 51.37, 152.57, 119.94, 148.93, 129.42, 152.01, 109.04, 165.91, 208.33, 143.74,
        197.50, 191.62, 266.13, 280.29, 178.37, 218.88, 149.68, 161.67, 174.45, 172.66, 179.97,
        105.36, 126.26, 150.55, 191.47, 248.38, 183.30, 226.76, 224.25, 154.10, 142.97, 102.68]u"m^3";

julia> report = completereplacementsampling(v1, v2, 1u"ha", 1500, 1500);

julia> occasion1(report).xm[1]
121.15222222222224 m^3

julia> change(report)
1×7 DataFrame
 Row │ dm           sdm          Ea           Er       D            Dl           Du
     │ Quantity…    Quantity…    Quantity…    Float64  Quantity…    Quantity…    Quantity…
─────┼───────────────────────────────────────────────────────────────────────────────────────
   1 │ 39.1896 m^3  1.88369 m^3  3.74343 m^3     9.55  58784.3 m^3  53169.2 m^3  64399.5 m^3
```

The book gets the means 121.15 and 160.34 m³, ``\\bar{d} = 39.19`` m³ with ``s_{xy} = 1\\,807.46``, ``s_{\\bar{d}}^2 = 3.548`` (``s_{\\bar{d}} = 1.884``), ``E_a = 3.749`` m³ and ``E_r = 9.56`` % with the tabulated ``t = 1.99`` (the package uses the exact 1.987), and ``\\hat{D} = 58\\,784`` m³.

# References
Péllico Netto, S. & Brena, D. A. (1997). *Inventário Florestal*. Curitiba, chapter 11.
"""
function completereplacementsampling(value1::AbstractVector, value2::AbstractVector, plotarea::Area,
  npop1::Integer, npop2::Integer; α::Real=0.95)

  length(value1) == length(value2) ||
    throw(DimensionMismatch("value1 and value2 must have the same length: the same plots are remeasured on both occasions."))

  o1 = _occasionmean(value1, npop1)
  o2 = _occasionmean(value2, npop2)
  occasion1 = _occasiontable(o1, npop1, plotarea, α)
  occasion2 = _occasiontable(o2, npop2, plotarea, α)

  n = o1.n
  growth = o2.x̅ - o1.x̅
  growthcovariance = ustrip(cov(value1, value2))
  # with a small sample and strong occasion-to-occasion correlation, the covariance
  # estimate can occasionally overshoot enough to make this formula dip below zero;
  # clamped at zero, since a variance estimate cannot legitimately be negative.
  growthvar = max(ustrip(o1.meanvar) + ustrip(o2.meanvar) - 2 * growthcovariance / n, 0.0)
  df = 2 * (n - 1)
  change = _changetable(growth, growthvar, df, npop2, α)

  return SamplingReport((; occasion1, occasion2, change))
end
