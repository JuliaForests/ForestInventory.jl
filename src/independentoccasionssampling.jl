"""
    independentoccasionssampling(value1::AbstractVector, value2::AbstractVector, plotarea::Area,
                                  npop1::Integer, npop2::Integer; α::Real=0.95)

Estimates the growth (change) of a forest between two inventory occasions with two independent samples, generic to any sample sizes.

# Description

The simplest of the sampling designs on successive occasions: occasion 1 and occasion 2 are sampled independently, with no attempt to remeasure the same plots. It is easy to carry out but the least efficient to detect growth, since the correlation between the plots of the two occasions is not used; see [`completereplacementsampling`](@ref) and [`partialreplacementsampling`](@ref) for the designs that remeasure plots. The notation and the formulas follow Péllico Netto & Brena (1997), chapter 11, Sec. 11.1.

# Arguments

- `value1`, `value2`: the independent vectors of plot measurements of occasion 1 and occasion 2 (the lengths may differ). A `Unitful` quantity keeps its unit; a plain number stays unitless.
- `plotarea`: the area of each plot, an `Area`, used only for the mean per hectare.
- `npop1`, `npop2`: the possible plots of the population on occasion 1 and on occasion 2.
- `α::Real=0.95`: the confidence level.

# Returns

A [`SamplingReport`](@ref) with `occasion1`, `occasion2` and `change`, each a one-row `DataFrame`:

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

Since the samples are independent, there is no covariance and the variance of the change is the sum of the variances of the two means (Eqs. 11.20 to 11.28):

```math
\\bar{d} = \\bar{y} - \\bar{x} \\qquad
s_{\\bar{d}}^2 = s_{\\bar{x}}^2 + s_{\\bar{y}}^2 \\qquad
s_{\\bar{d}} = \\sqrt{s_{\\bar{x}}^2 + s_{\\bar{y}}^2} \\qquad
E_a = t \\, s_{\\bar{d}} \\qquad E_r = \\frac{t \\, s_{\\bar{d}}}{\\bar{d}} \\times 100 \\qquad
\\hat{D} = N_2 \\, \\bar{d}
```

The confidence interval of the change is ``IC[\\hat{D} - N_2 t s_{\\bar{d}} \\le D \\le \\hat{D} + N_2 t s_{\\bar{d}}] = P``. The degrees of freedom of ``t`` are ``(n_1 - 1) + (n_2 - 1)`` (text after Eq. 11.27).

# Examples

The example of the book (chapter 11, Quadro 11.1): *Pinus elliottii* in Guarapuava, PR, with 44 independent plots in 1975 and 44 in 1976, in a population of 1 500 units. The volumes of the book are per hectare, so they are the volumes of units of 1 ha.

```julia-repl
julia> v1 = [146.7, 182.1, 258.6, 198.8, 196.7, 155.6, 201.0, 183.8, 212.4, 202.5, 217.5, 138.9, 242.7,
        200.0, 208.3, 193.0, 279.3, 181.0, 223.4, 210.4, 206.1, 171.2, 264.4, 278.6, 253.3, 41.0,
        242.9, 241.6, 251.2, 144.5, 225.9, 178.5, 182.0, 176.2, 199.5, 156.8, 136.2, 225.4,
        184.3, 245.8, 164.9, 198.8, 168.6, 228.7]u"m^3";
julia> v2 = [330.4, 222.3, 300.2, 293.8, 270.2, 230.4, 312.3, 341.6, 326.5, 46.3, 312.8, 273.8, 272.3, 225.2,
        282.9, 233.3, 213.4, 224.2, 229.1, 188.0, 176.7, 221.0, 247.5, 249.5, 230.3, 136.4,
        285.4, 265.9, 248.1, 250.2, 238.0, 265.9, 210.6, 280.0, 198.2, 290.3, 248.2, 235.1,
        265.2, 236.7, 280.4, 285.4, 198.3, 248.0]u"m^3";

julia> report = independentoccasionssampling(v1, v2, 1u"ha", 1500, 1500);

julia> occasion1(report).xm[1]
199.9795454545455 m^3

julia> change(report)
1×7 DataFrame
 Row │ dm           sdm          Ea           Er       D            Dl           Du
     │ Quantity…    Quantity…    Quantity…    Float64  Quantity…    Quantity…    Quantity…
─────┼────────────────────────────────────────────────────────────────────────────────────────
   1 │ 48.2091 m^3  10.2238 m^3  20.3243 m^3    42.16  72313.6 m^3  41827.3 m^3  102800.0 m^3
```

The book gets the means 199.98 and 248.19 m³, ``\\bar{d} = 48.21`` m³ with ``s_{\\bar{d}}^2 = 104.53`` (``s_{\\bar{d}} = 10.22``), ``E_a = 20.34`` m³ and ``E_r = 42.19`` % with the tabulated ``t = 1.99`` (the package uses the exact 1.988, ``E_r = 42.16`` %), and ``\\hat{D} = 72\\,315`` m³.

# References
Péllico Netto, S. & Brena, D. A. (1997). *Inventário Florestal*. Curitiba, chapter 11.
"""
function independentoccasionssampling(value1::AbstractVector, value2::AbstractVector, plotarea::Area,
  npop1::Integer, npop2::Integer; α::Real=0.95)

  o1 = _occasionmean(value1, npop1)
  o2 = _occasionmean(value2, npop2)
  occasion1 = _occasiontable(o1, npop1, plotarea, α)
  occasion2 = _occasiontable(o2, npop2, plotarea, α)

  growth = o2.x̅ - o1.x̅
  growthvar = ustrip(o1.meanvar) + ustrip(o2.meanvar)
  df = (o1.n - 1) + (o2.n - 1)
  change = _changetable(growth, growthvar, df, npop2, α)

  return SamplingReport((; occasion1, occasion2, change))
end
