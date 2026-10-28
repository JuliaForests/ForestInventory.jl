# A plot occupying a full reference area (1 ha in the metric system, 1 ac in the
# imperial one) is treated as effectively infinite regardless of the sampling fraction,
# matching the convention used by the original tool this design was ported from.
_isreferencearea(a::Area) = a == 1 * referencearea(a)

"""
    simplecasualsampling(value::AbstractVector, plotarea::Area, totalarea::Area; e::Real=10, α::Real=0.95)

Estimates the mean, the total and their confidence intervals of a forest by simple random sampling of fixed-area plots, and the number of plots needed for an admissible error.

# Description

In simple random sampling every combination of ``n`` plots out of the ``N`` possible ones has the same chance of being the sample. The forest is treated as a single population, with no division in strata and no grouping of plots, so this is the reference design against which the others are compared: stratified, cluster, systematic and the successive-occasions designs reduce to it in their degenerate cases. The notation and the formulas follow Péllico Netto & Brena (1997), chapter 5.

# Arguments

- `value`: the measurement of every plot (volume, basal area, number of trees, or any variable). A `Unitful` quantity keeps its unit in every result; a plain number stays a plain number, since a bare number carries no unit to assume.
- `plotarea`: the area of each plot, an `Area`. All plots have the same area.
- `totalarea`: the area of the whole forest, an `Area`.
- `e::Real=10`: the admissible sampling error, in percent of the mean (``e``).
- `α::Real=0.95`: the confidence level (``1 - \\alpha`` is the significance level of Student's ``t``).

# Returns

A one-row `DataFrame`. Every column that has a physical unit holds `Unitful` quantities, so `stripunits`/`parseunits` work on it. The columns:

- `xm`: the mean per plot, ``\\bar{x}``.
- `cv`: the coefficient of variation, ``cv`` (%).
- `s2xm`: the variance of the mean, ``s^2_{\\bar{x}}``.
- `sxm`: the standard error, ``s_{\\bar{x}}``.
- `Ea`: the absolute sampling error, ``E_a``.
- `Er`: the relative sampling error, ``E_r`` (%).
- `xha`: the mean per hectare (per acre if imperial), ``\\bar{x}`` expanded by the expansion factor of the plot.
- `X`: the estimated total of the forest, ``\\hat{X}``.
- `Xl`, `Xu`: the lower and upper limits of the confidence interval of the total.
- `population`: `"finite"` or `"infinite"`, the classification used for the variance.
- `f`: the sampling fraction, ``f = n/N``.
- `n`: the number of measured plots.
- `nreq`: the number of plots required for the admissible error `e`.
- `nadd`: the plots still to measure, ``\\max(n_{req} - n, 0)``.
- `N`: the number of possible plots in the population.
- `A`: the area of the forest.

# Mathematical basis

The number of possible plots is ``N = A / a``, with ``A`` the area of the forest and ``a`` the area of a plot, and the sampling fraction is ``f = n/N``. From the ``n`` measured values ``x_i``, the mean, the variance and the coefficient of variation (Eqs. 5.3, 5.5 and 5.14):

```math
\\bar{x} = \\frac{\\sum_{i=1}^{n} x_i}{n} \\qquad
s_x^2 = \\frac{\\sum_{i=1}^{n} (x_i - \\bar{x})^2}{n - 1} \\qquad
cv = \\frac{s_x}{\\bar{x}} \\times 100
```

The variance of the mean has the finite-population correction ``(1 - f)`` (Eqs. 5.9 and 5.10), which is dropped for an infinite population:

```math
s_{\\bar{x}}^2 = \\frac{s_x^2}{n} (1 - f) \\quad \\text{(finite)} \\qquad
s_{\\bar{x}}^2 = \\frac{s_x^2}{n} \\quad \\text{(infinite)} \\qquad
s_{\\bar{x}} = \\sqrt{s_{\\bar{x}}^2}
```

The population is taken as infinite when ``1 - f \\ge 0.98`` (Sec. 5.4) or when the plot is a full reference area (1 ha, or 1 ac). The sampling errors use Student's ``t`` with ``n - 1`` degrees of freedom (Eqs. 5.19 and 5.20), and the mean per hectare uses the expansion factor ``F`` of the plot (``F = 1\\,\\text{ha} / a``):

```math
E_a = t \\, s_{\\bar{x}} \\qquad E_r = \\frac{t \\, s_{\\bar{x}}}{\\bar{x}} \\times 100 \\qquad
x_{ha} = \\bar{x} F
```

The total of the population and its confidence interval (Eqs. 5.27 and 5.28):

```math
\\hat{X} = N \\bar{x} \\qquad
IC\\left[\\hat{X} - N t s_{\\bar{x}} \\le X \\le \\hat{X} + N t s_{\\bar{x}}\\right] = P
```

The sampling intensity for the admissible error ``E = (e/100)\\bar{x}``, written with the coefficient of variation in percent (Eqs. 5.40 and 5.42):

```math
n_{req} = \\frac{t^2 \\, cv^2}{e^2 + t^2 \\, cv^2 / N} \\quad \\text{(finite)} \\qquad
n_{req} = \\frac{t^2 \\, cv^2}{e^2} \\quad \\text{(infinite)}
```

The value of ``t`` depends on ``n_{req} - 1`` degrees of freedom, so the pair is solved iteratively until ``n_{req}`` stops changing (Sec. 5.4.3), and the result is rounded up.

# Technical description

With a single stratum, a single plot per cluster or a full census, every other design of this package gives the numbers of this function. The values are the plots, so the unit of the results is the unit of `value`, and the variance of the mean carries its square. A variable with no physical unit, such as the number of trees of each plot, is given as plain numbers: the mean per hectare comes in ``ha^{-1}``, trees per hectare, and the total is the number of trees.

# Examples

The example of the book (chapter 5): a population of 450 plots of 0.1 ha (45 ha) of *Pinus* sp., sampled with the 70 plots of Quadro 5.2, in m³ per plot.

```julia-repl
julia> v = [ 18.5, 22.8, 15.8, 9.8, 23.5, 20.7, 28.3, 6.9, 28.7, 26.4, 7.6, 10.3, 21.1, 26.3, 39.9, 12.7,
        11.5, 11.0, 18.9, 31.4, 9.6, 9.8, 7.1, 25.9, 25.7, 24.2, 19.0, 25.9, 37.7, 23.3, 24.6,
        32.3, 23.4, 17.1, 8.9, 32.7, 8.1, 29.2, 24.7, 25.6, 13.2, 25.7, 27.0, 20.4, 13.6, 32.4,
        31.3, 34.6, 9.5, 18.7, 28.9, 21.6, 8.3, 23.2, 26.6, 25.6, 8.0, 22.1, 21.9, 23.3, 9.9,
        30.7, 29.2, 28.0, 13.7, 27.6, 37.9, 27.3, 28.0, 30.7]u"m^3";

julia> report = simplecasualsampling(v, 0.1u"ha", 45u"ha"; e=10, α=0.95)
1×17 DataFrame
 Row │ xm           cv       s2xm          sxm           Ea           Er       xha                X            Xl          Xu           population  f        n      nreq   nadd   N      A
     │ Quantity…    Float64  Quantity…     Quantity…     Quantity…    Float64  Quantity…          Quantity…    Quantity…   Quantity…    String      Float64  Int64  Int64  Int64  Int64  Quantity…
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 21.7971 m^3    39.58  0.897789 m^6  0.947518 m^3  1.89025 m^3     8.67  217.971 m^3 ha^-1  9808.71 m^3  8958.1 m^3  10659.3 m^3  finite        0.156     70     56      0    450      45 ha

julia> report.xm[1]
21.79714285714286 m^3
```

The book gets ``\\bar{x} = 21.797`` m³, ``cv = 39.58`` %, ``s_{\\bar{x}} = 0.947`` m³ and ``\\hat{X} = 9\\,808.71`` m³. It takes the tabulated ``t = 1.99`` for 69 degrees of freedom (``E_r = 8.65`` %), while the package uses the exact ``t = 1.995`` (``E_r = 8.67`` %). The pilot inventory of the book, the first 20 plots of Quadro 5.1 (mean 19.605 m³, ``cv = 45.64`` %), requires for ``e = 10`` the 70 plots of the definitive inventory:

```julia-repl
julia> pilot = [18.5, 22.8, 15.8, 9.8, 23.5, 20.7, 28.3, 6.9, 28.7, 26.4, 7.6, 10.3, 21.1, 26.3, 39.9, 12.7,
        11.5, 11.0, 18.9, 31.4]u"m^3";

julia> simplecasualsampling(pilot, 0.1u"ha", 45u"ha").nreq[1]
70
```

# References
Péllico Netto, S. & Brena, D. A. (1997). *Inventário Florestal*. Curitiba, chapter 5.
"""
function simplecasualsampling(value::AbstractVector, plotarea::Area, totalarea::Area;
  e::Real=10, α::Real=0.95)

  n = length(value)
  n > 1 || throw(ArgumentError("At least two sampled plots are required."))

  N = round(Int, NoUnits(totalarea / plotarea))
  f = 1 - n / N
  x̅, s² = mean_and_var(value)
  cv = sqrt(s²) / x̅ * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  infinite = _isinfinitepopulation(n, N) || _isreferencearea(plotarea)
  if infinite
    population = "infinite"
    s²x̅ = s² / n
    requiredplots = _requiredsamplesize(cv^2, n, e, α)
  else
    population = "finite"
    s²x̅ = (s² / n) * f
    requiredplots = _requiredsamplesize(cv^2, n, e, α, N)
  end
  sx̅ = sqrt(s²x̅)
  missingplots = n > requiredplots ? 0 : requiredplots - n

  absoluteerror = Ttab * sx̅
  relativeerror = absoluteerror / x̅ * 100
  hectarevalue = x̅ * expansionfactor(plotarea)
  totalvalue = x̅ * N
  ciupper = totalvalue + N * absoluteerror
  cilower = totalvalue - N * absoluteerror

  return DataFrame(
    xm=x̅, cv=round(cv, digits=2), s2xm=s²x̅, sxm=sx̅,
    Ea=absoluteerror, Er=round(relativeerror, digits=2), xha=hectarevalue,
    X=totalvalue, Xl=cilower, Xu=ciupper, population=population,
    f=round(n / N, digits=3), n=n, nreq=requiredplots,
    nadd=missingplots, N=N, A=totalarea,
  )
end
