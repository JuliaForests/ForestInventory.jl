# Per-point sample: sums the individual tree contributions -- basal area/ha and
# value/ha -- for every point that had at least one "in" tree. Points with zero "in"
# trees are added back as explicit zero observations by the caller (see `nzero` below) --
# omitting them would systematically inflate the mean, a well-known pitfall when
# processing angle-count data.
function _pointtable(point::AbstractVector, ef::AbstractVector{<:Quantity}, g::AbstractVector{<:Area}, value::AbstractVector)
  working = DataFrame(point=point, ef=ef, g=g, value=value)
  return combine(groupby(working, :point, sort=true)) do df
    (n=nrow(df), G=sum(df.ef .* df.g), xha=sum(df.ef .* df.value))
  end
end

"""
    horizontalpointsampling(point::AbstractVector, diameter::AbstractVector{<:Length}, value::AbstractVector, baf::Quantity,
                             npoints::Integer, pointarea::Area, totalarea::Area; e::Real=10, α::Real=0.95)

Estimates the mean, the total and their confidence intervals of a forest by horizontal point sampling (Bitterlich's angle-count method).

# Description

An observer stands at each sample point and, with an angle gauge of fixed basal area factor (`baf`), decides tree by tree whether it is counted: a tree is in when its diameter, seen from the point, subtends an angle at least as wide as the critical angle of the gauge. Every counted tree represents exactly `baf` of basal area per hectare, **whatever its size**, and this is why the method needs no plot radius: the "plot" of a tree grows with the tree. It is much faster than a fixed-area design of comparable precision, at the cost of weighting large trees more (a selection with probability proportional to size, corrected exactly by the expansion factor below). The method is due to Bitterlich (1948).

# Arguments

- `point::AbstractVector`: the sample point of every counted tree, in the order of `diameter`/`value`.
- `diameter::AbstractVector`: the diameter of every counted tree, `Length` quantities (any unit, imperial included).
- `value::AbstractVector`: the measurement of every counted tree (typically the volume). A `Unitful` quantity keeps its unit; a plain number stays unitless.
- `baf::Quantity`: the basal area factor of the gauge, the basal area per area that every counted tree represents, e.g. `2.0u"m^2/ha"` or `10.0u"ft^2/ac"`. It is converted to the system of the diameters.
- `npoints::Integer`: the number of points visited, **the ones without counted trees included**. It must not be smaller than the number of distinct points in `point`; the difference is added as observations of zero, which they must be for the mean not to be inflated.
- `pointarea`: the area that each point represents, typically the cell of the grid of a systematic layout (there is no fixed-radius plot), an `Area`.
- `totalarea`: the area of the forest, an `Area`.
- `e::Real=10`: the admissible sampling error, in percent of the mean.
- `α::Real=0.95`: the confidence level.

# Returns

A [`SamplingReport`](@ref) with `pointtable`, one row per point with at least one counted tree, with `n` (trees counted), `G` (basal area per hectare) and `xha` (value per hectare), and `resulttable`, one row with:

- `xha`: the estimated mean per hectare of `value` (per acre if imperial), ``\\bar{x}_{ha}``. It is the primary estimator: the errors and the sample size are computed from it, as `xm` drives them in [`simplecasualsampling`](@ref).
- `cv`, `s2xm`, `sxm`, `Ea`, `Er`: as in [`simplecasualsampling`](@ref), over the `n` per-point values (the points without counted trees included).
- `X`, `Xl`, `Xu`: the estimated total of the forest and the limits of its confidence interval.
- `G`, `Gt`: the mean basal area per hectare, ``G``, and the total basal area of the forest, ``G A``.
- `population`, `f`: the classification of the population and the sampling fraction ``f = n/N``.
- `n`, `nreq`, `nadd`: the points measured, required and still to measure.
- `N`: the possible points of the population, ``N = A / a``.
- `A`: the area of the forest.

# Mathematical basis

For a counted tree of diameter ``D_i``, its basal area and its expansion factor, in trees per hectare (the number of trees per hectare that it stands for), are

```math
g_i = \\frac{\\pi D_i^2}{4} \\qquad EF_i = \\frac{BAF}{g_i}
```

Large trees have a large ``g_i`` and so stand for few trees per hectare, which exactly offsets how easy they are to count from afar. At point ``j``, with ``n_j`` counted trees, the basal area and the value per hectare are the sums of the contribution of every counted tree:

```math
G_j = \\sum_i EF_i \\, g_i = n_j \\cdot BAF \\qquad
x_{ha,j} = \\sum_i EF_i \\, x_i
```

``G_j = n_j \\, BAF`` is exact, since ``EF_i g_i = BAF`` for every tree: the basal area comes from the **count** of trees, without measuring a diameter. The mean basal area per hectare is ``G = \\sum G_j / n`` over the ``n`` points, those without counted trees included, and the total basal area is ``G_t = G A``.

The ``n`` per-point values ``x_{ha,j}`` (zero for the points without counted trees) are the observations of a simple random sample, and each one is already a rate per hectare. With ``N = A / a`` the possible points (``a`` the area of the point) and ``f = n/N``, the mean per hectare, the variance, the coefficient of variation and the variance of the mean (Eqs. 5.3, 5.5, 5.14 and 5.9):

```math
x_{ha} = \\bar{x}_{ha} = \\frac{\\sum_{j=1}^{n} x_{ha,j}}{n} \\qquad
s^2 = \\frac{\\sum_{j=1}^{n} (x_{ha,j} - \\bar{x}_{ha})^2}{n - 1} \\qquad
cv = \\frac{s}{\\bar{x}_{ha}} \\times 100 \\qquad
s_{\\bar{x}}^2 = \\frac{s^2}{n} (1 - f)
```

The correction ``(1 - f)`` is dropped when ``1 - f \\ge 0.98`` (infinite population). With Student's ``t`` of ``n - 1`` degrees of freedom, the errors, the total and its interval; the total is the mean per hectare times the area of the forest, not times ``N``, because the observations are rates per hectare (Eqs. 5.19, 5.20, 5.27 and 5.28):

```math
s_{\\bar{x}} = \\sqrt{s_{\\bar{x}}^2} \\qquad
E_a = t \\, s_{\\bar{x}} \\qquad E_r = \\frac{t \\, s_{\\bar{x}}}{\\bar{x}_{ha}} \\times 100 \\qquad
\\hat{X} = A \\, \\bar{x}_{ha} \\qquad
IC\\left[\\hat{X} - A t s_{\\bar{x}} \\le X \\le \\hat{X} + A t s_{\\bar{x}}\\right] = P
```

The required number of points for the admissible error ``e`` (percent of the mean), solved iteratively like that of [`simplecasualsampling`](@ref) and rounded up:

```math
n_{req} = \\frac{t^2 \\, cv^2}{e^2 + t^2 \\, cv^2 / N} \\quad \\text{(finite)} \\qquad
n_{req} = \\frac{t^2 \\, cv^2}{e^2} \\quad \\text{(infinite)}
```

# Technical description

The population is ``N = A / a`` points, with the area of the point `pointarea` as the cell of the grid that each point controls, and the classification of the population follows [`simplecasualsampling`](@ref). A point visited and without counted trees is an observation of zero, and leaving it out inflates the mean, a well known error of the processing of angle counts. See Husch, Beers & Kershaw (2003, ch. 15) or Avery & Burkhart (2001, ch. 8) for the derivation of the angle-count sampling.

# Examples

```julia-repl
julia> point    = [1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 5, 5, 5];

julia> diameter = [25.0, 30.0, 20.0, 28.0, 22.0, 35.0, 30.0, 25.0, 20.0, 22.0, 30.0, 28.0, 26.0]u"cm";

julia> value    = [0.35, 0.55, 0.22, 0.48, 0.28, 0.85, 0.55, 0.35, 0.22, 0.28, 0.55, 0.48, 0.40]u"m^3";
       # point 6 was visited but had no counted trees: it is not in `point`

julia> report = horizontalpointsampling(point, diameter, value, 2.0u"m^2/ha", 6, 0.1u"ha", 10u"ha");

julia> resulttable(report).xha[1]
32.766571189807216 m^3 ha^-1

julia> resulttable(report).G[1]   # 13 trees x 2 m²/ha / 6 points
4.333333333333333 m^2 ha^-1

julia> pointtable(report)
5×4 DataFrame
 Row │ point  n      G              xha
     │ Int64  Int64  Quantity…      Quantity…
─────┼────────────────────────────────────────────────
   1 │     1      3  6.0 m^2 ha^-1  43.8277 m^3 ha^-1
   2 │     2      2  4.0 m^2 ha^-1  30.3224 m^3 ha^-1
   3 │     3      4  8.0 m^2 ha^-1  61.4972 m^3 ha^-1
   4 │     4      1  2.0 m^2 ha^-1  14.7317 m^3 ha^-1
   5 │     5      3  6.0 m^2 ha^-1  46.2204 m^3 ha^-1
```

Point 6 never appears in `pointtable`, but it is one of the `npoints = 6` observations of `resulttable`: dropping it would inflate `xha`.

# References
Bitterlich, W. (1948). Die Winkelzählprobe. *Allgemeine Forst- und Holzwirtschaftliche Zeitung*, 59, 4-5.
Husch, B., Beers, T. W. & Kershaw, J. A. (2003). *Forest Mensuration*, 4th ed. Wiley, ch. 15.
"""
function horizontalpointsampling(point::AbstractVector, diameter::AbstractVector{<:Length}, value::AbstractVector, baf::Quantity,
  npoints::Integer, pointarea::Area, totalarea::Area; e::Real=10, α::Real=0.95)

  checklengths(:point => point, :diameter => diameter, :value => value)
  ustrip(baf) > 0 || throw(ArgumentError("baf must be a positive value."))

  d = diameter
  g = sectionalarea.(d)
  bafunit = unit(g[1]) / referencearea(unit(d[1]))
  bafQ = uconvert(bafunit, baf)
  ef = bafQ ./ g

  pointtable = _pointtable(point, ef, g, value)
  ncounted = nrow(pointtable)
  nzero = npoints - ncounted
  nzero >= 0 || throw(ArgumentError(
    "npoints ($npoints) is smaller than the number of distinct points with at least " *
    "one tree in point ($ncounted)."))

  n = npoints
  n > 1 || throw(ArgumentError("At least two sampled points are required."))

  valueperha = nzero == 0 ? pointtable.xha : vcat(pointtable.xha, fill(0.0, nzero) .* unit(pointtable.xha[1]))
  basalareaperha = nzero == 0 ? pointtable.G : vcat(pointtable.G, fill(0.0, nzero) .* unit(pointtable.G[1]))

  N = round(Int, NoUnits(totalarea / pointarea))
  f = 1 - n / N
  x̅, s² = mean_and_var(valueperha)
  cv = sqrt(s²) / x̅ * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  infinite = _isinfinitepopulation(n, N)
  if infinite
    population = "infinite"
    s²x̅ = s² / n
    requiredpoints = _requiredsamplesize(cv^2, n, e, α)
  else
    population = "finite"
    s²x̅ = (s² / n) * f
    requiredpoints = _requiredsamplesize(cv^2, n, e, α, N)
  end
  sx̅ = sqrt(s²x̅)
  missingpoints = n > requiredpoints ? 0 : requiredpoints - n

  absoluteerror = Ttab * sx̅
  relativeerror = absoluteerror / x̅ * 100
  totalvalue = x̅ * totalarea
  ciupper = totalvalue + totalarea * absoluteerror
  cilower = totalvalue - totalarea * absoluteerror
  basalareamean = mean(basalareaperha)
  basalareatotal = basalareamean * totalarea

  resulttable = DataFrame(
    xha=x̅, cv=round(cv, digits=2), s2xm=s²x̅, sxm=sx̅,
    Ea=absoluteerror, Er=round(relativeerror, digits=2), X=totalvalue,
    Xl=cilower, Xu=ciupper, G=basalareamean, Gt=basalareatotal,
    population=population, f=round(n / N, digits=3), n=n, nreq=requiredpoints,
    nadd=missingpoints, N=N, A=totalarea,
  )

  return SamplingReport((; pointtable=pointtable, resulttable=resulttable))
end
