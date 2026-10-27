# Per-line Strand sample: sums each tallied tree's per-tree factor (trees/ha, basal
# area/ha, value/ha -- see `strandsampling`'s docstring for why, unlike Bitterlich's
# point sampling, every one of these depends on the tree's own diameter) for every line
# that had at least one "in" tree. Lines with zero "in" trees are added back as explicit
# zero observations by the caller, exactly as `_pointtable`/`horizontalpointsampling`
# already does -- omitting them would inflate the mean the same way here.
function _linetable(line::AbstractVector, tf::AbstractVector{<:Quantity}, g::AbstractVector{<:Area}, value::AbstractVector)
  working = DataFrame(line=line, tf=tf, g=g, value=value)
  return combine(groupby(working, :line, sort=true)) do df
    (n=nrow(df), N=sum(df.tf), G=sum(df.tf .* df.g), xha=sum(df.tf .* df.value))
  end
end

"""
    strandsampling(line::AbstractVector, diameter::AbstractVector{<:Length}, value::AbstractVector, baf::Quantity, linelength::Length,
                   nlines::Integer, totalarea::Area; e::Real=10, α::Real=0.95)

Estimates the mean, the total and their confidence intervals of a forest by horizontal line sampling (Strand's method).

# Description

An observer walks a straight line of fixed length `linelength` and projects, at every point of it, the same angle gauge of [`horizontalpointsampling`](@ref) perpendicular to the line, counting every tree, on either side, whose diameter subtends at least the critical angle. The inclusion zone of a tree of diameter ``D_i`` is a strip of width proportional to ``D_i`` along the whole line (instead of a circle), so the reach of a counted tree grows linearly with its diameter and not quadratically as in point sampling: no per-tree factor is independent of size. The factors of tree density and of basal area depend on the diameter of the tree, and **it is the sum of the diameters, not of the basal areas, that the method estimates directly**. Each line gives one observation per hectare, as each point does in [`horizontalpointsampling`](@ref), and the ``n`` lines are a simple random sample with the same exact confidence-interval calculations: it is a proper design with selection proportional to size, unlike the approximate [`prodansampling`](@ref). The method is due to Strand (1957).

# Arguments

- `line::AbstractVector`: the sample line of every counted tree, in the order of `diameter`/`value`.
- `diameter::AbstractVector`: the diameter of every counted tree, `Length` quantities. **Only metric results are supported**: diameters and `linelength` are converted to `cm` and `m` and the results are per hectare, even when the input is imperial.
- `value::AbstractVector`: the measurement of every counted tree (typically the volume). A `Unitful` quantity keeps its unit; a plain number stays unitless.
- `baf::Quantity`: the basal area factor of the gauge, e.g. `2.0u"m^2/ha"`, the same constant that [`horizontalpointsampling`](@ref) takes.
- `linelength::Length`: the length of every sample line.
- `nlines::Integer`: the number of lines walked, **the ones without counted trees included**, treated as the `npoints` of [`horizontalpointsampling`](@ref).
- `totalarea`: the area of the forest, an `Area`.
- `e::Real=10`: the admissible sampling error, in percent of the mean.
- `α::Real=0.95`: the confidence level.

# Returns

A [`SamplingReport`](@ref) with `pointtable`, one row per line with at least one counted tree (it shares the accessor [`pointtable`](@ref) with [`horizontalpointsampling`](@ref)), with `n` (trees counted), `N` (trees per hectare), `G` (basal area per hectare) and `xha` (value per hectare), and `resulttable`, one row with:

- `xha`: the estimated mean per hectare of `value` (per acre if imperial), ``\\bar{x}_{ha}``. It is the primary estimator: the errors and the sample size are computed from it, as `xm` drives them in [`simplecasualsampling`](@ref).
- `cv`, `s2xm`, `sxm`, `Ea`, `Er`: as in [`simplecasualsampling`](@ref), over the `n` per-line values (the lines without counted trees included).
- `X`, `Xl`, `Xu`: the estimated total of the forest and the limits of its confidence interval.
- `G`, `Gt`: the mean basal area per hectare, ``G``, and the total basal area of the forest, ``G A``.
- `population`, `f`: the classification of the population and the sampling fraction ``f = n/N``.
- `n`, `nreq`, `nadd`: the lines measured, required and still to measure.
- `N`: the possible lines of the population, ``N = A / 1\\,\\text{ha}``.
- `A`: the area of the forest.

Here the possible lines are counted as hectares: ``N = A / 1\\,\\text{ha}``.

# Mathematical basis

The dimensionless constant of the gauge, ``k = 2\\sin(\\theta/2)``, relates to the basal area factor by ``BAF = 2500 k^2`` (Kershaw et al., 2017, Table A.6.B), so ``k = \\sqrt{BAF/2500}``. For a counted tree of diameter ``D_i`` (cm) on a line of length ``L`` (m), the tree factor, in trees per hectare per tree counted, and the factors of basal area and of any measurement ``x_i`` (typically the volume) are

```math
TF_i = \\frac{500\\,000\\,k}{D_i \\, L} \\qquad
G_i = TF_i \\, g_i = \\frac{12.5\\pi k D_i}{L} \\qquad
x_{ha,i} = TF_i \\, x_i = \\frac{500\\,000\\,k\\,x_i}{D_i \\, L}
```

while the diameter factor does not depend on ``D_i``, ``D_i \\, TF_i = 500\\,000\\,k/L``: this is the sense in which the sum of the diameters per unit area is estimated by the count alone (Kershaw et al., 2017, p. 298). The values of a line ``j`` are the sums over its counted trees, ``N_j = \\sum TF_i``, ``G_j = \\sum G_i`` and ``x_{ha,j} = \\sum x_{ha,i}``, and the mean basal area per hectare is ``G = \\sum G_j / n`` over the ``n`` lines, those without counted trees included, with ``G_t = G A``.

The ``n`` per-line values ``x_{ha,j}`` (zero for the lines without counted trees) are the observations of a simple random sample, and each one is already a rate per hectare. With ``N = A / a`` the possible lines (``a`` the area of one hectare) and ``f = n/N``, the mean per hectare, the variance, the coefficient of variation and the variance of the mean (Eqs. 5.3, 5.5, 5.14 and 5.9):

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

The required number of lines for the admissible error ``e`` (percent of the mean), solved iteratively like that of [`simplecasualsampling`](@ref) and rounded up:

```math
n_{req} = \\frac{t^2 \\, cv^2}{e^2 + t^2 \\, cv^2 / N} \\quad \\text{(finite)} \\qquad
n_{req} = \\frac{t^2 \\, cv^2}{e^2} \\quad \\text{(infinite)}
```

# Technical description

Due to Strand (1957); see Kershaw et al. (2017, ch. 9, Sec. 9.3, "Horizontal Line Sampling", and Table A.6.B) for the factors above. A line walked without counted trees is an observation of zero, as a point is in [`horizontalpointsampling`](@ref).

# Examples

```julia-repl
julia> line     = [1, 1, 1, 2, 2, 3];

julia> diameter = [28.1, 33.2, 31.0, 40.1, 27.9, 34.6]u"cm";

julia> value    = [0.65, 1.01, 0.88, 1.78, 0.60, 1.19]u"m^3";

julia> report = strandsampling(line, diameter, value, 2.0u"m^2/ha", 20.0u"m", 4, 10u"ha");
       # line 4 was walked but had no counted trees: it is not in `line`

julia> resulttable(report).xha[1]
32.21365010327445 m^3 ha^-1

julia> pointtable(report)
3×5 DataFrame
 Row │ line   n      N              G                  xha
     │ Int64  Int64  Quantity…      Quantity…          Quantity…
─────┼───────────────────────────────────────────────────────────────────
   1 │     1      3  69.2722 ha^-1  5.12598 m^2 ha^-1  57.9407 m^3 ha^-1
   2 │     2      2  42.9779 ha^-1  3.77645 m^2 ha^-1  46.5944 m^3 ha^-1
   3 │     3      1  20.4366 ha^-1  1.92155 m^2 ha^-1  24.3196 m^3 ha^-1
```

# References
Strand, L. (1957). Relascopic height and volume determination by horizontal line sampling. *Norsk Skogbruk*.
Kershaw, J. A., Ducey, M. J., Beers, T. W. & Husch, B. (2017). *Forest Mensuration*, 5th ed. Wiley, ch. 9.
"""
function strandsampling(line::AbstractVector, diameter::AbstractVector{<:Length}, value::AbstractVector, baf::Quantity, linelength::Length,
  nlines::Integer, totalarea::Area; e::Real=10, α::Real=0.95)

  length(line) == length(diameter) == length(value) ||
    throw(DimensionMismatch("line, diameter, and value must have the same length."))
  ustrip(baf) > 0 || throw(ArgumentError("baf must be a positive value."))

  d = ustrip.(uconvert.(u"cm", diameter))
  L = ustrip(uconvert(u"m", linelength))
  g = basalarea.(diameter)

  k = sqrt(ustrip(uconvert(u"m^2/ha", baf)) / 2500)
  tf = (500_000k ./ (d .* L)) .* u"ha^-1"

  linetable = _linetable(line, tf, g, value)
  ncounted = nrow(linetable)
  nzero = nlines - ncounted
  nzero >= 0 || throw(ArgumentError(
    "nlines ($nlines) is smaller than the number of distinct lines with at least " *
    "one tree in line ($ncounted)."))

  n = nlines
  n > 1 || throw(ArgumentError("At least two sampled lines are required."))

  valueperha = nzero == 0 ? linetable.xha : vcat(linetable.xha, fill(0.0, nzero) .* unit(linetable.xha[1]))
  basalareaperha = nzero == 0 ? linetable.G : vcat(linetable.G, fill(0.0, nzero) .* unit(linetable.G[1]))

  N = round(Int, ustrip(uconvert(NoUnits, totalarea / 1u"ha")))
  f = 1 - n / N
  x̅, s² = mean_and_var(valueperha)
  cv = sqrt(s²) / x̅ * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  infinite = _isinfinitepopulation(n, N)
  if infinite
    population = "infinite"
    s²x̅ = s² / n
    requiredlines = _requiredsamplesize(cv^2, n, e, α)
  else
    population = "finite"
    s²x̅ = (s² / n) * f
    requiredlines = _requiredsamplesize(cv^2, n, e, α, N)
  end
  sx̅ = sqrt(s²x̅)
  missinglines = n > requiredlines ? 0 : requiredlines - n

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
    population=population, f=round(n / N, digits=3), n=n, nreq=requiredlines,
    nadd=missinglines, N=N, A=totalarea,
  )

  return SamplingReport((; pointtable=linetable, resulttable=resulttable))
end
