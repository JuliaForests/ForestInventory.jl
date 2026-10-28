"""
    multistartsystematicsampling(start::AbstractVector, value::AbstractVector, plotarea::Area, totalarea::Area;
                                  e::Real=10, α::Real=0.95)

Estimates the mean, the total and their confidence intervals of a forest by systematic sampling with multiple random starts, generic to any number of starts.

# Description

A single systematic line (see [`systematicsampling`](@ref)) has no replication to estimate its own precision. Laying out several systematic lines, each from an independent random start, makes the design a true probability sample: every line is a systematic subsample of ``M`` plots and the ``n`` lines play the role of clusters drawn from the ``N`` possible starts (Shiue, 1960). The calculations are those of one-stage cluster sampling, with each random start treated as a cluster. The notation and the formulas follow Péllico Netto & Brena (1997), chapter 10, which sends the calculations to chapter 9.

# Arguments

- `start::AbstractVector`: the random start (line) of every plot, in the order of `value`.
- `value::AbstractVector`: the measurement of every plot. A `Unitful` quantity keeps its unit; a plain number stays unitless.
- `plotarea`: the area of each plot, an `Area`.
- `totalarea`: the area of the forest, an `Area`.
- `e::Real=10`: the admissible sampling error, in percent of the mean.
- `α::Real=0.95`: the confidence level.

# Returns

A [`SamplingReport`](@ref) with `starttable`, one row per start with `n` (plots measured), `xmi` (mean ``\\bar{x}_i``) and `s2i` (variance of the plots), and `resulttable`, one row with the same columns as that of [`clustersampling`](@ref):

- `xm`: the mean per plot, ``\\bar{x}``.
- `cv`, `s2d`, `s2e`, `s2x`: the coefficient of variation (%), the variance within clusters ``s_d^2``, between clusters ``s_e^2`` and the total ``s_x^2`` (per plot).
- `r`: the intracluster correlation coefficient.
- `sxm`, `Ea`, `Er`, `xha`, `X`, `Xl`, `Xu`, `population`, `f`: as in [`simplecasualsampling`](@ref), with the mean per plot (the total is ``N M \\bar{x}``).
- `M`: the plots per cluster.
- `n`, `nreq`, `nadd`: the clusters measured, required and still to measure.
- `N`: the possible clusters of the population, ``N = A / (a M)``.
- `A`: the area of the forest.

# Mathematical basis

Each start is a cluster of ``M`` plots, so the estimators are those of [`clustersampling`](@ref), with ``N`` the number of non-overlapping systematic subsamples that tile the population, ``N = A / (a M)``:

With ``n`` clusters of ``M`` plots out of ``N`` possible clusters, ``x_{ij}`` the plot ``j`` of cluster ``i``, the mean per plot, the mean of each cluster and the variance of the plots (Eqs. 9.2, 9.3 and 9.5):

```math
\\bar{x} = \\frac{\\sum_{i=1}^{n}\\sum_{j=1}^{M} x_{ij}}{n M} \\qquad
\\bar{x}_i = \\frac{\\sum_{j=1}^{M} x_{ij}}{M} \\qquad
s_x^2 = \\frac{\\sum_{i=1}^{n}\\sum_{j=1}^{M} (x_{ij} - \\bar{x})^2}{n M - 1}
```

The total variance splits in a part between clusters, ``s_e^2``, and a part within clusters (between the plots of a cluster), ``s_d^2``, by the analysis of variance (Eqs. 9.6 to 9.12):

```math
s_x^2 = s_e^2 + s_d^2 \\qquad
MQ_{d} = \\frac{\\sum_{i}\\sum_{j} (x_{ij} - \\bar{x}_i)^2}{n(M - 1)} = s_d^2 \\qquad
MQ_{e} = \\frac{\\sum_{i} M (\\bar{x}_i - \\bar{x})^2}{n - 1} \\qquad
s_e^2 = \\frac{MQ_{e} - MQ_{d}}{M}
```

The intracluster correlation coefficient measures how alike the plots of a cluster are (Eq. 9.15); ``r = 0`` means a cluster is no more homogeneous than the forest, and the variance of the mean equals that of simple random sampling:

```math
r = \\frac{s_e^2}{s_e^2 + s_d^2}
```

The variance of the mean has the finite-population correction ``(1 - f)``, ``f = n/N``, only in the term between clusters (Eq. 9.17), and it is dropped (Eq. 9.18) when ``1 - f \\ge 0.98``. The standard error and the coefficient of variation follow (Eqs. 9.29 and 9.30):

```math
s_{\\bar{x}}^2 = (1 - f) \\frac{s_e^2}{n} + \\frac{s_d^2}{n M} \\qquad
s_{\\bar{x}} = \\sqrt{s_{\\bar{x}}^2} \\qquad
cv = \\frac{s_x}{\\bar{x}} \\times 100
```

The sampling errors, the total and its interval (Eqs. 9.31, 9.32, 9.34 and 9.35), with ``t`` of ``n - 1`` degrees of freedom:

```math
E_a = t \\, s_{\\bar{x}} \\qquad E_r = \\frac{t \\, s_{\\bar{x}}}{\\bar{x}} \\times 100 \\qquad
\\hat{X} = N M \\bar{x} \\qquad
IC\\left[\\hat{X} - N M t s_{\\bar{x}} \\le X \\le \\hat{X} + N M t s_{\\bar{x}}\\right] = P
```

The required number of clusters for the admissible error ``E = (e/100)\\bar{x}`` (Eq. 9.37, equal to ``t^2 s_x^2 [1 + r(M - 1)] / (M E^2)``), solved iteratively and rounded up; for a finite population the correction is added as in the other designs:

```math
n_{req} = \\frac{t^2 \\left(s_e^2 + s_d^2/M\\right)}{E^2} \\quad \\text{(infinite)} \\qquad
n_{req} = \\frac{t^2 \\left(s_e^2 + s_d^2/M\\right)}{E^2 + t^2 \\left(s_e^2 + s_d^2/M\\right) / N} \\quad \\text{(finite)}
```

# Technical description

It is a separate function from [`clustersampling`](@ref) because of the field procedure, not the statistics: a start is a subsample spread over the whole area, and a cluster is a group of contiguous plots. The variance between starts captures the variation between lines, which is largest when the lines run across the slope or the trend of the terrain. The book recommends at least 5 or 6 random starts, since ``t`` is 12.706 for two starts and 4.303 for three at 95%.

# Examples

The example of the book (chapter 10, Quadro 10.1): 9 random starts of 8 plots of 0.25 ha in a stand of 13 000 ha, so ``N = 6\\,500`` possible starts, in m³ per plot. Each row of the matrix is a start.

```julia-repl
julia> rows = [25.97 34.62 25.97 60.27 49.75 81.64 85.90 72.43;
        54.02 56.77 63.78 61.94 61.28 38.93 44.68 48.30;
        66.43 50.46 43.81 39.63 42.47 38.84 37.13 35.93;
        59.33 56.18 39.01 32.93 95.67 32.82 65.09 46.29;
        62.69 68.88 65.84 85.93 48.13 55.98 70.91 37.90;
        47.85 69.66 78.11 36.71 28.09 53.57 40.34 56.22;
        36.64 39.23 42.81 35.51 27.54 54.78 42.56 88.45;
        34.72 41.21 41.79 34.97 41.47 29.11 25.60 21.99;
        43.13 59.13 92.98 49.14 46.68 42.60 26.17 30.14];
julia> value = vec(permutedims(rows)) * u"m^3";
julia> start = repeat(1:9, inner=8);

julia> report = multistartsystematicsampling(start, value, 0.25u"ha", 13000u"ha");
julia> rt = resulttable(report);

julia> rt.xm[1]
49.76986111111111 m^3

julia> rt.s2e[1]
27.486528540426598 m^6

julia> rt.r[1]
0.0885

julia> rt.s2xm[1]
6.9869658613040135 m^6

julia> rt.nreq[1]
13
```

The book gets ``\\bar{x} = 49.769`` m³, ``s_e^2 = 27.487`` m⁶, ``r = 0.088``, ``s_{\\bar{x}}^2 = 6.987`` m⁶, ``\\hat{X} = 2\\,587\\,988`` m³ and 13 starts. For the errors it takes ``t = 1.996``, while the package uses the ``t`` of the 8 degrees of freedom of the starts (2.306), so ``E_a = 6.10`` m³ instead of 5.28 m³.

# References
Péllico Netto, S. & Brena, D. A. (1997). *Inventário Florestal*. Curitiba, chapter 10.
Shiue, C. J. (1960). Systematic sampling with multiple random starts. *Forest Science*, 6(1), 42-50.
"""
function multistartsystematicsampling(start::AbstractVector, value::AbstractVector, plotarea::Area, totalarea::Area;
  e::Real=10, α::Real=0.95)

  checklengths(:start => start, :value => value)

  table = _clustertable(start, value; groupname=:start)
  M = _equalclustersize(table.n)
  n = nrow(table)

  N = round(Int, NoUnits(totalarea / (plotarea * M)))
  x̅ = mean(table.xmi)
  vunit = unit(x̅)
  v2unit = vunit^2

  variance = _clustervariance(ustrip.(table.xmi), ustrip.(table.s2i), M, N)
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
  hectarevalue = x̅ * expansionfactor(plotarea)
  totalvalue = x̅ * M * N
  ciupper = totalvalue + M * N * absoluteerror
  cilower = totalvalue - M * N * absoluteerror

  resulttable = DataFrame(
    xm=x̅, cv=round(cv, digits=2), s2d=variance.withinvar * v2unit,
    s2e=variance.betweenvar * v2unit, s2x=variance.totalvar * v2unit,
    r=round(variance.icc, digits=4), s2xm=variance.meanvar * v2unit, sxm=sx̅, Ea=absoluteerror,
    Er=round(relativeerror, digits=2), xha=hectarevalue, X=totalvalue,
    Xl=cilower, Xu=ciupper, population=population, f=round(n / N, digits=3),
    M=M, n=n, nreq=requiredplots, nadd=missingplots, N=N, A=totalarea,
  )

  return SamplingReport((; starttable=table, resulttable=resulttable))
end
