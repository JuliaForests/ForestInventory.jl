# Per-primary-unit sample size, mean and internal variance, one row per primary unit
# (cluster/start/primary), generic to any number of them via `groupby`. `groupname` is
# the column name to expose in the returned table (`:cluster`, `:start`, `:primary`) --
# shared by `clustersampling`, `multistartsystematicsampling`, and `twostagesampling`,
# which differ only in what the primary unit is called and how the population size is
# derived from it, not in this statistics step.
function _clustertable(group::AbstractVector, value::AbstractVector; groupname::Symbol=:cluster)
  data = DataFrame(; groupname => group, value=value)
  return combine(groupby(data, groupname, sort=true)) do df
    x̅, s² = mean_and_var(df.value)
    (n=nrow(df), xmi=x̅, s2i=s²)
  end
end

# Validates that every primary unit has the same number of secondary units and returns
# it; `_clustervariance`'s decomposition assumes equal-size primary units.
function _equalclustersize(n::AbstractVector{<:Integer})
  sizes = unique(n)
  length(sizes) == 1 ||
    throw(ArgumentError("Every cluster must contain the same number of secondary units; found sizes $(sort(sizes))."))
  return only(sizes)
end

"""
    clustersampling(cluster::AbstractVector, value::AbstractVector, plotarea::Area, totalarea::Area;
                     e::Real=10, α::Real=0.95)

Estimates the mean, the total and their confidence intervals of a forest by one-stage cluster sampling, generic to any number of clusters.

# Description

In cluster sampling the sampling unit is a group of ``M`` neighboring plots (a cluster, or "conglomerado") instead of one plot. A cluster is cheaper to lay out in the field, but the plots of a cluster are alike, so the variance of the mean is larger than that of a simple random sample of the same number of plots. Every cluster has the same number of plots ``M``, and every plot of a drawn cluster is measured. The notation and the formulas follow Péllico Netto & Brena (1997), chapter 9.

# Arguments

- `cluster::AbstractVector`: the cluster of every plot, in the order of `value`.
- `value::AbstractVector`: the measurement of every plot. A `Unitful` quantity keeps its unit; a plain number stays unitless.
- `plotarea`: the area of each plot, an `Area`.
- `totalarea`: the area of the forest, an `Area`.
- `e::Real=10`: the admissible sampling error, in percent of the mean.
- `α::Real=0.95`: the confidence level.

# Returns

A [`SamplingReport`](@ref) with `clustertable`, one row per cluster with `n` (plots measured), `xmi` (mean ``\\bar{x}_i``) and `s2i` (variance of the plots), and `resulttable`, one row with the columns:

- `xm`: the mean per plot, ``\\bar{x}``.
- `cv`, `s2d`, `s2e`, `s2x`: the coefficient of variation (%), the variance within clusters ``s_d^2``, between clusters ``s_e^2`` and the total ``s_x^2`` (per plot).
- `r`: the intracluster correlation coefficient.
- `sxm`, `Ea`, `Er`, `xha`, `X`, `Xl`, `Xu`, `population`, `f`: as in [`simplecasualsampling`](@ref), with the mean per plot (the total is ``N M \\bar{x}``).
- `M`: the plots per cluster.
- `n`, `nreq`, `nadd`: the clusters measured, required and still to measure.
- `N`: the possible clusters of the population, ``N = A / (a M)``.
- `A`: the area of the forest.

# Mathematical basis

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

With ``M = 1`` (a plot per cluster) it is numerically identical to [`simplecasualsampling`](@ref): there is nothing left to cluster. The population is made of ``N`` possible clusters, not plots, so ``N = A / (a M)``. [`multistartsystematicsampling`](@ref) uses the same estimator, since each random start is a cluster. The ``t`` of the errors has ``n - 1`` degrees of freedom, since the clusters are the sampling units.

# Examples

The example of the book (chapter 9, Quadro 9.1): 25 clusters of 4 plots of 0.25 ha in a stand of 4 000 ha, so ``N = 4\\,000`` possible clusters, in m³ per plot. Each row of the matrix is a cluster.

```julia-repl
julia> rows = [84.7 78.2 63.0 88.5;
        70.1 65.4 59.6 58.1;
        75.6 77.7 66.1 43.5;
        81.3 38.7 62.8 31.1;
        99.5 47.0 82.1 77.4;
        97.3 66.1 57.4 113.9;
        70.2 43.5 48.3 69.3;
        81.7 68.8 22.6 78.2;
        60.9 71.4 75.2 81.1;
        67.3 90.5 72.6 67.1;
        75.4 88.5 104.6 105.5;
        71.0 57.0 79.3 71.8;
        61.8 100.9 69.4 61.0;
        116.9 91.0 36.9 57.3;
        56.0 22.9 94.8 130.2;
        68.2 32.5 68.2 58.6;
        55.4 61.3 76.3 40.0;
        55.0 63.0 56.8 36.8;
        78.8 59.6 68.1 46.8;
        33.6 44.4 34.1 79.6;
        48.8 26.6 67.5 51.0;
        82.4 61.8 23.4 60.2;
        69.8 49.5 53.9 0.0;
        58.4 61.5 45.5 50.2;
        27.4 63.6 69.4 60.8];
julia> value = vec(permutedims(rows)) * u"m^3";
julia> cluster = repeat(1:25, inner=4);

julia> report = clustersampling(cluster, value, 0.25u"ha", 4000u"ha");
julia> rt = resulttable(report);

julia> rt.xm[1]
64.548 m^3

julia> rt.s2e[1]
46.249385416666684 m^6

julia> rt.r[1]
0.0966

julia> rt.s2xm[1]
6.175822750000001 m^6

julia> rt.nreq[1]
17
```

The book gets ``\\bar{x} = 64.548`` m³, ``s_d^2 = 432.58`` m⁶, ``s_e^2 = 46.25`` m⁶, ``r = 0.0966``, ``s_{\\bar{x}}^2 = 6.176`` m⁶, ``\\hat{X} = 1\\,032\\,768`` m³ and 17 clusters. For the errors it takes ``t = 1.99``, while the package uses the ``t`` of the 24 degrees of freedom of the clusters (2.064), so ``E_a = 5.13`` m³ instead of 4.94 m³.

# References
Péllico Netto, S. & Brena, D. A. (1997). *Inventário Florestal*. Curitiba, chapter 9.
Cochran, W. G. (1977). *Sampling Techniques*, 3rd ed. New York: Wiley.
"""
function clustersampling(cluster::AbstractVector, value::AbstractVector, plotarea::Area, totalarea::Area;
  e::Real=10, α::Real=0.95)

  length(cluster) == length(value) || throw(DimensionMismatch("cluster and value must have the same length."))

  table = _clustertable(cluster, value; groupname=:cluster)
  M = _equalclustersize(table.n)
  n = nrow(table)

  N = round(Int, ustrip(uconvert(NoUnits, totalarea / (plotarea * M))))
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

  return SamplingReport((; clustertable=table, resulttable=resulttable))
end
