# Two-stage sampling's variance decomposition, unlike one-stage cluster sampling's
# `_clustervariance`, has two independent sources of sampling error: which `n` (of `npop`)
# primary units were drawn, and which `m` (of `mpop`) secondary units were measured within
# each of them (Péllico Netto & Brena, Eq. 8.14, equal-size primaries, SRS at both stages). The
# correction of each stage is dropped when its sampling fraction is negligible (Eqs. 8.15, 8.16).
function _twostagevariance(means::AbstractVector{<:Real}, vars::AbstractVector{<:Real}, m::Real, mpop::Real, npop::Real)
  n = length(means)
  ȳ = mean(means)
  withinvar = mean(vars)
  msb = m * sum(abs2, means .- ȳ) / (n - 1)
  betweenvar = (msb - withinvar) / m
  fpc1 = _isinfinitepopulation(n, npop) ? one(npop) : (npop - n) / npop
  fpc2 = _isinfinitepopulation(m, mpop) ? one(mpop) : (mpop - m) / mpop
  meanvar = fpc1 * betweenvar / n + fpc2 * withinvar / (n * m)
  varterm = betweenvar + withinvar / m
  correctionvarterm = betweenvar + withinvar / mpop
  return (; withinvar, betweenvar, meanvar, varterm, correctionvarterm)
end

"""
    twostagesampling(primary::AbstractVector, value::AbstractVector, plotarea::Area, npop::Integer, mpop::Integer;
                      e::Real=10, α::Real=0.95)

Estimates the mean, the total and their confidence intervals of a forest by two-stage sampling, generic to any number of primary units.

# Description

Two-stage sampling draws ``n`` primary units (stands, blocks, strips) out of ``N`` and, inside each drawn one, measures ``m`` secondary units (plots) out of its ``M`` possible ones. Unlike [`clustersampling`](@ref), which measures all the plots of a drawn cluster, the sub-sampling gives up some precision for less fieldwork. The notation and the formulas follow Péllico Netto & Brena (1997), chapter 8.

# Arguments

- `primary::AbstractVector`: the primary unit of every plot, in the order of `value`.
- `value::AbstractVector`: the measurement of every plot. A `Unitful` quantity keeps its unit; a plain number stays unitless.
- `plotarea`: the area of each plot, an `Area`, used only for the mean per hectare.
- `npop::Integer`: the number of primary units of the population.
- `mpop::Integer`: the number of secondary units of each primary unit of the population (equal for all).
- `e::Real=10`: the admissible sampling error, in percent of the mean.
- `α::Real=0.95`: the confidence level.

# Returns

A [`SamplingReport`](@ref) with `primarytable`, one row per primary unit with `n` (plots measured, ``m``), `xmi` (mean ``\\bar{x}_i``) and `s2i` (variance of the plots), and `resulttable`, one row with:

- `xm`: the mean per secondary unit, ``\\bar{x}``.
- `cv`, `s2d`, `s2e`: the coefficient of variation (%), the variance within the primary units ``s_d^2`` and between them ``s_e^2``.
- `s2xm`, `sxm`, `Ea`, `Er`, `xha`, `X`, `Xl`, `Xu`, `population`, `f`: as in [`simplecasualsampling`](@ref), with the total ``N M \\bar{x}``.
- `m`, `M`: the secondary units measured and possible in each primary unit.
- `n`, `nreq`, `nadd`: the primary units measured, required and still to measure.
- `N`: the primary units of the population.

# Mathematical basis

With ``x_{ij}`` the secondary unit ``j`` of the primary unit ``i``, the mean per secondary unit, the mean of each primary unit and the variance (Eqs. 8.2, 8.4 and 8.6):

```math
\\bar{x} = \\frac{\\sum_{i=1}^{n}\\sum_{j=1}^{m} x_{ij}}{n m} \\qquad
\\bar{x}_i = \\frac{\\sum_{j=1}^{m} x_{ij}}{m} \\qquad
s_x^2 = \\frac{\\sum_{i=1}^{n}\\sum_{j=1}^{m} (x_{ij} - \\bar{x})^2}{n m - 1}
```

The total variance splits in a part between the primary units, ``s_e^2``, and a part within them, ``s_d^2``, by the analysis of variance (Eqs. 8.8 to 8.12):

```math
s_x^2 = s_e^2 + s_d^2 \\qquad
MQ_{d} = \\frac{\\sum_{i}\\sum_{j} (x_{ij} - \\bar{x}_i)^2}{n(m - 1)} = s_d^2 \\qquad
MQ_{e} = \\frac{\\sum_{i} m (\\bar{x}_i - \\bar{x})^2}{n - 1} \\qquad
s_e^2 = \\frac{MQ_{e} - MQ_{d}}{m}
```

The variance of the mean has the correction of each stage, ``(N - n)/N`` and ``(M - m)/M`` (Eq. 8.14), each one dropped (Eqs. 8.15 and 8.16) when ``1 - f \\ge 0.98``; the standard error and the coefficient of variation are

```math
s_{\\bar{x}}^2 = \\frac{N - n}{N} \\frac{s_e^2}{n} + \\frac{M - m}{M} \\frac{s_d^2}{n m} \\qquad
s_{\\bar{x}} = \\sqrt{s_{\\bar{x}}^2} \\qquad
cv = \\frac{\\sqrt{s_e^2 + s_d^2}}{\\bar{x}} \\times 100
```

With ``m = M`` the second stage vanishes and the formula becomes that of [`clustersampling`](@ref). The sampling errors, the total and its interval (Eqs. 8.18 to 8.22), with ``t`` of ``n - 1`` degrees of freedom:

```math
E_a = t \\, s_{\\bar{x}} \\qquad E_r = \\frac{t \\, s_{\\bar{x}}}{\\bar{x}} \\times 100 \\qquad
\\hat{X} = N M \\bar{x} \\qquad
IC\\left[\\hat{X} - N M t s_{\\bar{x}} \\le X \\le \\hat{X} + N M t s_{\\bar{x}}\\right] = P
```

The required number of primary units for ``E = (e/100)\\bar{x}``, solved iteratively and rounded up (Eqs. 8.25 and 8.26):

```math
n_{req} = \\frac{t^2 \\left(s_e^2 + s_d^2/m\\right)}{E^2} \\quad \\text{(infinite)} \\qquad
n_{req} = \\frac{t^2 \\left(s_e^2 + s_d^2/m\\right)}{E^2 + \\frac{t^2}{N}\\left(s_e^2 + s_d^2/M\\right)} \\quad \\text{(finite)}
```

# Technical description

The correction of the finite population in the sampling intensity uses the secondary units of the population, ``M``, while the estimate uses those measured, ``m``: it is the only design of this package where the two differ. The same number of secondary units must be measured in every primary unit.

# Examples

The pilot inventory of the book (chapter 8, Quadro 8.2): a forest of 1 000 000 ha, with ``N = 100`` primary units of 10 000 ha, each of ``M = 40\\,000`` plots of 0.25 ha, and 18 primary units with 6 plots each, in m³ per plot. Each row of the matrix is a primary unit.

```julia-repl
julia> rows = [29.07 72.37 96.03 55.90 106.91 54.51;
        113.55 53.14 63.08 125.64 113.10 27.21;
        37.26 58.57 40.27 58.81 57.30 72.97;
        52.66 58.97 37.86 74.61 69.56 70.67;
        60.81 37.76 63.43 75.40 46.24 44.19;
        72.98 65.00 94.62 55.52 57.81 48.97;
        120.71 73.83 108.75 47.76 101.33 98.52;
        62.48 70.30 41.59 49.85 65.64 65.05;
        73.35 38.46 138.26 106.78 154.13 145.69;
        56.35 19.51 20.25 23.54 38.34 25.06;
        73.65 51.87 62.15 43.96 30.82 41.77;
        65.58 36.48 38.76 26.37 83.33 47.14;
        34.31 29.64 30.69 23.87 23.63 37.50;
        45.98 50.52 44.55 72.07 42.39 86.99;
        55.58 67.72 44.15 107.90 51.05 101.61;
        26.79 27.88 69.14 36.53 41.43 55.25;
        52.11 47.06 80.57 20.37 166.67 118.35;
        44.40 117.04 82.26 130.69 123.53 131.59;
        41.48 90.94 60.43 47.60 58.71 44.05;
        59.69 65.41 72.89 55.52 41.42 76.90;
        50.27 52.25 31.06 79.50 88.10 60.19;
        62.24 27.18 69.97 72.08 60.97 46.94;
        36.71 31.28 34.94 37.25 64.99 40.83;
        63.36 26.13 25.89 33.80 25.12 60.89;
        43.08 41.36 30.66 32.78 28.10 30.42;
        38.29 37.51 29.58 62.67 56.40 73.09;
        32.34 57.92 89.36 16.09 52.05 54.40;
        75.36 59.68 36.97 22.02 40.25 69.87;
        53.63 21.01 46.76 84.00 88.49 45.35;
        48.85 45.06 30.19 30.51 27.35 48.72;
        21.05 79.97 20.08 62.15 16.20 56.69;
        29.13 27.72 35.92 80.50 98.08 25.65;
        24.77 20.32 62.24 27.18 69.97 72.08;
        47.21 43.56 25.74 23.53 28.29 50.16;
        43.11 40.80 53.10 28.08 54.64 19.82;
        40.25 33.02 32.24 23.87 64.41 18.85];
julia> value = vec(permutedims(rows[1:18, :])) * u"m^3";
julia> primary = repeat(1:18, inner=6);

julia> report = twostagesampling(primary, value, 0.25u"ha", 100, 40000);
julia> rt = resulttable(report);

julia> rt.xm[1]
64.74046296296297 m^3

julia> rt.s2d[1]
697.2313972222222 m^6

julia> rt.s2e[1]
389.50438209876535 m^6

julia> rt.nreq[1]
36
```

The book gets ``\\bar{x} = 64.740`` m³, ``MQ_d = s_d^2 = 697.231``, ``s_e^2 = 389.505`` and 36 primary units, which it then measures (the 36 rows of the matrix, the first 18 being the pilot). The definitive inventory gives ``\\bar{x} = 55.877`` m³ and ``s_d^2 = 535.637``; the book prints ``s_e^2 = 280.555`` and ``s_{\\bar{x}}^2 = 7.467``, but its own means give ``s_e^2 = 276.228`` and ``s_{\\bar{x}}^2 = 7.390``, which the package computes. The cell of the third plot of the 32nd primary unit is printed as 125.77 in the book, but the means and sums of squares of the book need 35.92, the value of the matrix.

# References
Péllico Netto, S. & Brena, D. A. (1997). *Inventário Florestal*. Curitiba, chapter 8.
Cochran, W. G. (1977). *Sampling Techniques*, 3rd ed. New York: Wiley, Eq. 10.4.
"""
function twostagesampling(primary::AbstractVector, value::AbstractVector, plotarea::Area, npop::Integer, mpop::Integer;
  e::Real=10, α::Real=0.95)

  checklengths(:primary => primary, :value => value)

  table = _clustertable(primary, value; groupname=:primary)
  m = _equalclustersize(table.n)
  n = nrow(table)
  m <= mpop || throw(ArgumentError("The number of measured secondary units per primary unit (m=$m) cannot exceed the population size mpop=$mpop."))

  x̅ = mean(table.xmi)
  vunit = unit(x̅)
  v2unit = vunit^2

  variance = _twostagevariance(ustrip.(table.xmi), ustrip.(table.s2i), m, mpop, npop)
  cv = sqrt(variance.withinvar + variance.betweenvar) / ustrip(x̅) * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  infinite = _isinfinitepopulation(n, npop)
  E = (e / 100) * x̅
  requiredprimaries = _requiredsamplesize(variance.varterm, n, ustrip(E), α, infinite ? nothing : npop;
    correctionvarterm=variance.correctionvarterm)
  missingprimaries = n > requiredprimaries ? 0 : requiredprimaries - n

  population = infinite ? "infinite" : "finite"
  sx̅ = sqrt(variance.meanvar) * vunit
  absoluteerror = Ttab * sx̅
  relativeerror = absoluteerror / x̅ * 100
  hectarevalue = x̅ * expansionfactor(plotarea)
  totalvalue = x̅ * npop * mpop
  ciupper = totalvalue + npop * mpop * absoluteerror
  cilower = totalvalue - npop * mpop * absoluteerror

  resulttable = DataFrame(
    xm=x̅, cv=round(cv, digits=2), s2d=variance.withinvar * v2unit,
    s2e=variance.betweenvar * v2unit, s2xm=variance.meanvar * v2unit, sxm=sx̅, Ea=absoluteerror,
    Er=round(relativeerror, digits=2), xha=hectarevalue, X=totalvalue,
    Xl=cilower, Xu=ciupper, population=population, f=round(n / npop, digits=3),
    m=m, M=mpop, n=n, nreq=requiredprimaries, nadd=missingprimaries, N=npop,
  )

  return SamplingReport((; primarytable=table, resulttable=resulttable))
end
