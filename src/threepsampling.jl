"""
    threepsampling(predicted::AbstractVector, measured::AbstractVector, totalarea::Area; e::Real=10, α::Real=0.95)

Estimates the total of a forest and its confidence interval by 3P sampling (probability proportional to prediction; Grosenbaugh, 1965, 1967).

# Description

Every tree of the tract is visited once and given a quick and consistent, not necessarily accurate, visual estimate of the quantity of interest, its ``KPI`` (the key predicted index, typically the volume). In the field the estimate is compared with a random number between ``0`` and a ceiling ``Z`` larger than the ``KPI`` of any tree: the tree is selected, and measured carefully, whenever the estimate is at least the number drawn. So the probability of selection of a tree is proportional to its own prediction, which gives the method its name, and only a small subsample needs the expensive measurement.

This function does not do the selection in the field, which cannot be replaced by the table of random numbers: it takes the cruise already done, the ``predicted`` value of every tree and the ``measured`` value of the selected ones (`missing` for the others) and computes Grosenbaugh's ratio estimator. Prediction needs only to be fast and consistent, while measurement is exact, and so the method typically brings the coefficient of variation from the ~60% of the visual cruise to ~20% (Grosenbaugh, 1965). Since it is a 100% cruise, `predicted` is a census of the population and there is no plot or expansion factor: `totalarea` only turns the total into a rate per hectare.

# Arguments

- `predicted::AbstractVector`: the visual estimate (``KPI``) of every tree of the population (``N`` values). Typically a `Unitful` volume; a plain number stays unitless.
- `measured::AbstractVector`: the exact measurement of the selected trees, in the unit of `predicted`, and `missing` for the other trees. Same length and order as `predicted`.
- `totalarea`: the area of the tract, an `Area`, used only for the rate per hectare.
- `e::Real=10`: the admissible sampling error, in percent, used only for `nreq` and `nadd`.
- `α::Real=0.95`: the confidence level.

# Returns

A one-row `DataFrame` with:

- `X`: the bias-corrected total of the population, ``\\hat{X} = K \\bar{R}``.
- `xha`: the total per hectare, ``\\hat{X} / A``.
- `Rm`: the mean ratio ``\\bar{R}`` between the measurement and the prediction of the ``n`` measured trees, Grosenbaugh's correction ratio.
- `K`: the sum ``\\sum KPI`` of `predicted` over **all** the trees of the population.
- `cv`, `sxm`, `Ea`, `Er`, `Xl`, `Xu`: of the total, treating the ``n`` ratios as a simple random sample of the ``N`` trees (`cv` is that of the ratios, `sxm` is that of the total).
- `population`, `f`: as in [`simplecasualsampling`](@ref), with ``f = n/N``.
- `n`, `N`: the trees measured and the trees predicted.
- `nreq`, `nadd`: the trees required for the admissible error and still to measure, from the coefficient of variation of the ratios.
- `A`: the area of the tract.

# Mathematical basis

For the ``n`` measured trees, with the ratio ``R_i = x_i / KPI_i`` between the measurement and the prediction, and the sum of the predictions of the ``N`` trees ``K = \\sum_{j=1}^{N} KPI_j``:

```math
\\bar{R} = \\frac{\\sum_{i=1}^{n} R_i}{n} \\qquad
s_R^2 = \\frac{\\sum_{i=1}^{n} (R_i - \\bar{R})^2}{n - 1} \\qquad
cv = \\frac{s_R}{\\bar{R}} \\times 100 \\qquad
\\hat{X} = K \\bar{R}
```

The variance of the mean ratio is that of a simple random sample mean, with the correction ``(1 - f)``, ``f = n/N``, dropped when ``1 - f \\ge 0.98``, and it scales with ``K`` to the standard error of the total. With Student's ``t`` of ``n - 1`` degrees of freedom:

```math
s_{\\bar{R}}^2 = \\frac{s_R^2}{n} (1 - f) \\qquad
s_{\\hat{X}} = K \\, s_{\\bar{R}} \\qquad
E_a = t \\, s_{\\hat{X}} \\qquad E_r = \\frac{E_a}{\\hat{X}} \\times 100 \\qquad
IC\\left[\\hat{X} - E_a \\le X \\le \\hat{X} + E_a\\right] = P
```

The required number of trees for the admissible error ``e``, solved iteratively as in [`simplecasualsampling`](@ref), with the coefficient of variation of the ratios:

```math
n_{req} = \\frac{t^2 \\, cv^2}{e^2 + t^2 \\, cv^2 / N} \\quad \\text{(finite)} \\qquad
n_{req} = \\frac{t^2 \\, cv^2}{e^2} \\quad \\text{(infinite)}
```

# Technical description

Treating the ratios as an ordinary simple random sample, as this function and the usual teaching of the method do, ignores that the trees were not selected with equal probability; there are more refined estimators (the bias-adjusted form of Grosenbaugh, and the cosmetically calibrated estimator of Brewer) that correct for it, mainly when ``n`` is small relative to ``N``, which are out of the scope here. Due to Grosenbaugh (1965, 1967); Rennie (1979) describes the "Point-3P", the two-stage variant that uses [`horizontalpointsampling`](@ref) as the first stage.

# Examples

```julia-repl
julia> predicted = [0.20, 0.45, 0.30, 0.55, 0.15, 0.60, 0.25, 0.40, 0.35, 0.50]u"m^3";

julia> measured  = [missing, 0.48u"m^3", missing, 0.58u"m^3", missing,
                     0.63u"m^3", missing, missing, missing, 0.52u"m^3"];

julia> report = threepsampling(predicted, measured, 5u"ha");

julia> report.X[1]
3.9480113636363634 m^3

julia> report.Rm[1]   # the measured trees were on average 5.3% above their prediction
1.0528030303030302

julia> report.n[1], report.N[1]   # 4 of the 10 trees were selected for measurement
(4, 10)
```

# References
Grosenbaugh, L. R. (1965). Three-pee sampling theory and program 'THRP' for computer generation of selection criteria. USDA Forest Service Research Paper PSW-21.
Grosenbaugh, L. R. (1967). The gains from sample-tree selection with unequal probabilities. *Journal of Forestry*, 65(3), 203-206.
Rennie, J. C. (1979). Point-3P sampling. *Forest Science*, 25(3), 441-447.
"""
function threepsampling(predicted::AbstractVector, measured::AbstractVector, totalarea::Area; e::Real=10, α::Real=0.95)
  length(predicted) == length(measured) || throw(DimensionMismatch("predicted and measured must have the same length."))

  measuredidx = findall(!ismissing, measured)
  n = length(measuredidx)
  n >= 2 || throw(ArgumentError("At least two measured trees are required."))
  N = length(predicted)
  n <= N || throw(ArgumentError("measured cannot select more trees than predicted has ($n > $N)."))

  K = sum(predicted)
  R = [measured[i] / predicted[i] for i in measuredidx]
  any(<(0), R) && throw(ArgumentError("measured/predicted ratios must be non-negative."))

  ratiostats = _occasionmean(R, N)
  R̄ = ratiostats.x̅
  cv = sqrt(ratiostats.s²) / R̄ * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  infinite = _isinfinitepopulation(n, N)
  population = infinite ? "infinite" : "finite"
  requiredn = infinite ? _requiredsamplesize(cv^2, n, e, α) : _requiredsamplesize(cv^2, n, e, α, N)
  missingn = n > requiredn ? 0 : requiredn - n

  total = K * R̄
  se = K * sqrt(ratiostats.meanvar)
  abserr = Ttab * se
  relerr = abserr / total * 100

  DataFrame(
    X=total, xha=total / totalarea,
    Rm=R̄, K=K,
    cv=round(cv, digits=2), sxm=se,
    Ea=abserr, Er=round(relerr, digits=2),
    Xl=total - abserr, Xu=total + abserr,
    population=population, f=round(n / N, digits=3),
    n=n, N=N, nreq=requiredn, nadd=missingn,
    A=totalarea,
  )
end
