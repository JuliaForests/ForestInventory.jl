# Per-stratum descriptive statistics: sample size, mean, variance, standard deviation,
# and the stratum's relative weight (its share of the total population area) together
# with weighted standard deviation/variance, the building blocks of every stratified-
# sampling formula below. Generic to any number of strata via `groupby`.
function _auxiliarytable(stratum::AbstractVector, value::AbstractVector, areas::AbstractVector{<:Area})
  data = DataFrame(stratum=stratum, value=value)
  table = combine(groupby(data, :stratum, sort=true)) do df
    x̅, s² = mean_and_var(df.value)
    (nh=nrow(df), xmh=x̅, s2h=s², sh=sqrt(s²))
  end
  p = areas ./ sum(areas)
  table.Wh = p
  table.Whsh = p .* table.sh
  table.Whs2h = p .* table.s2h
  return table
end

# One-way ANOVA testing whether the strata means differ, via the standard sum-of-squares
# decomposition -- mathematically identical to fitting `value ~ 1 + stratum` and reading
# off its F-test, without needing a modeling/design-matrix dependency for what is, here,
# a single categorical predictor.
function _anova(value::AbstractVector, table::AbstractDataFrame)
  n, k = length(value), nrow(table)
  ȳ = mean(value)
  sst = sum(abs2, value .- ȳ)
  ssb = sum(table.nh .* abs2.(table.xmh .- ȳ))
  ssw = sst - ssb
  dfb, dfw = k - 1, n - k
  msb, msw = ssb / dfb, ssw / dfw
  # a single stratum has no between-strata comparison to make: the F-test is undefined
  # (0 numerator degrees of freedom), so the between-strata row is left blank instead.
  fstat = dfb > 0 ? msb / msw : missing
  pval = dfb > 0 ? ccdf(FDist(dfb, dfw), fstat) : missing
  DataFrame(
    source=["Between strata", "Within strata", "Total"],
    gl=[dfb, dfw, n - 1],
    SQ=[ssb, ssw, sst],
    QM=[dfb > 0 ? msb : missing, msw, missing],
    F=[fstat, missing, missing],
    p=[pval, missing, missing],
  )
end

# Elementwise max(required - measured, 0): how many additional plots each stratum still
# needs, never negative when a stratum already exceeds its required allocation.
_missingplots(measured::AbstractVector{<:Integer}, required::AbstractVector{<:Integer}) = max.(required .- measured, 0)

# Resolves a stratum => area dict into a plain vector, in the ascending-sorted stratum
# order `groupby(..., sort=true)` will itself produce -- a real key-based lookup instead
# of trusting the caller to have supplied the areas in that order already, and errors
# clearly if any stratum present in the data has no matching entry in `strataarea`.
function _resolvestrataareas(stratum::AbstractVector, strataarea::AbstractDict{<:Any,<:Area})
  sortedstrata = sort(unique(stratum))
  missingstrata = setdiff(sortedstrata, keys(strataarea))
  isempty(missingstrata) ||
    throw(ArgumentError("strataarea has no entry for stratum(s) $(missingstrata)."))
  return [strataarea[s] for s in sortedstrata]
end

"""
    stratifiedsampling(stratum::AbstractVector, value::AbstractVector, plotarea::Area,
                        strataarea::AbstractDict{<:Any,<:Area}; allocation::Symbol=:proportional, e::Real=10, α::Real=0.95)

Estimates the mean, the total and their confidence intervals of a forest divided in strata, by stratified random sampling, generic to any number of strata.

# Description

Stratified random sampling divides the population in non-overlapping subpopulations (strata), each more homogeneous than the whole, and samples each one independently. When the means of the strata differ, the variance of the estimated mean is smaller than that of simple random sampling with the same number of plots; when they do not, the two designs are equivalent. The price is knowing the area of each stratum in advance. The notation and the formulas follow Péllico Netto & Brena (1997), chapter 6.

# Arguments

- `stratum::AbstractVector`: the stratum of every plot, in the same order as `value`.
- `value::AbstractVector`: the measurement of every plot. A `Unitful` quantity keeps its unit; a plain number stays unitless.
- `plotarea`: the area of each plot, an `Area`.
- `strataarea::AbstractDict`: the area of each stratum, keyed by the labels of `stratum`, e.g. `Dict(1 => 12.0u"ha", 2 => 8.0u"ha")`. With a dictionary the plots need not be sorted by stratum.
- `allocation::Symbol=:proportional`: how the required number of plots is distributed among the strata: `:proportional` to the weight of the stratum (the first stratification of the book, Eqs. 6.27 to 6.29), or `:neyman`, the optimum allocation with equal costs (Eqs. 6.33 to 6.35).
- `e::Real=10`: the admissible sampling error, in percent of the mean.
- `α::Real=0.95`: the confidence level.

# Returns

A [`SamplingReport`](@ref) with three tables:

- `auxiliarytable`: one row per stratum `h`, with `nh` (plots measured), `xmh` (mean ``\\bar{x}_h``), `s2h` (variance ``s_h^2``), `sh` (standard deviation ``s_h``), `Wh` (weight ``W_h = A_h/A``), `Whsh` (``W_h s_h``) and `Whs2h` (``W_h s_h^2``).
- `anova`: the analysis of variance between and within the strata, with `source`, `gl` (degrees of freedom), `SQ` (sum of squares), `QM` (mean square), `F` and `p`. With a single stratum the line between strata is blank.
- `resulttable`: one row with
  - `xm`, `cv`, `s2xm`, `sxm`, `Ea`, `Er`, `X`, `Xl`, `Xu`, `population`, `f`: as in [`simplecasualsampling`](@ref), for the stratified mean ``\\bar{x}_{st}`` (`cv` is the coefficient of variation of the stratified population, ``s_{st}/\\bar{x}_{st} \\times 100``);
  - `nh`, `nhreq`, `nhadd`: the plots measured, required and still to measure in each stratum, as tuples in the order of the strata;
  - `n`, `nreq`: the total of plots measured and required, and `N` the possible plots of the population;
  - `Ah`: the area of each stratum, a tuple of plain numbers in the unit of `strataarea` (a tuple has no single type that `stripunits`/`parseunits` can strip a unit from, so it is stored without unit).

# Mathematical basis

With ``L`` strata, ``N_h = A_h / a`` possible plots in stratum ``h``, ``N = \\sum N_h``, ``n_h`` measured plots, ``n = \\sum n_h``, the weight of the stratum ``W_h = N_h/N = A_h/A`` and its sampling fraction ``f_h = n_h/N_h``. The mean and the variance of each stratum (Eqs. 6.2 and 6.7), the stratified mean (Eq. 6.4) and the variance of the stratified population (Eq. 6.9):

```math
\\bar{x}_h = \\frac{\\sum_{i=1}^{n_h} x_{hi}}{n_h} \\qquad
s_h^2 = \\frac{\\sum_{i=1}^{n_h} (x_{hi} - \\bar{x}_h)^2}{n_h - 1} \\qquad
\\bar{x}_{st} = \\sum_{h=1}^{L} W_h \\bar{x}_h \\qquad
s_{st}^2 = \\sum_{h=1}^{L} W_h s_h^2
```

The variance of the stratified mean in its expanded form (Eq. 6.12), which drops the second term when the sampling fraction is negligible (``1 - f \\ge 0.98``, Eq. 6.13), the standard error and the coefficient of variation:

```math
s_{\\bar{x}(st)}^2 = \\sum_{h=1}^{L} \\frac{W_h^2 s_h^2}{n_h} - \\sum_{h=1}^{L} \\frac{W_h s_h^2}{N} \\qquad
s_{\\bar{x}(st)} = \\sqrt{s_{\\bar{x}(st)}^2} \\qquad
cv = \\frac{s_{st}}{\\bar{x}_{st}} \\times 100
```

The degrees of freedom of Student's ``t`` are not ``n - 1``: they are the effective number of degrees of freedom of Satterthwaite (1946), Eq. 6.26, with ``g_h = N_h(N_h - n_h)/n_h``, rounded up to the next integer as in the book:

```math
n_e = \\frac{\\left(\\sum_{h=1}^{L} g_h s_h^2\\right)^2}{\\sum_{h=1}^{L} \\frac{g_h^2 s_h^4}{n_h - 1}}
```

The sampling errors, the total and its confidence interval (Eqs. 6.18, 6.19, 6.23 and 6.25):

```math
E_a = t \\, s_{\\bar{x}(st)} \\qquad E_r = \\frac{t \\, s_{\\bar{x}(st)}}{\\bar{x}_{st}} \\times 100 \\qquad
\\hat{X} = N \\bar{x}_{st} \\qquad
IC\\left[\\hat{X} - N t s_{\\bar{x}(st)} \\le X \\le \\hat{X} + N t s_{\\bar{x}(st)}\\right] = P
```

The sampling intensity for the admissible error ``E = (e/100)\\bar{x}_{st}``, solved iteratively like that of [`simplecasualsampling`](@ref) and rounded up. With the proportional allocation, ``n_h = W_h n`` (Eqs. 6.27 to 6.29):

```math
n_{req} = \\frac{t^2 \\sum_h W_h s_h^2}{E^2 + t^2 \\sum_h W_h s_h^2 / N} \\quad \\text{(finite)} \\qquad
n_{req} = \\frac{t^2 \\sum_h W_h s_h^2}{E^2} \\quad \\text{(infinite)}
```

With Neyman's optimum allocation for equal costs, ``n_h = W_h s_h / \\sum_h W_h s_h \\; n`` (Eqs. 6.33 to 6.35):

```math
n_{req} = \\frac{t^2 \\left(\\sum_h W_h s_h\\right)^2}{E^2 + t^2 \\sum_h W_h s_h^2 / N} \\quad \\text{(finite)} \\qquad
n_{req} = \\frac{t^2 \\left(\\sum_h W_h s_h\\right)^2}{E^2} \\quad \\text{(infinite)}
```

The analysis of variance (Quadro 6.1) tests whether the means of the strata differ, with ``\\bar{x}`` the mean of all the plots:

```math
SQ_e = \\sum_{h} n_h (\\bar{x}_h - \\bar{x})^2 \\qquad
SQ_d = \\sum_{h}\\sum_{i} (x_{hi} - \\bar{x}_h)^2 \\qquad
SQ_t = \\sum_{h}\\sum_{i} (x_{hi} - \\bar{x})^2 \\qquad
F = \\frac{SQ_e / (L - 1)}{SQ_d / (n - L)}
```

# Technical description

With a single stratum this reduces to [`simplecasualsampling`](@ref): the line between strata of the analysis of variance has no degrees of freedom and the stratified mean is the sample mean. The columns `Wh`, `Whsh` and `Whs2h` are the pieces of the variance, of the required number of plots and of its allocation, kept in the table so that they can be audited. Every stratum needs at least two plots.

# Examples

The example of the book (chapter 6): 3 strata of the population of 450 plots of 0.1 ha (45 ha), with 144, 164 and 142 plots, and the 22 plots of Quadro 6.2, in m³ per plot.

```julia-repl
julia> stratum = vcat(fill(1, 7), fill(2, 8), fill(3, 7));

julia> value = [15.8, 7.6, 8.8, 12.5, 11.1, 16.2, 12.2, 20.4, 30.5, 26.8, 27.2, 28.4, 19.7, 20.4, 23.1, 21.3,
        24.3, 29.2, 21.8, 33.1, 35.8, 26.7]u"m^3";

julia> area = Dict(1 => 14.4u"ha", 2 => 16.4u"ha", 3 => 14.2u"ha");

julia> report = stratifiedsampling(stratum, value, 0.1u"ha", area);

julia> resulttable(report).xm[1]
21.46506349206349 m^3

julia> resulttable(report).X[1]
9659.278571428571 m^3

julia> resulttable(report).nreq[1]
19

julia> auxiliarytable(report)
3×8 DataFrame
 Row │ stratum  nh     xmh          s2h          sh           Wh        Whsh         Whs2h
     │ Int64    Int64  Quantity…    Quantity…    Quantity…    Float64   Quantity…    Quantity…
─────┼───────────────────────────────────────────────────────────────────────────────────────────
   1 │       1      7  12.0286 m^3  10.4624 m^6  3.23456 m^3  0.32      1.03506 m^3  3.34796 m^6
   2 │       2      8  24.5625 m^3  17.4827 m^6  4.18123 m^3  0.364444  1.52383 m^3  6.37147 m^6
   3 │       3      7  27.4571 m^3  30.8229 m^6  5.55183 m^3  0.315556  1.75191 m^3  9.72632 m^6

julia> anova(report)
3×6 DataFrame
 Row │ source          gl     SQ           QM           F             p
     │ String          Int64  Quantity…    Quantity…?   Float64?      Float64?
─────┼─────────────────────────────────────────────────────────────────────────────────
   1 │ Between strata      2  951.399 m^6    475.7 m^6       24.4219        5.60811e-6
   2 │ Within strata      19   370.09 m^6  19.4784 m^6  missing       missing
   3 │ Total              21  1321.49 m^6      missing  missing       missing
```

The book gets ``\\bar{x}_{st} = 21.465`` m³, ``s_{st}^2 = 19.452`` m⁶, ``s_{\\bar{x}} = 0.916`` m³, ``E_r = 9.0`` %, ``\\hat{X} = 9\\,659`` m³, ``F = 24.42`` and 18 plots (it rounds the 18.05 plots that the iteration gives down, the package up). The analysis of variance of the book has ``SQ_e = 951.355``, the package 951.399, since the book rounds the means.

# References
Péllico Netto, S. & Brena, D. A. (1997). *Inventário Florestal*. Curitiba, chapter 6.
Satterthwaite, F. E. (1946). An approximate distribution of estimates of variance components. *Biometrics Bulletin*, 2(6), 110-114.
"""
function stratifiedsampling(stratum::AbstractVector, value::AbstractVector, plotarea::Area,
  strataarea::AbstractDict{<:Any,<:Area}; allocation::Symbol=:proportional, e::Real=10, α::Real=0.95)

  length(stratum) == length(value) || throw(DimensionMismatch("stratum and value must have the same length."))
  allocation in (:proportional, :neyman) || throw(ArgumentError("allocation must be :proportional or :neyman."))

  areas = _resolvestrataareas(stratum, strataarea)
  table = _auxiliarytable(stratum, value, areas)
  all(>=(2), table.nh) || throw(ArgumentError("Every stratum needs at least two plots to estimate its variance."))
  anova = _anova(value, table)

  N = round(Int, ustrip(uconvert(NoUnits, sum(areas) / plotarea)))
  measured = table.nh
  n = sum(measured)
  infinite = _isinfinitepopulation(n, N)

  # the stratified mean, the variance of the population and the variance of the mean:
  # Eqs. 6.4, 6.8 and 6.12, or 6.13 when the sampling fraction is negligible
  x̅st = sum(table.Wh .* table.xmh)
  s²st = sum(table.Whs2h)
  s²x̅st = sum(table.Whsh .^ 2 ./ measured) - (infinite ? zero(s²st) : s²st / N)
  sx̅st = sqrt(s²x̅st)
  cv = sqrt(s²st) / x̅st * 100
  E = (e / 100) * x̅st

  # the sampling intensity of the proportional allocation (Eqs. 6.28 and 6.29) or of Neyman's
  # optimum allocation (Eqs. 6.34 and 6.35), and its distribution among the strata (Eqs. 6.27, 6.33)
  s²w = ustrip(s²st)
  ΣWhsh = ustrip(sum(table.Whsh))
  requiredtotal =
    allocation == :proportional ?
    _requiredsamplesize(s²w, n, ustrip(E), α, infinite ? nothing : N) :
    _requiredsamplesize(ΣWhsh^2, n, ustrip(E), α, infinite ? nothing : N; correctionvarterm=s²w)
  share = allocation == :proportional ? table.Wh : ustrip.(table.Whsh) ./ ΣWhsh
  required = ceil.(Int, requiredtotal .* share)
  missingplots = _missingplots(measured, required)

  # the degrees of freedom of Student's t are the effective ones of Satterthwaite (Eq. 6.26),
  # rounded up to the next integer as in the book
  Nh = ustrip.(uconvert.(NoUnits, areas ./ plotarea))
  gh = Nh .* (Nh .- measured) ./ measured
  dfeff = ceil(Int, sum(gh .* ustrip.(table.s2h))^2 / sum((gh .^ 2 .* ustrip.(table.s2h) .^ 2) ./ (measured .- 1)))
  Ttab = -quantile(TDist(dfeff), (1 - α) / 2)

  population = infinite ? "infinite" : "finite"
  absoluteerror = Ttab * sx̅st
  relativeerror = absoluteerror / x̅st * 100
  totalvalue = x̅st * N
  ciupper = totalvalue + N * absoluteerror
  cilower = totalvalue - N * absoluteerror

  resulttable = DataFrame(
    xm=x̅st, cv=round(cv, digits=2), s2xm=s²x̅st, sxm=sx̅st,
    Ea=absoluteerror, Er=round(relativeerror, digits=2), X=totalvalue,
    Xl=cilower, Xu=ciupper, population=population, f=round(n / N, digits=3),
    nh=Tuple(measured), nhreq=Tuple(required),
    nhadd=Tuple(missingplots), n=n, nreq=requiredtotal, N=N,
    Ah=Tuple(ustrip.(areas)),
  )

  return SamplingReport((; anova, auxiliarytable=table, resulttable=resulttable))
end
