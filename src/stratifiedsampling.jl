# Per-stratum descriptive statistics: sample size, mean, variance, standard deviation,
# and the stratum's relative weight (its share of the total population area) together
# with weighted standard deviation/variance, the building blocks of every stratified-
# sampling formula below. Generic to any number of strata via `groupby`.
function _auxiliarytable(stratum::AbstractVector, value::AbstractVector, areas::AbstractVector{<:Area})
  data = DataFrame(stratum=stratum, value=value)
  table = combine(groupby(data, :stratum, sort=true)) do df
    (n=nrow(df), mean=mean(df.value), variance=var(df.value), standardDeviation=std(df.value))
  end
  p = areas ./ sum(areas)
  table.areaShare = p
  table.weightedStandardDeviation = p .* table.standardDeviation
  table.weightedVariance = p .* table.variance
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
  ssb = sum(table.n .* abs2.(table.mean .- ȳ))
  ssw = sst - ssb
  dfb, dfw = k - 1, n - k
  msb, msw = ssb / dfb, ssw / dfw
  # a single stratum has no between-strata comparison to make: the F-test is undefined
  # (0 numerator degrees of freedom), so the between-strata row is left blank instead.
  fstat = dfb > 0 ? msb / msw : missing
  pval = dfb > 0 ? ccdf(FDist(dfb, dfw), fstat) : missing
  DataFrame(
    "Source" => ["Between strata", "Within strata", "Total"],
    "DOF" => [dfb, dfw, n - 1],
    "SS" => [ssb, ssw, sst],
    "MS" => [dfb > 0 ? msb : missing, msw, missing],
    "F" => [fstat, missing, missing],
    "p(>F)" => [pval, missing, missing],
  )
end

# Elementwise max(required - measured, 0): how many additional plots each stratum still
# needs, never negative when a stratum already exceeds its required allocation.
_missingplots(measured::AbstractVector{<:Integer}, required::AbstractVector{<:Integer}) = max.(required .- measured, 0)

# Resolves a stratum => area dict into a plain vector, in the ascending-sorted stratum
# order `groupby(..., sort=true)` will itself produce -- a real key-based lookup instead
# of trusting the caller to have supplied the areas in that order already, and errors
# clearly if any stratum present in the data has no matching entry in `strata_area`.
function _resolvestrataareas(stratum::AbstractVector, strata_area::AbstractDict{<:Any,<:Area})
  sortedstrata = sort(unique(stratum))
  missingstrata = setdiff(sortedstrata, keys(strata_area))
  isempty(missingstrata) ||
    throw(ArgumentError("strata_area has no entry for stratum(s) $(missingstrata)."))
  return [strata_area[s] for s in sortedstrata]
end

"""
    stratifiedSampling(stratum::AbstractVector, value::AbstractVector, plot_area::Area,
                        strata_area::AbstractDict{<:Any,<:Area}; e::Real=10, α::Real=0.95)
    stratifiedSampling(stratum::AbstractVector, value::AbstractVector, plot_area::Real,
                        strata_area::AbstractDict{<:Any,<:Real}; kwargs...)

Performs stratified random sampling for forest inventory analysis, generic to any number of strata.

# Description

Stratified random sampling divides the population into non-overlapping subdivisions
(strata) that are each more internally homogeneous than the population as a whole, then
samples independently within each. This usually shrinks the variance of the estimated
mean relative to simple random sampling of the same total size, at the cost of needing
each stratum's area up front.

# Arguments

- `stratum::AbstractVector`: each plot's stratum label, same order as `value`.
- `value::AbstractVector`: measurement for each plot, same order as `stratum`. If given as a `Unitful` quantity, its unit is preserved; plain numbers are kept unitless throughout.
- `plot_area`: the area of each sample plot, as an `Area` quantity. A plain number is taken to be hectares.
- `strata_area::AbstractDict`: each stratum's total area, keyed by the same labels used in `stratum` — e.g. `Dict(1 => 12.0u"ha", 2 => 8.0u"ha", 3 => 20.0u"ha")`. A dict of plain numbers is taken to be hectares. Using a dictionary rather than a plainly-ordered vector means the plots in `stratum`/`value` don't need to be pre-sorted by stratum — they can be given in whatever order they were collected in the field.
- `e::Real=10`: desired relative error margin as a percentage (default 10%).
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- [`SamplingReport`](@ref) with three tables: `anova` (test for a difference between
  strata means), `auxiliaryTable` (per-stratum descriptive statistics and allocation
  weights), and `resultTable` (the final stratified estimates, one row, one column per
  statistic — same shape as [`simpleCasualSampling`](@ref)'s return value). `resultTable`
  columns:
  - `mean`, `coefficientOfVariation`, `varianceOfMean`, `standardError`, `absoluteError`, `relativeError`, `total`, `lowerLimit`, `upperLimit`, `population`, `samplingFraction`: as in [`simpleCasualSampling`](@ref), using the stratified mean/variance.
  - `sampleSizeByStratum`, `requiredSampleSizeByStratum`, `additionalSampleSizeByStratum`: per-stratum measured/required/missing plot counts, as tuples in stratum order.
  - `n`, `requiredSampleSize`: total measured and required plot counts across all strata.
  - `N`: total number of possible plots across all strata.
  - `areaByStratum`: each stratum's area, as a tuple of plain numbers in the same unit `strata_area`
    was supplied in (hectares by default). Unlike every other column here, this one is
    **not** `Unitful`-tagged and is left untouched by `removeunits`/`restoreunits` — a
    `Tuple`-typed column has no single `eltype` those functions can strip/restore a unit
    from, so the value is stored already-unitless instead of silently carrying a hidden
    `Unitful` type through the round trip.

# Mathematical basis

The stratified mean and its variance, weighting each stratum by its area share `pₕ = Aₕ/ΣAₕ`:
```math
\\bar{x}_{st} = \\sum_h p_h \\bar{x}_h \\qquad
s^2_{\\bar{x}_{st}} = \\sum_h \\frac{(p_h s_h)^2}{n_h} - \\frac{\\sum_h p_h s_h^2}{N}
```

Optimal (Neyman) allocation distributes the required total sample size across strata in
proportion to `pₕsₕ`, solved iteratively the same way as [`simpleCasualSampling`](@ref)
(see `_requiredsamplesize`). The confidence interval uses a Satterthwaite-approximated
effective degrees of freedom that accounts for the strata's differing sizes and variances
(Cochran, 1977, eq. 5.35):
```math
df = \\frac{\\left(\\sum_h g_h s_h^2\\right)^2}{\\sum_h \\dfrac{g_h^2 s_h^4}{n_h - 1}}
\\qquad \\text{where} \\qquad g_h = \\frac{N_h(N_h - n_h)}{n_h}
```

# Technical description

With a single stratum this reduces exactly to [`simpleCasualSampling`](@ref) — the
`anova` table's between-strata line becomes meaningless (0 degrees of freedom) and the
weighted mean collapses to the plain sample mean. The `auxiliaryTable`'s
`weightedStandardDeviation`/`weightedVariance` columns are the building blocks of every
downstream formula (allocation, variance, required sample size) and are exposed directly
so they can be audited.

# Examples

```julia-repl
julia> stratum = [1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3];

julia> value = [18.2, 21.4, 19.8, 20.1, 32.5, 35.1, 30.8, 12.4, 11.9, 13.6, 12.8, 13.1]u"m^3";

julia> report = stratifiedSampling(stratum, value, 0.1, Dict(1 => 12.0, 2 => 8.0, 3 => 20.0));

julia> resultTable(report)
julia> auxiliaryTable(report)
julia> anova(report)
```
"""
function stratifiedSampling(stratum::AbstractVector, value::AbstractVector, plot_area::Area,
  strata_area::AbstractDict{<:Any,<:Area}; e::Real=10, α::Real=0.95)

  length(stratum) == length(value) || throw(DimensionMismatch("stratum and value must have the same length."))

  areas = _resolvestrataareas(stratum, strata_area)
  table = _auxiliarytable(stratum, value, areas)
  anova = _anova(value, table)

  N = round(Int, ustrip(uconvert(NoUnits, sum(areas) / plot_area)))
  x̅st = sum(table.areaShare .* table.mean)
  s²x̅st = sum(table.weightedStandardDeviation .^ 2 ./ table.n) - sum(table.weightedVariance) / N
  sx̅st = sqrt(s²x̅st)
  cv = sx̅st / x̅st * 100
  E = (e / 100) * x̅st
  measured = table.n
  n = sum(measured)

  infinite = _isinfinitepopulation(n, N)
  requiredtotal = _requiredsamplesize(ustrip(sum(table.weightedVariance)), n, ustrip(E), α, infinite ? nothing : N)
  required = round.(Int, requiredtotal .* table.weightedStandardDeviation ./ sum(table.weightedStandardDeviation))
  missingplots = _missingplots(measured, required)

  Nh = ustrip.(uconvert.(NoUnits, areas ./ plot_area))
  gh = Nh .* (Nh .- measured) ./ measured
  dfeff = sum(gh .* table.variance)^2 / sum((gh .^ 2 .* table.variance .^ 2) ./ (measured .- 1))
  Ttab = -quantile(TDist(dfeff), (1 - α) / 2)

  population = infinite ? "infinite" : "finite"
  f = 1 - n / N
  absoluteerror = Ttab * sx̅st
  relativeerror = absoluteerror / x̅st * 100
  totalvalue = x̅st * N
  ciupper = totalvalue + N * absoluteerror
  cilower = totalvalue - N * absoluteerror

  resulttable = DataFrame(
    mean=x̅st, coefficientOfVariation=round(cv, digits=2), varianceOfMean=s²x̅st, standardError=sx̅st,
    absoluteError=absoluteerror, relativeError=round(relativeerror, digits=2), total=totalvalue,
    lowerLimit=cilower, upperLimit=ciupper, population=population, samplingFraction=round(f, digits=3),
    sampleSizeByStratum=Tuple(measured), requiredSampleSizeByStratum=Tuple(required),
    additionalSampleSizeByStratum=Tuple(missingplots), n=n, requiredSampleSize=requiredtotal, N=N,
    areaByStratum=Tuple(ustrip.(areas)),
  )

  return SamplingReport((; anova, auxiliaryTable=table, resultTable=resulttable))
end

function stratifiedSampling(stratum::AbstractVector, value::AbstractVector, plot_area::Real,
  strata_area::AbstractDict{<:Any,<:Real}; kwargs...)
  stratifiedSampling(stratum, value, plot_area * AUNIT, Dict(k => v * AUNIT for (k, v) in strata_area); kwargs...)
end
