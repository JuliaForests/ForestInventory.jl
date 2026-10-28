# Internal helpers shared by every sampling design in this submodule. None of these are
# exported; each public function calls into them with the specific quantities its own
# variance formula produces.

# Required sample size for an infinite population, given a t-value, a *squared*
# variance-like term (`varterm`, e.g. `cv^2` or a variance in volume² units -- not `cv` or
# a standard deviation) and an admissible error (`e`) on the matching linear scale (both
# relative, e.g. percent of the mean, or both absolute, e.g. volume units).
_infinitesamplesize(t::Real, varterm::Real, e::Real) = (t^2 * varterm) / e^2

# Required sample size for a finite population of size npop, on the same scale as above.
# `correctionvarterm` lets the finite-population correction term use a different
# variance-like quantity than the numerator -- needed by `twostagesampling`, where the
# correction reflects the population size of secondary units while the numerator reflects
# how many were actually measured; it defaults to `varterm` for every other design, where
# the two coincide.
_finitesamplesize(t::Real, varterm::Real, e::Real, npop::Real, correctionvarterm::Real=varterm) =
  (t^2 * varterm) / (e^2 + (t^2 * correctionvarterm) / npop)

# Iteratively solves for the required sample size. The t quantile depends on the degrees
# of freedom of the very sample size it is used to compute, so the pair is refined
# together until the estimate stabilizes (matching the fixed-point iteration used
# throughout the legacy ForestInventory.jl sampling routines). `npop=nothing` treats the
# population as infinite; passing a population size switches to the finite-population
# formula. Shared by every sampling design below instead of each reimplementing the loop.
function _requiredsamplesize(varterm::Real, n::Real, e::Real, α::Real, npop::Union{Real,Nothing}=nothing;
  correctionvarterm::Real=varterm, maxiterations::Int=1000)
  samplesize(t::Real) = isnothing(npop) ? _infinitesamplesize(t, varterm, e) : _finitesamplesize(t, varterm, e, npop, correctionvarterm)
  for _ in 1:maxiterations
    t1 = -quantile(TDist(n - 1), (1 - α) / 2)
    # a near-zero variance (exactly the outcome a well-executed design aims for) can
    # solve the formula down to under 1 "required" observation; TDist needs at least 1
    # degree of freedom, and no valid design needs fewer than 2 observations anyway, so
    # the candidate size is floored at 2 before it feeds the next iteration's quantile.
    iter1 = max(samplesize(t1), 2.0)
    t2 = -quantile(TDist(iter1 - 1), (1 - α) / 2)
    iter2 = max(samplesize(t2), 2.0)
    if round(n, digits=3) == round(iter1, digits=3) || round(n, digits=3) == round(iter2, digits=3)
      return ceil(Int, n)
    end
    n = iter1
  end
  error("Maximum number of iterations ($maxiterations) reached without convergence.")
end

# Whether a sample of size `n` drawn without replacement from a population of size `npop`
# can be treated as effectively infinite: the standard forest-mensuration rule of thumb
# is that the finite-population correction is negligible once the sampling fraction stays
# under 2%.
_isinfinitepopulation(n::Real, npop::Real) = (1 - n / npop) >= 0.98

# Between/within primary-unit variance decomposition shared by every design built on
# equal-size primary sampling units (clusters, systematic starts, or first-stage units of
# a two-stage design): `means`/`vars` are the per-primary-unit mean and internal variance
# (one entry per sampled primary unit), `mpop` is the number of secondary units per primary
# unit (assumed equal), and `npop` is the number of possible primary units in the population.
#
# Returns a named tuple with the pooled within-unit variance, the between-unit variance
# (per secondary unit), their sum (the total per-secondary-unit variance), the intraclass
# correlation coefficient, and the variance of the overall mean.
function _clustervariance(means::AbstractVector{<:Real}, vars::AbstractVector{<:Real}, mpop::Real, npop::Real)
  n = length(means)
  ȳ = mean(means)
  # a cluster of exactly one secondary unit has no internal spread to estimate (`var` of a
  # single observation is NaN, not 0); with mpop = 1 there simply is no within-unit variance,
  # and this reduces the whole decomposition to simple random sampling of the means.
  withinvar = mpop == 1 ? zero(eltype(vars)) : mean(vars)  # pooled within-unit variance: Σ(M-1)sⱼ² / (n(M-1)) simplifies to mean(sⱼ²) for equal M
  msb = mpop * sum(abs2, means .- ȳ) / (n - 1)  # between-unit mean square
  betweenvar = (msb - withinvar) / mpop
  totalvar = betweenvar + withinvar
  icc = totalvar > 0 ? betweenvar / totalvar : zero(betweenvar)
  # the finite-population correction multiplies only the between-clusters term (Péllico Netto &
  # Brena, Eq. 9.17); it is dropped, Eq. 9.18, when the sampling fraction is negligible
  f = _isinfinitepopulation(n, npop) ? one(npop) : (npop - n) / npop
  meanvar = f * betweenvar / n + withinvar / (n * mpop)
  varterm = betweenvar + withinvar / mpop  # feeds `_requiredsamplesize` directly
  return (; withinvar, betweenvar, totalvar, icc, meanvar, varterm)
end

# Simple-random-sampling point estimate and variance of the mean for one occasion of a
# successive-occasions design (independent/complete-replacement/partial-replacement/
# double sampling), reused instead of recomputing it for every occasion in every design.
function _occasionmean(value::AbstractVector, npop::Real)
  n = length(value)
  x̅, s² = mean_and_var(value)
  meanvar = (s² / n) * (1 - n / npop)
  return (; n, x̅, s², meanvar)
end

# Builds one occasion's result table (mean, standard error, confidence interval, total)
# as a single-row wide `DataFrame`, shared by every successive-occasions design below.
# Column names are full words and unit-homogeneous per column -- exactly the shape
# `stripunits`/`parseunits` (ForestFoundations.jl) expect.
function _occasiontable(stats::NamedTuple, npop::Real, plotarea::Area, α::Real)
  t = -quantile(TDist(stats.n - 1), (1 - α) / 2)
  se = sqrt(stats.meanvar)
  abserr = t * se
  relerr = abserr / stats.x̅ * 100
  total = stats.x̅ * npop
  # expansionfactor expresses the plot area in its reference area (ha/ac): dividing straight
  # by a bare m^2 quantity would let Unitful silently cancel same-base powers (m^3/m^2 ->
  # m^1) instead of keeping the m^3/ha shape callers expect.
  perhectare = stats.x̅ * expansionfactor(plotarea)
  return DataFrame(
    xm=stats.x̅, sxm=se, Ea=abserr, Er=round(relerr, digits=2),
    xha=perhectare, X=total, Xl=total - npop * abserr, Xu=total + npop * abserr,
    n=stats.n, N=npop,
  )
end

# Builds the change-between-occasions result table (mean growth, standard error,
# confidence interval, total growth) as a single-row wide `DataFrame`, shared by every
# successive-occasions design below.
function _changetable(growth, growthvar::Real, df::Real, npop::Real, α::Real)
  t = -quantile(TDist(df), (1 - α) / 2)
  se = sqrt(growthvar) * unit(growth)
  abserr = t * se
  relerr = abserr / growth * 100
  total = growth * npop
  return DataFrame(
    dm=growth, sdm=se, Ea=abserr, Er=round(relerr, digits=2),
    D=total, Dl=total - npop * abserr, Du=total + npop * abserr,
  )
end

# Ordinary least squares slope of `y` on `x` via the normal equations, solved through the
# same transposed-design-matrix approach as the ANOVA helper in `stratifiedsampling.jl`
# (`X'X`, `X'y`, Cholesky solve) rather than a dedicated regression model — this module
# only ever needs the slope of a single predictor, used as the regression coefficient in
# the double-sampling and partial-replacement estimators.
function _olsslope(y::AbstractVector{<:Real}, x::AbstractVector{<:Real})
  checklengths(:y => y, :x => x)
  X = [ones(length(x)) x]
  β = X'y
  chol = cholesky!(Symmetric(X'X))
  ldiv!(chol, β)
  return β[2]
end
