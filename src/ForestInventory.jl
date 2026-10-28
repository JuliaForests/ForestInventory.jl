"""
Description:

ForestInventory.jl facilitates the processing of forest-inventory sampling procedures: it
estimates a stand-level total (and its confidence interval) — volume, basal area, tree
count, or any other measured quantity — from a set of sample plots or points, following the
sampling designs described in Péllico Netto & Brena's *Inventário Florestal* (1997). It
provides 11 classic designs, each its own function: simple ([`simplecasualsampling`](@ref)),
stratified ([`stratifiedsampling`](@ref)), systematic ([`systematicsampling`](@ref),
[`multistartsystematicsampling`](@ref)), cluster/two-stage ([`clustersampling`](@ref),
[`twostagesampling`](@ref)), horizontal point/Bitterlich angle-count sampling
([`horizontalpointsampling`](@ref)), and the four sampling-on-successive-occasions designs
used in continuous forest inventory ([`independentoccasionssampling`](@ref),
[`completereplacementsampling`](@ref), [`partialreplacementsampling`](@ref),
[`doublesampling`](@ref)) — every one generic to any number of strata/clusters/plots and
returning a [`SamplingReport`](@ref). Three more designs, due to Strand, Prodan, and Grosenbaugh rather than Péllico Netto &
Brena, round out the package: horizontal line sampling ([`strandsampling`](@ref),
Strand, 1957) is a proper PPS design like [`horizontalpointsampling`](@ref) and gets the
same full confidence-interval treatment; Prodan's k-tree nearest-neighbor sampling
([`prodansampling`](@ref), Prodan, 1968b) is a distance-based plotless method with no
exactly-unbiased variance theory in the literature, so it reports descriptive spread
instead; and 3P sampling ([`threepsampling`](@ref), Grosenbaugh, 1965) selects trees for
exact measurement with probability proportional to a quick field prediction, correcting
the bias with a ratio estimator over a 100%-cruised population (see each one's own
docstring for the caveats specific to it).

# Notation

The symbols follow Péllico Netto & Brena (*Inventário Florestal*, 1997), in ASCII and in the form of the other packages of the ecosystem: lowercase for the estimates of the sample, uppercase for the totals and the values of the population, the index written next to the symbol (`xmh`, the mean ``\\bar{x}_h`` of the stratum ``h``), and `m` for a mean, `l` and `u` for the lower and upper limits. The columns of the tables of every design use them:

| Column | Quantity |
|:-------|:---------|
| `n`, `N` | sampled units, possible units of the population |
| `f` | sampling fraction, ``f = n/N`` (the finite-population correction is ``1 - f``) |
| `xm`, `xha` | mean ``\\bar{x}`` per sampling unit, mean per hectare (or acre) |
| `cv` | coefficient of variation, ``cv = s/\\bar{x} \\times 100`` (%) |
| `s2xm`, `sxm` | variance of the mean ``s^2_{\\bar{x}}``, standard error ``s_{\\bar{x}}`` |
| `Ea`, `Er` | absolute sampling error ``E_a = t s_{\\bar{x}}``, relative error ``E_r`` (%) |
| `X`, `Xl`, `Xu` | estimated total ``\\hat{X} = N \\bar{x}`` and the limits of its confidence interval |
| `nreq`, `nadd` | sampling intensity required for the admissible error `e`, and the units still to measure |
| `A` | area of the forest |
| `s2e`, `s2d`, `r` | variance between and within clusters (`e` = *entre*, `d` = *dentro*), intracluster correlation |
| `nh`, `xmh`, `s2h`, `sh`, `Wh`, `Ah` | sample size, mean, variance, standard deviation, weight ``W_h = A_h/A`` and area of the stratum ``h`` |
| `dm`, `D` | mean change per unit and total change between two occasions |
| `G`, `Gt` | basal area per hectare, total basal area |

Every sampling error uses Student's ``t``, and the number of units that an admissible error `e` (in percent of the mean) requires is solved iteratively, since ``t`` depends on the degrees of freedom of that very number.

# Units of measurement

Every function is built on `Unitful` quantities. Diameters, distances and areas must carry
units and may be given in `cm`, `ha`, `acre`, ... or any compatible unit: a plain number
says nothing about its unit, so none is assumed. The **sampled quantity
itself** (usually named `value`, `value1`/`value2`) is the one deliberate exception: since
it might be volume, basal area, a tree count, or anything else being sampled, this package
never assumes a unit for it — a `Unitful` quantity is preserved and propagated as given,
and a plain number stays a plain, unitless number all the way through to the result table,
rather than being treated as any particular unit.

Results carry whatever unit the input did (or none); use `stripunits`/`parseunits` from
[ForestFoundations.jl](https://github.com/JuliaForests/ForestFoundations.jl) to move between
unitful DataFrames and plain numeric ones for export.
"""
module ForestInventory

using DataFrames, Distributions, LinearAlgebra, Reexport, StatsBase

# the units (`@u_str`, `ustrip`, `uconvert`, ...) and the helpers of ForestFoundations, so
# `using ForestInventory` is enough to write quantities
@reexport using ForestFoundations

include("common.jl")
include("inventoryreport.jl")
include("simplecasualsampling.jl")
include("stratifiedsampling.jl")
include("systematicsampling.jl")
include("clustersampling.jl")
include("horizontalpointsampling.jl")
include("prodansampling.jl")
include("strandsampling.jl")
include("threepsampling.jl")
include("multistartsystematicsampling.jl")
include("twostagesampling.jl")
include("independentoccasionssampling.jl")
include("completereplacementsampling.jl")
include("partialreplacementsampling.jl")
include("doublesampling.jl")

export
  simplecasualsampling,
  stratifiedsampling,
  systematicsampling,
  multistartsystematicsampling,
  clustersampling,
  horizontalpointsampling,
  prodansampling,
  strandsampling,
  threepsampling,
  twostagesampling,
  independentoccasionssampling,
  completereplacementsampling,
  partialreplacementsampling,
  doublesampling,
  SamplingReport,
  resulttable,
  clustertable,
  auxiliarytable,
  pointtable,
  starttable,
  primarytable,
  anova,
  occasion1,
  occasion2,
  change

end
