"""
Description:

ForestInventory.jl facilitates the processing of forest-inventory sampling procedures: it
estimates a stand-level total (and its confidence interval) — volume, basal area, tree
count, or any other measured quantity — from a set of sample plots or points, following the
sampling designs described in Péllico Netto & Brena's *Inventário Florestal* (1997). It
provides 11 classic designs, each its own function: simple ([`simpleCasualSampling`](@ref)),
stratified ([`stratifiedSampling`](@ref)), systematic ([`systematicSampling`](@ref),
[`multistartSystematicSampling`](@ref)), cluster/two-stage ([`clusterSampling`](@ref),
[`twoStageSampling`](@ref)), horizontal point/Bitterlich angle-count sampling
([`horizontalPointSampling`](@ref)), and the four sampling-on-successive-occasions designs
used in continuous forest inventory ([`independentOccasionsSampling`](@ref),
[`completeReplacementSampling`](@ref), [`partialReplacementSampling`](@ref),
[`doubleSampling`](@ref)) — every one generic to any number of strata/clusters/plots and
returning a [`SamplingReport`](@ref) or plain `DataFrame`.

# Units of measurement

Every function is built on `Unitful` quantities. Diameters and areas may be given in `cm`,
`ha`, `acre`, ... or any compatible unit; plain numbers for these default to `cm`/`ha`
respectively, matching the rest of the `JuliaForests` ecosystem. The **sampled quantity
itself** (usually named `value`, `value1`/`value2`) is the one deliberate exception: since
it might be volume, basal area, a tree count, or anything else being sampled, this package
never assumes a unit for it — a `Unitful` quantity is preserved and propagated as given,
and a plain number stays a plain, unitless number all the way through to the result table,
rather than being silently treated as `m^3` the way it would be elsewhere in
`JuliaForests`.

Results carry whatever unit the input did (or none); use `removeunits`/`restoreunits` from
[ForestFoundations.jl](https://github.com/JuliaForests/ForestFoundations.jl) to move between
unitful DataFrames and plain numeric ones for export.
"""
module ForestInventory

using DataFrames, Distributions, ForestFoundations, LinearAlgebra, StatsBase

include("common.jl")
include("inventoryreport.jl")
include("simplecasualsampling.jl")
include("stratifiedsampling.jl")
include("systematicsampling.jl")
include("clustersampling.jl")
include("horizontalpointsampling.jl")
include("multistartsystematicsampling.jl")
include("twostagesampling.jl")
include("independentoccasionssampling.jl")
include("completereplacementsampling.jl")
include("partialreplacementsampling.jl")
include("doublesampling.jl")

export
  simpleCasualSampling,
  stratifiedSampling,
  systematicSampling,
  multistartSystematicSampling,
  clusterSampling,
  horizontalPointSampling,
  twoStageSampling,
  independentOccasionsSampling,
  completeReplacementSampling,
  partialReplacementSampling,
  doubleSampling,
  SamplingReport,
  resultTable,
  clusterTable,
  auxiliaryTable,
  pointTable,
  startTable,
  primaryTable,
  anova,
  occasion1,
  occasion2,
  change

end
