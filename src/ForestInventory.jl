"""
Description:

ForestInventory.jl is a Julia package for forest-inventory sampling: estimating stand-level
volume (or any other measured variable) and its confidence interval from a set of sample
plots or points. A single entry point, [`sampling`](@ref), dispatches on a
[`SamplingDesign`](@ref) type to select among 11 classic designs — simple
([`SimpleCasualSampling`](@ref)), stratified ([`StratifiedSampling`](@ref)), systematic
([`SystematicSampling`](@ref), [`MultistartSystematicSampling`](@ref)), cluster/two-stage
([`ClusterSampling`](@ref), [`TwoStageSampling`](@ref)), horizontal point/Bitterlich
angle-count sampling ([`HorizontalPointSampling`](@ref)), and the four
sampling-on-successive-occasions designs used in continuous forest inventory
([`IndependentOccasionsSampling`](@ref), [`CompleteReplacementSampling`](@ref),
[`PartialReplacementSampling`](@ref), [`DoubleSampling`](@ref)) — every one generic to any
number of strata/clusters/plots and returning a [`SamplingReport`](@ref) or plain
`DataFrame`. `subtypes(SamplingDesign)` lists all 11 at once.

# Units of measurement

Every function is built on `Unitful` quantities, so measurements carry their unit and are
converted automatically. Diameters may be given in `cm`, `mm`, `inch`, ..., and areas/volumes
in any compatible unit; results are normalised to a single coherent unit.

For convenience every function also accepts plain numbers, in which case the package assumes
its default units:

| Measurement                  | Assumed unit |
|:------------------------------|:-------------|
| diameter                      | `cm`         |
| plot/point area, total area   | `ha`         |
| volume                        | `m^3`        |

So `sampling(SimpleCasualSampling, volume, 0.05, 10)` reads `volume` as `m^3`, the plot area
as `0.05 ha`, and the total area as `10 ha`, returning exactly the same result as the
equivalent call with explicit units.

Results always carry units; use `removeunits`/`restoreunits` from
[ForestFoundations.jl](https://github.com/JuliaForests/ForestFoundations.jl) to move between
unitful DataFrames and plain numeric ones for export.
"""
module ForestInventory

using DataFrames, Distributions, ForestFoundations, LinearAlgebra, StatsBase

include("common.jl")
include("designs.jl")
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
  sampling,
  SamplingDesign,
  SimpleCasualSampling,
  StratifiedSampling,
  SystematicSampling,
  MultistartSystematicSampling,
  ClusterSampling,
  HorizontalPointSampling,
  TwoStageSampling,
  IndependentOccasionsSampling,
  CompleteReplacementSampling,
  PartialReplacementSampling,
  DoubleSampling,
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
