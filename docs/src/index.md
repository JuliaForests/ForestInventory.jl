```@meta
CurrentModule = ForestInventory
```

# ForestInventory

Documentation for [ForestInventory](https://github.com/JuliaForests/ForestInventory.jl) —
the forest-inventory sampling engine behind
[ForestMensuration.jl](https://github.com/JuliaForests/ForestMensuration.jl).

## Installation

```julia-repl
using Pkg
Pkg.add("ForestInventory")
```

## Quickstart

```@example quickstart
using ForestInventory

v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1]

sampling(SimpleCasualSampling, v, 0.05, 10; e=10, α=0.95)
```

All 11 classic forest inventory sampling designs go through a single entry point,
[`sampling`](@ref): pass a [`SamplingDesign`](@ref) subtype as the first argument to select
the design, then that design's own arguments — `subtypes(SamplingDesign)` lists all 11 at
once. Every design accepts plain numbers (volume defaults to `m^3`, plot/total areas to
`ha`) or explicit `Unitful` quantities. Designs that produce a single table return a plain
`DataFrame`; designs that produce several related tables (one per stratum/cluster, or one
per inventory occasion) return a [`SamplingReport`](@ref), whose tables are reachable
through accessor functions (`resultTable(report)`, `clusterTable(report)`, ...) or directly
as properties (`report.resultTable`).

The examples below walk through every design.

## Simple Random Sampling

[`SimpleCasualSampling`](@ref) is the reference design every other method in this section
is compared against: each plot has an equal chance of being selected, with no further
structure. See the [Quickstart](#Quickstart) example above.

## Stratified Random Sampling

[`StratifiedSampling`](@ref) divides the population into non-overlapping strata before
sampling within each — usually a variance reduction over simple random sampling of the
same total size. It returns a [`SamplingReport`](@ref) with an ANOVA table (is there
really a difference between strata?), an auxiliary table (per-stratum descriptive
statistics and allocation weights), and the final result table.

```@example inv_stratified
using ForestInventory, DataFrames

data = DataFrame(
    stratum=[1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3],
    volume=[18.2, 21.4, 19.8, 20.1, 32.5, 35.1, 30.8, 12.4, 11.9, 13.6, 12.8, 13.1],
)

# strata_area is given in the same order as the sorted strata: 1, 2, 3
report = sampling(StratifiedSampling, :stratum, :volume, 0.1, [12.0, 8.0, 20.0], data)
resultTable(report)
```

```@example inv_stratified
auxiliaryTable(report)
```

```@example inv_stratified
anova(report)
```

With a single stratum, [`StratifiedSampling`](@ref) reduces exactly to
[`SimpleCasualSampling`](@ref) — useful as a sanity check when strata are added
incrementally to a growing dataset.

## Systematic Sampling

[`SystematicSampling`](@ref) estimates the variance of the mean from the method of
successive differences between consecutive plots, since plots laid out at a fixed
interval tend to be more alike than a true random sample.

```@example inv_systematic
using ForestInventory

v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1]

sampling(SystematicSampling, v, 0.05, 10)
```

An optional `line` vector groups plots into several independent transects, so the
difference between the last plot of one line and the first plot of the next — which
aren't actually adjacent on the ground — is excluded from the variance estimate:

```@example inv_systematic
line = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
sampling(SystematicSampling, v, 0.05, 10; line=line)
```

## Systematic Sampling with Multiple Random Starts

[`MultistartSystematicSampling`](@ref) lays out several independent systematic lines,
each starting at its own random point, turning the design into a genuine probability
sample: statistically identical to [`ClusterSampling`](@ref)'s between/within
decomposition, with each line playing the role of one cluster.

```@example inv_multistart
using ForestInventory, DataFrames

data = DataFrame(
    start=repeat(1:4, inner=5),
    volume=[18.2, 19.1, 17.8, 18.9, 19.4, 22.4, 23.1, 21.9, 22.8, 23.4,
            15.1, 16.0, 14.8, 15.6, 15.9, 20.5, 21.2, 19.8, 20.6, 21.0],
)

report = sampling(MultistartSystematicSampling, :start, :volume, 0.02, 15, data)
resultTable(report)
```

## One-Stage Cluster Sampling

[`ClusterSampling`](@ref) samples clusters of `M` neighboring plots ("conglomerates")
instead of individual plots — cheaper to lay out in the field, at the cost of
within-cluster homogeneity typically inflating the variance of the mean.

```@example inv_cluster
using ForestInventory, DataFrames

data = DataFrame(
    cluster=repeat(1:6, inner=4),
    volume=[18.2, 19.1, 17.8, 18.9, 22.4, 23.1, 21.9, 22.8, 15.1, 16.0, 14.8, 15.6,
            27.3, 28.1, 26.9, 27.8, 19.8, 20.5, 19.1, 20.0, 24.5, 25.2, 23.9, 24.8],
)

report = sampling(ClusterSampling, :cluster, :volume, 0.02, 15, data)
resultTable(report)
```

```@example inv_cluster
clusterTable(report)
```

## Horizontal Point Sampling (Bitterlich)

[`HorizontalPointSampling`](@ref) estimates volume and basal area per hectare from
angle-count ("Bitterlich") data: at each point, every tree that is "in" for the chosen
basal area factor (`baf`) represents exactly `baf` m²/ha of basal area, regardless of its
own size — no fixed plot radius is needed at all. Point 6 below was visited but had no
"in" trees; it is not a row of `data`, but is still counted through `npoints=6` so the
zero observation is not silently dropped from the mean.

```@example inv_hps
using ForestInventory, DataFrames

data = DataFrame(
    point=[1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 5, 5, 5],
    diameter=[25.0, 30.0, 20.0, 28.0, 22.0, 35.0, 30.0, 25.0, 20.0, 22.0, 30.0, 28.0, 26.0],
    volume=[0.35, 0.55, 0.22, 0.48, 0.28, 0.85, 0.55, 0.35, 0.22, 0.28, 0.55, 0.48, 0.40],
)

# baf=2 m²/ha, 6 points visited (5 counted trees + 1 empty), 0.1 ha per point, 10 ha stand
report = sampling(HorizontalPointSampling, :point, :diameter, :volume, 2.0, 6, 0.1, 10, data)
resultTable(report)
```

```@example inv_hps
pointTable(report)
```

## Two-Stage Sampling

[`TwoStageSampling`](@ref) draws `n` primary units (e.g. stands) from a population of `N`,
then sub-samples `m` secondary units (plots) from within each drawn primary out of `M`
possible — unlike cluster sampling, not every secondary unit inside a drawn primary needs
to be measured.

```@example inv_twostage
using ForestInventory, DataFrames

data = DataFrame(
    primary=repeat(1:5, inner=3),
    volume=[18.2, 19.1, 17.8, 22.4, 23.1, 21.9, 15.1, 16.0, 14.8, 27.3, 28.1, 26.9, 19.8, 20.5, 19.1],
)

# N=40 possible primary units, M=6 possible secondary units per primary
report = sampling(TwoStageSampling, :primary, :volume, 0.02, 40, 6, data)
resultTable(report)
```

## Successive Occasions

The remaining four designs estimate growth between two inventory occasions, differing in
how much of the plot network is remeasured at the second occasion. All four return a
[`SamplingReport`](@ref) with `occasion1`, `occasion2`, and `change` tables.

### Independent Samples

[`IndependentOccasionsSampling`](@ref) samples each occasion completely independently —
the simplest design, but the least efficient at detecting growth, since it exploits none
of the natural plot-to-plot correlation between occasions.

```@example inv_independent
using ForestInventory

v1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0]
v2 = [24.1, 23.8, 22.9, 25.6, 24.0]

report = sampling(IndependentOccasionsSampling, v1, v2, 0.05, 200, 200)
change(report)
```

### Complete Replacement

[`CompleteReplacementSampling`](@ref) remeasures the exact same plots at both occasions,
exploiting their positive correlation to tighten the growth estimate.

```@example inv_complete
using ForestInventory

v1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0, 24.6, 17.3]
v2 = [22.1, 25.8, 23.4, 21.0, 26.6, 20.9, 26.0, 22.5]

report = sampling(CompleteReplacementSampling, v1, v2, 0.05, 200, 200)
change(report)
```

### Partial Replacement

[`PartialReplacementSampling`](@ref) is the middle ground: a matched subset of plots is
remeasured, some temporary plots are dropped, and new temporary plots are added at the
second occasion. `missing` marks a plot that wasn't measured on a given occasion.

```@example inv_partial
using ForestInventory

volume1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0, 23.1, missing, missing]
volume2 = [missing, missing, 23.4, 24.0, 26.6, 22.9, 27.5, 25.2, 21.8]

report = sampling(PartialReplacementSampling, volume1, volume2, 0.05, 200)
occasion2(report)
```

```@example inv_partial
change(report)
```

### Double Sampling

[`DoubleSampling`](@ref) measures a large first-occasion sample but only remeasures a
smaller "permanent" subset at the second occasion, estimating the rest via a regression of
the permanent subset's second-occasion volume on its first-occasion volume — solved
through the normal equations `(X'X)β = X'y`, the same transposed-design-matrix approach
used for the ANOVA in [`StratifiedSampling`](@ref).

```@example inv_double
using ForestInventory

volume1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0, 23.1, 17.6, 20.8, 21.9]
volume2 = [22.1, 25.8, 23.4, missing, 26.6, missing, 27.5, missing, missing, 26.0]

report = sampling(DoubleSampling, volume1, volume2, 0.05, 200)
occasion2(report)
```

```@example inv_double
change(report)
```

## Index

```@index
```

## API Reference

```@autodocs
Modules = [ForestInventory]
```
