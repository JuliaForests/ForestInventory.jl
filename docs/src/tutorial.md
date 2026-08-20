# Getting Started

Estimate a stand's total — with its confidence interval, coefficient of variation, and
required sample size — from any of the 11 classic sampling designs, each its own function
taking plain vectors. Every design accepts `plot_area`/`total_area` as a `Unitful` area
(`ha`, `acre`, `m^2`, ...) or, for either one given as a plain number, assumes it to be
`ha`. This never extends to the sampled quantity itself (usually named `value`): it might
be volume, basal area, a tree count, or anything else, so a `Unitful` quantity is preserved
as given and a plain number stays unitless all the way through. Designs that produce a single table return a plain `DataFrame`; designs that
produce several related tables (one per stratum/cluster, or one per inventory occasion)
return a [`SamplingReport`](@ref), whose tables are reachable through accessor functions
(`resultTable(report)`, `clusterTable(report)`, ...) or directly as properties
(`report.resultTable`).

The examples below walk through every design, from the simplest — simple random sampling —
to horizontal point (Bitterlich) sampling and beyond. Most of them (every design except
horizontal point sampling, which has no equivalent in the original 10) reuse the exact
plot-level measurements from the classic exercises in Nardini's dissertation — the
Master's work this package's lineage traces back to (see [About the
Authors](@ref)) — pulled from the original desktop tool's own repository
([claitonnardini/Inventario-Florestal.jl](https://github.com/claitonnardini/Inventario-Florestal.jl),
`Dados/*.csv`), which are themselves classic exercises from Péllico Netto & Brena's
*Inventário Florestal* (1997). Where the dissertation's spreadsheet reports enough detail
(plot/population size, area), this package's result reproduces its published mean and
total exactly; small differences elsewhere (e.g. `requiredSampleSize`) are expected, since
this package solves that equation with its own iterative loop.

## Simple Random Sampling

[`simpleCasualSampling`](@ref) is the reference design every other method in this section
is compared against: plots are laid out with no further structure at all — no strata, no
clusters, no systematic grid — so each one has an equal, independent chance of being
selected. It is the simplest design to lay out in the field, and the baseline every other,
more elaborate design here is measured against for how much variance reduction its extra
structure actually buys.

Give it a vector of the sampled quantity (`value` — volume, basal area, tree count, or
anything else), the average plot area, and the total area. The optional keyword arguments
`e` (desired relative sampling error, %) and `α` (confidence level) drive the
required-sample-size calculation reported alongside the estimate.

The 70 plot volumes below (1000 m² plots, a 45 ha stand) come from Nardini's own
simple-random-sampling exercise — our `mean`/`total` reproduce the dissertation's
published 217.97 m³/ha and 9808.71 m³ exactly:

```@example inv_simple
using ForestInventory, ForestFoundations

v = [
    18.5, 22.8, 15.8, 9.8, 23.5, 20.7, 28.3, 6.9, 28.7, 26.4, 7.6, 10.3, 21.1, 26.3, 39.9,
    12.7, 11.5, 11.0, 18.9, 31.4, 9.6, 9.8, 7.1, 25.9, 25.7, 24.2, 19.0, 25.9, 37.7, 23.3,
    24.6, 32.3, 23.4, 17.1, 8.9, 32.7, 8.1, 29.2, 24.7, 25.6, 13.2, 25.7, 27.0, 20.4, 13.6,
    32.4, 31.3, 34.6, 9.5, 18.7, 28.9, 21.6, 8.3, 23.2, 26.6, 25.6, 8.0, 22.1, 21.9, 23.3,
    9.9, 30.7, 29.2, 28.0, 13.7, 27.6, 37.9, 27.3, 28.0, 30.7,
]u"m^3"

simpleCasualSampling(v, 1000u"m^2", 45u"ha")
```

The result reports, among other statistics, the sample mean, its standard error and
confidence interval, the population total estimate, and how many additional plots would be
needed to reach the target relative error `e` at confidence level `α`. Tightening `e` to 5%
raises `requiredSampleSize` accordingly — same data, `plot_area`/`total_area` now given as
plain numbers, which is equivalent to the call above:

```@example inv_simple
simpleCasualSampling(v, 0.1, 45; e=5)
```

## Stratified Random Sampling

[`stratifiedSampling`](@ref) divides the population into non-overlapping strata before
sampling within each — usually a variance reduction over simple random sampling of the
same total size. It returns a [`SamplingReport`](@ref) with an ANOVA table (is there
really a difference between strata?), an auxiliary table (per-stratum descriptive
statistics and allocation weights), and the final result table. Each stratum's area is
given as a `Dict` keyed by the same stratum labels used in `stratum` — not a plainly
ordered vector — so the plots don't need to be pre-sorted by stratum; they can be listed
in whatever order they were collected in the field.

The data below comes from Nardini's stratified-sampling exercise (3 strata of 7, 8, and 7
plots, 0.1 ha each, again in a 45 ha stand). The dissertation's spreadsheet doesn't record
each stratum's individual area — only the 45 ha total — so it's split evenly here (15 ha
each); consequently the mean/total below land close to, but not exactly on, the
dissertation's published 214.95 m³/ha / 9672.95 m³, which used the book's actual
per-stratum weights:

```@example inv_stratified
using ForestInventory, ForestFoundations

stratum = vcat(fill(1, 7), fill(2, 8), fill(3, 7))
value = [
    15.8, 7.6, 8.8, 12.5, 11.1, 16.2, 12.2,
    20.4, 30.5, 26.8, 27.2, 28.4, 19.7, 20.4, 23.1,
    21.3, 24.3, 29.2, 21.8, 33.1, 35.8, 26.7,
]u"m^3"
strata_area = Dict(1 => 15.0, 2 => 15.0, 3 => 15.0)   # hectares, one entry per stratum label

report = stratifiedSampling(stratum, value, 0.1, strata_area)
resultTable(report)
```

```@example inv_stratified
auxiliaryTable(report)
```

```@example inv_stratified
anova(report)
```

With a single stratum, [`stratifiedSampling`](@ref) reduces exactly to
[`simpleCasualSampling`](@ref) — useful as a sanity check when strata are added
incrementally to a growing dataset.

## Systematic Sampling

[`systematicSampling`](@ref) estimates the variance of the mean from the method of
successive differences between consecutive plots, since plots laid out at a fixed
interval tend to be more alike than a true random sample.

An optional `line` vector groups plots into several independent transects, so the
difference between the last plot of one line and the first plot of the next — which
aren't actually adjacent on the ground — is excluded from the variance estimate. Nardini's
systematic-sampling exercise lays out 5 lines of 15 plots each (0.1 ha, the same 45 ha
stand as the simple-random-sampling example above) — our `mean`/`total` reproduce the
dissertation's published 221.11 m³/ha and 9949.92 m³ exactly:

```@example inv_systematic
using ForestInventory, ForestFoundations

line1 = [8.0, 8.6, 13.1, 16.6, 21.6, 23.6, 19.7, 30.6, 20.4, 28.0, 30.1, 22.6, 26.7, 31.8, 25.5]
line2 = [9.4, 12.7, 7.6, 19.0, 15.1, 23.7, 18.4, 29.4, 24.6, 38.6, 27.6, 24.7, 24.8, 28.7, 33.1]
line3 = [7.3, 9.8, 9.5, 18.6, 15.8, 23.3, 22.8, 27.8, 19.2, 28.3, 23.1, 26.7, 26.9, 33.9, 35.1]
line4 = [10.1, 11.6, 8.5, 10.8, 18.0, 23.0, 27.1, 28.7, 28.2, 34.9, 28.4, 31.7, 31.8, 28.3, 30.1]
line5 = [8.7, 8.8, 10.4, 22.9, 16.22, 24.0, 19.5, 30.5, 19.9, 22.9, 21.4, 28.7, 21.8, 30.7, 28.3]
v = vcat(line1, line2, line3, line4, line5)u"m^3"
line = repeat(1:5, inner=15)

systematicSampling(v, 0.1, 45; line=line)
```

Omitting `line` treats all 75 plots as one continuous transect instead, so the differences
straddling each pair of lines are no longer excluded from the variance estimate:

```@example inv_systematic
systematicSampling(v, 0.1, 45)
```

## Systematic Sampling with Multiple Random Starts

[`multistartSystematicSampling`](@ref) lays out several independent systematic lines,
each starting at its own random point, turning the design into a genuine probability
sample: statistically identical to [`clusterSampling`](@ref)'s between/within
decomposition, with each line playing the role of one cluster.

Nardini's exercise lays out 9 random starts of 8 plots each (0.25 ha, a 13000 ha stand)
— our `mean`/`total` reproduce the dissertation's published 199.08 m³/ha and
2,588,062 m³ exactly:

```@example inv_multistart
using ForestInventory, ForestFoundations

start = repeat(1:9, inner=8)
value = [
    25.97, 34.62, 25.97, 60.27, 49.75, 81.64, 85.90, 72.43,
    54.02, 56.77, 63.78, 61.94, 61.28, 38.93, 44.68, 48.30,
    66.43, 50.46, 43.81, 39.63, 42.47, 38.84, 37.13, 35.93,
    59.37, 56.18, 39.01, 32.93, 95.67, 32.82, 65.09, 46.29,
    62.69, 68.88, 65.84, 85.93, 48.13, 55.98, 70.91, 37.90,
    47.85, 69.66, 78.11, 36.71, 28.09, 53.57, 40.34, 56.22,
    36.64, 39.23, 42.81, 35.51, 27.54, 54.78, 42.56, 88.45,
    34.72, 41.21, 41.79, 34.97, 41.47, 29.11, 25.60, 21.99,
    43.13, 59.13, 92.98, 49.14, 46.68, 42.60, 26.17, 30.14,
]u"m^3"

report = multistartSystematicSampling(start, value, 0.25, 13000)
resultTable(report)
```

## One-Stage Cluster Sampling

[`clusterSampling`](@ref) samples clusters of `M` neighboring plots ("conglomerates")
instead of individual plots — cheaper to lay out in the field, at the cost of
within-cluster homogeneity typically inflating the variance of the mean.

Nardini's exercise samples 25 clusters ("conglomerados") of 4 plots each (0.25 ha, a
4000 ha stand) — our `mean`/`total` reproduce the dissertation's published 256.99 m³/ha
and 1,027,968 m³ exactly:

```@example inv_cluster
using ForestInventory, ForestFoundations

cluster = repeat(1:25, inner=4)
value = [
    84.7, 78.2, 63.0, 88.5, 70.1, 65.4, 59.6, 58.1, 75.6, 77.7, 66.1, 43.5, 81.3, 38.7, 62.8, 31.1,
    99.5, 47.0, 82.1, 77.4, 97.3, 66.1, 57.4, 113.9, 70.2, 43.5, 48.3, 69.3, 81.7, 68.8, 22.6, 78.2,
    60.9, 71.4, 75.2, 81.1, 67.3, 90.5, 72.6, 67.1, 75.4, 88.5, 104.6, 105.5, 71.0, 57.0, 79.3, 71.8,
    61.8, 100.9, 69.4, 61.0, 116.9, 91.0, 36.9, 57.3, 56.0, 22.9, 94.8, 130.2, 68.2, 32.5, 68.2, 58.6,
    55.4, 31.3, 76.3, 40.0, 55.0, 63.0, 56.8, 36.8, 78.8, 59.6, 68.1, 46.8, 33.6, 44.4, 34.1, 79.6,
    48.8, 26.6, 67.5, 51.0, 82.4, 61.8, 23.4, 60.2, 69.8, 49.5, 53.9, 0.0, 58.4, 61.5, 45.5, 50.2,
    27.4, 63.6, 69.4, 60.8,
]u"m^3"

report = clusterSampling(cluster, value, 0.25, 4000)
resultTable(report)
```

```@example inv_cluster
clusterTable(report)
```

## Horizontal Point Sampling (Bitterlich)

[`horizontalPointSampling`](@ref) estimates the value and basal area per hectare from
angle-count ("Bitterlich") data: at each point, every tree that is "in" for the chosen
basal area factor (`baf`) represents exactly `baf` m²/ha of basal area, regardless of its
own size — no fixed plot radius is needed at all. Point 6 below was visited but had no
"in" trees; it is not a row of `point`/`diameter`/`value`, but is still counted through
`npoints=6` so the zero observation is not silently dropped from the mean. Unlike every
other design in this section, this one has no counterpart in Nardini's dissertation — it
wasn't one of the original 10 designs — so the data below is illustrative rather than
sourced from that exercise set.

```@example inv_hps
using ForestInventory, ForestFoundations

point    = [1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 5, 5, 5]
diameter = [25.0, 30.0, 20.0, 28.0, 22.0, 35.0, 30.0, 25.0, 20.0, 22.0, 30.0, 28.0, 26.0]u"cm"
value    = [0.35, 0.55, 0.22, 0.48, 0.28, 0.85, 0.55, 0.35, 0.22, 0.28, 0.55, 0.48, 0.40]u"m^3"

# baf=2 m²/ha, 6 points visited (5 counted trees + 1 empty), 0.1 ha per point, 10 ha stand
report = horizontalPointSampling(point, diameter, value, 2.0, 6, 0.1, 10)
resultTable(report)
```

```@example inv_hps
pointTable(report)
```

## Two-Stage Sampling

[`twoStageSampling`](@ref) draws `n` primary units (e.g. stands) from a population of `N`,
then sub-samples `m` secondary units (plots) from within each drawn primary out of `M`
possible — unlike cluster sampling, not every secondary unit inside a drawn primary needs
to be measured.

Nardini's exercise sub-samples 36 primary units with 6 secondary units (0.25 ha plots)
each, out of a population of `N=100` possible primaries with `M=40000` possible secondary
units each — our `mean`/`total` reproduce the dissertation's published 223.51 m³/ha and
223,509,259 m³ exactly:

```@example inv_twostage
using ForestInventory, ForestFoundations

primary = repeat(1:36, inner=6)
value = [
    29.07, 72.37, 96.03, 55.90, 106.91, 54.51,
    113.55, 53.14, 63.08, 125.64, 113.10, 27.21,
    37.26, 58.57, 40.27, 58.81, 57.30, 72.97,
    52.66, 58.97, 37.86, 74.61, 69.56, 70.67,
    60.81, 37.76, 63.43, 75.40, 46.24, 44.19,
    72.98, 65.00, 94.62, 55.52, 57.81, 48.97,
    120.71, 73.83, 108.75, 47.76, 101.33, 98.52,
    62.48, 70.30, 41.59, 49.85, 65.64, 65.05,
    73.35, 38.46, 138.26, 106.78, 154.13, 145.69,
    56.35, 19.51, 20.25, 23.54, 38.34, 25.06,
    73.65, 51.87, 62.15, 43.96, 30.82, 41.77,
    65.58, 36.48, 38.76, 26.37, 83.33, 47.14,
    34.31, 29.64, 30.69, 23.87, 23.63, 37.50,
    45.98, 50.52, 44.55, 72.07, 42.39, 86.99,
    55.58, 67.72, 44.15, 107.90, 51.05, 101.61,
    26.79, 27.88, 69.14, 36.53, 41.43, 55.25,
    52.11, 47.06, 80.57, 20.37, 166.67, 118.35,
    44.40, 117.04, 82.26, 130.69, 123.53, 131.59,
    41.48, 90.94, 60.43, 47.60, 58.71, 44.05,
    59.69, 65.41, 72.89, 55.52, 41.42, 76.90,
    50.27, 52.25, 31.06, 79.50, 88.10, 60.19,
    62.24, 27.18, 69.97, 72.08, 60.97, 46.94,
    36.71, 31.28, 34.94, 37.25, 64.99, 40.83,
    63.36, 26.13, 25.89, 33.80, 25.12, 60.89,
    43.08, 41.36, 30.66, 32.78, 28.10, 30.42,
    38.29, 37.51, 29.58, 62.67, 56.40, 73.09,
    32.34, 57.92, 89.36, 16.09, 52.05, 54.40,
    75.36, 59.68, 36.97, 22.02, 40.25, 69.87,
    53.63, 21.01, 46.76, 84.00, 88.49, 45.35,
    48.85, 45.06, 30.19, 30.51, 27.35, 48.72,
    21.05, 79.97, 20.08, 62.15, 16.20, 56.69,
    29.13, 27.72, 35.92, 80.50, 98.08, 25.65,
    24.77, 20.32, 62.24, 27.18, 69.97, 72.08,
    47.21, 43.56, 25.74, 23.53, 28.29, 50.16,
    43.11, 40.80, 53.10, 28.08, 54.64, 19.82,
    40.25, 33.02, 32.24, 23.87, 64.41, 18.85,
]u"m^3"

report = twoStageSampling(primary, value, 0.25, 100, 40000)
resultTable(report)
```

## Successive Occasions

The remaining four designs estimate growth between two inventory occasions, differing in
how much of the plot network is remeasured at the second occasion. All four return a
[`SamplingReport`](@ref) with `occasion1`, `occasion2`, and `change` tables.

### Independent Samples

[`independentOccasionsSampling`](@ref) samples each occasion completely independently —
the simplest design, but the least efficient at detecting growth, since it exploits none
of the natural plot-to-plot correlation between occasions.

Nardini's exercise measures 44 plots at each occasion (0.05 ha plots) — our `meanGrowth`
reproduces the dissertation's published 48.21 m³ exactly (the dissertation's result sheet
doesn't record each occasion's population size, so error/CI figures here can differ):

```@example inv_independent
using ForestInventory, ForestFoundations

v1 = [
    146.7, 182.1, 258.6, 198.8, 196.7, 155.6, 201.0, 183.8, 212.4, 202.5, 217.5, 138.9, 242.7,
    200.0, 208.3, 193.0, 279.3, 181.0, 223.4, 210.4, 206.1, 171.2, 264.4, 278.6, 253.3, 41.0,
    242.9, 241.6, 251.2, 144.5, 225.9, 178.5, 182.0, 176.2, 199.5, 156.8, 136.2, 225.4, 184.3,
    245.8, 164.9, 198.8, 168.6, 228.7,
]u"m^3"
v2 = [
    330.4, 222.3, 300.2, 293.8, 270.2, 230.4, 312.3, 341.6, 326.5, 46.3, 312.8, 273.8, 272.3,
    225.2, 282.9, 233.3, 213.4, 224.2, 229.1, 188.0, 176.7, 221.0, 247.5, 249.5, 230.3, 136.4,
    285.4, 265.9, 248.1, 250.2, 238.0, 265.9, 210.6, 280.0, 198.2, 290.3, 248.2, 235.1, 265.2,
    236.7, 280.4, 285.4, 198.3, 248.0,
]u"m^3"

report = independentOccasionsSampling(v1, v2, 0.05, 30000, 30000)
change(report)
```

### Complete Replacement

[`completeReplacementSampling`](@ref) remeasures the exact same plots at both occasions,
exploiting their positive correlation to tighten the growth estimate.

Nardini's exercise remeasures the same 45 plots at both occasions (0.05 ha plots) — our
`meanGrowth` reproduces the dissertation's published 39.19 m³ exactly:

```@example inv_complete
using ForestInventory, ForestFoundations

v1 = [
    141.89, 147.74, 126.26, 137.08, 39.30, 133.19, 110.21, 126.79, 132.09, 73.79, 76.22, 63.42,
    77.47, 37.24, 98.34, 87.04, 116.52, 83.75, 110.99, 77.87, 139.73, 182.11, 94.68, 181.69,
    161.44, 172.58, 191.29, 137.93, 167.08, 110.10, 123.59, 107.45, 137.35, 137.77, 79.56,
    80.79, 118.52, 133.61, 174.91, 140.90, 176.79, 185.32, 135.84, 109.97, 73.65,
]u"m^3"
v2 = [
    174.99, 187.58, 153.81, 176.28, 61.76, 160.87, 150.23, 164.15, 170.99, 101.36, 114.85, 87.89,
    102.06, 51.37, 152.57, 119.94, 148.93, 129.42, 152.01, 109.04, 165.91, 208.33, 143.74, 197.50,
    191.62, 266.13, 280.29, 178.37, 218.88, 149.68, 161.67, 174.45, 172.66, 179.97, 105.36,
    126.26, 150.55, 191.47, 248.38, 183.30, 226.76, 224.25, 154.10, 142.97, 102.68,
]u"m^3"

report = completeReplacementSampling(v1, v2, 0.05, 30000, 30000)
change(report)
```

### Partial Replacement

[`partialReplacementSampling`](@ref) is the middle ground: a matched subset of plots is
remeasured, some temporary plots are dropped, and new temporary plots are added at the
second occasion. `missing` marks a plot that wasn't measured on a given occasion.

Nardini's exercise measures 17 plots only at the first occasion, remeasures 13 matched
plots at both, and adds 13 new plots only at the second (0.05 ha plots) — our
`meanGrowth` (44.99 m³) lands close to, but not exactly on, the dissertation's published
44.32 m³/ha, since this package solves the regression-weighted estimator with its own
iterative sample-size loop:

```@example inv_partial
using ForestInventory, ForestFoundations

value1 = [
    141.9, 147.7, 126.3, 137.1, 133.2, 110.2, 126.8, 132.1, 73.8, 76.2, 63.4, 77.4, 98.3, 87.4,
    116.5, 83.8, 111.0,
    77.9, 139.7, 182.1, 94.7, 181.7, 161.4, 172.6, 191.3, 137.9, 167.1, 110.1, 123.6, 107.5,
    fill(missing, 13)...,
]u"m^3"
value2 = [
    fill(missing, 17)...,
    109.0, 165.9, 208.3, 143.7, 197.5, 191.6, 266.1, 280.3, 178.4, 218.9, 149.7, 161.7, 174.5,
    172.7, 180.0, 105.4, 126.3, 150.6, 191.5, 248.4, 183.3, 226.8, 224.3, 154.1, 142.9, 102.7,
]u"m^3"

report = partialReplacementSampling(value1, value2, 0.05, 10000)
occasion2(report)
```

```@example inv_partial
change(report)
```

### Double Sampling

[`doubleSampling`](@ref) measures a large first-occasion sample but only remeasures a
smaller "permanent" subset at the second occasion, estimating the rest via a regression of
the permanent subset's second-occasion value on its first-occasion value — solved
through the normal equations `(X'X)β = X'y`, the same transposed-design-matrix approach
used for the ANOVA in [`stratifiedSampling`](@ref).

Nardini's exercise measures 56 plots at the first occasion, only remeasuring the last 26
("permanent" subset) at the second (0.05 ha plots) — our `meanGrowth` reproduces the
dissertation's published 49.89 m³ exactly:

```@example inv_double
using ForestInventory, ForestFoundations

value1 = [
    132.3, 122.9, 133.3, 124.8, 127.5, 116.4, 120.0, 140.6, 116.8, 113.9, 138.7, 144.5, 104.0,
    157.4, 116.2, 100.1, 164.4, 195.7, 129.8, 102.1, 140.7, 187.0, 64.6, 129.6, 147.3, 151.9,
    110.7, 147.1, 101.8, 146.5,
    139.1, 190.3, 100.7, 102.7, 98.4, 68.2, 42.2, 91.6, 91.6, 77.4, 91.4, 91.6, 123.5, 108.7,
    108.7, 93.6, 42.0, 53.8, 138.0, 173.1, 121.1, 81.8, 64.0, 145.9, 168.9, 186.8,
]u"m^3"
value2 = [
    fill(missing, 30)...,
    177.6, 232.9, 165.5, 181.8, 164.3, 103.6, 68.0, 145.5, 137.6, 123.7, 146.6, 141.3, 174.4,
    164.9, 164.9, 134.0, 57.4, 84.4, 202.0, 207.9, 171.8, 154.0, 101.2, 189.2, 226.9, 239.3,
]u"m^3"

report = doubleSampling(value1, value2, 0.05, 10000)
occasion2(report)
```

```@example inv_double
change(report)
```
