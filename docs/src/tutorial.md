# Getting Started

Estimate a stand's total — with its confidence interval, coefficient of variation, and
required sample size — from any of the 14 sampling designs, each its own function
taking plain vectors. Every design accepts `plotarea`/`totalarea` as a `Unitful` area
(`ha`, `acre`, `m^2`, ...); a plain number is not accepted for them. This never extends to the sampled quantity itself (usually named `value`): it might
be volume, basal area, a biomass, a tree count, or anything else, so a `Unitful` quantity is preserved
as given and a plain number stays unitless all the way through: a count of trees has no physical unit, and the
mean per hectare of plain numbers comes in `ha^-1`, trees per hectare. Every design returns a [`SamplingReport`](@ref): the result table and, when the design
produces them, the related tables (one per stratum/cluster, or one per inventory occasion). The tables are reachable through accessor functions
(`resulttable(report)`, `clustertable(report)`, ...) or directly as properties
(`report.resulttable`).

The examples below walk through every design, from the simplest — simple random sampling —
to the plotless methods (horizontal point, horizontal line, Prodan) and 3P sampling. The examples of the
designs of Péllico Netto & Brena's *Inventário Florestal* (1997), chapters 5 to 11, are the worked
examples of the book, with its data (checked against the sums and statistics that the book prints), and
the package reproduces the means, the variances, the totals and the sampling intensities of the book.
The sampling errors and the limits differ in the second or third digit, since the book reads Student's *t*
from a table, with the degrees of freedom rounded, and the package computes it (see the CHANGELOG for the
few differences that remain). The volumes of the book of the successive occasions are per hectare, so they
are given here as the volumes of units of 1 ha.

## Notation

The columns of every table use the symbols of Péllico Netto & Brena (1997) in ASCII, and each docstring gives the complete formulas of its design. The common ones: `n` and `N` are the sampled and the possible units, `f = n/N` the sampling fraction, `xm` the mean ``\bar{x}`` per unit and `xha` the mean per hectare, `cv` the coefficient of variation (%), `s2xm` and `sxm` the variance and the standard error of the mean, `Ea` and `Er` the absolute and relative sampling errors, `X` the estimated total ``\hat{X} = N\bar{x}`` with the limits `Xl` and `Xu` of its confidence interval, and `nreq` and `nadd` the sampling intensity required for the admissible error `e` and the units still to measure. For several strata, clusters or occasions the index follows the symbol: `nh`, `xmh`, `s2h` and `Wh` are the size, the mean, the variance and the weight of the stratum ``h``, `s2e` and `s2d` the variances between and within clusters, and `dm` and `D` the mean and the total change between occasions.

## Simple Random Sampling

[`simplecasualsampling`](@ref) is the reference design every other method in this section
is compared against: plots are laid out with no further structure at all — no strata, no
clusters, no systematic grid — so each one has an equal, independent chance of being
selected. It is the simplest design to lay out in the field, and the baseline every other,
more elaborate design here is measured against for how much variance reduction its extra
structure actually buys.

Give it a vector of the sampled quantity (`value` — volume, basal area, tree count, or
anything else), the area of a plot, and the total area. The optional keyword arguments
`e` (admissible sampling error, % of the mean) and `α` (confidence level) drive the
required-sample-size calculation reported alongside the estimate.

The example is that of chapter 5 of the book: a population of 450 plots of 0.1 ha (a stand of 45 ha) of
*Pinus* sp., sampled with the 70 plots of Quadro 5.2, in m³ per plot. The book gets a mean of 21.797 m³,
`cv` = 39.58 %, a standard error of 0.947 m³ and a total of 9 808.71 m³:

```@example inv_simple
using ForestInventory, ForestFoundations

v = [
    18.5, 22.8, 15.8, 9.8, 23.5, 20.7, 28.3, 6.9, 28.7, 26.4, 7.6, 10.3, 21.1, 26.3, 39.9, 12.7,
    11.5, 11.0, 18.9, 31.4, 9.6, 9.8, 7.1, 25.9, 25.7, 24.2, 19.0, 25.9, 37.7, 23.3, 24.6,
    32.3, 23.4, 17.1, 8.9, 32.7, 8.1, 29.2, 24.7, 25.6, 13.2, 25.7, 27.0, 20.4, 13.6, 32.4,
    31.3, 34.6, 9.5, 18.7, 28.9, 21.6, 8.3, 23.2, 26.6, 25.6, 8.0, 22.1, 21.9, 23.3, 9.9,
    30.7, 29.2, 28.0, 13.7, 27.6, 37.9, 27.3, 28.0, 30.7
]u"m^3"

simplecasualsampling(v, 0.1u"ha", 45u"ha")
```

The result reports, among other statistics, the sample mean, its standard error and
confidence interval, the population total estimate, and how many additional plots would be
needed to reach the admissible error `e` at confidence level `α`. The book's pilot inventory, the first 20 plots
(Quadro 5.1), gives a mean of 19.605 m³ and `cv` = 45.64 %, and requires the 70 plots of the definitive inventory
(`nreq`):

```@example inv_simple
pilot = v[1:20]
simplecasualsampling(pilot, 0.1u"ha", 45u"ha")
```

The result table has a column per statistic, which suits putting many inventories together
(`vcat`) and exporting with `stripunits`. For the vertical layout of the reports, a row per
statistic with its unit apart from the value, use `longtable`:

```@example inv_simple
longtable(simplecasualsampling(v, 0.1u"ha", 45u"ha"))
```

A variable with no unit, such as the number of trees of each plot, is given as plain numbers; the mean per
hectare (`xha`) comes in trees per hectare:

```@example inv_simple
trees = [23, 31, 27, 19, 35, 28, 30, 25, 22, 29]
simplecasualsampling(trees, 0.05u"ha", 10u"ha")
```

## Stratified Random Sampling

[`stratifiedsampling`](@ref) divides the population into non-overlapping strata before
sampling within each — usually a variance reduction over simple random sampling of the
same total size. It returns a [`SamplingReport`](@ref) with an ANOVA table (is there
really a difference between strata?), an auxiliary table (per-stratum descriptive
statistics and allocation weights), and the final result table. Each stratum's area is
given as a `Dict` keyed by the same stratum labels used in `stratum` — not a plainly
ordered vector — so the plots don't need to be pre-sorted by stratum; they can be listed
in whatever order they were collected in the field.

The example is that of chapter 6 of the book: the same population of 450 plots of 0.1 ha, in 3 strata of 144,
164 and 142 plots (14.4, 16.4 and 14.2 ha), sampled with the 22 plots of Quadro 6.2. The book gets a stratified
mean of 21.465 m³, a standard error of 0.916 m³, a total of 9 659 m³ and ``F = 24.42``:

```@example inv_stratified
using ForestInventory, ForestFoundations

stratum = vcat(fill(1, 7), fill(2, 8), fill(3, 7))
value = [
    15.8, 7.6, 8.8, 12.5, 11.1, 16.2, 12.2, 20.4, 30.5, 26.8, 27.2, 28.4, 19.7, 20.4, 23.1, 21.3,
    24.3, 29.2, 21.8, 33.1, 35.8, 26.7
]u"m^3"
strataarea = Dict(1 => 14.4u"ha", 2 => 16.4u"ha", 3 => 14.2u"ha")   # one entry per stratum label

report = stratifiedsampling(stratum, value, 0.1u"ha", strataarea)
resulttable(report)
```

```@example inv_stratified
auxiliarytable(report)
```

```@example inv_stratified
anova(report)
```

The plots are distributed among the strata proportionally to their weights `Wh`, as in the book's first
stratification; `allocation=:neyman` uses Neyman's optimum allocation instead:

```@example inv_stratified
resulttable(stratifiedsampling(stratum, value, 0.1u"ha", strataarea; allocation=:neyman))
```

With a single stratum, [`stratifiedsampling`](@ref) reduces exactly to
[`simplecasualsampling`](@ref) — useful as a sanity check when strata are added
incrementally to a growing dataset.

## Systematic Sampling

[`systematicsampling`](@ref) estimates the variance of the mean from the method of
successive differences between consecutive plots, since plots laid out at a fixed
interval tend to be more alike than a true random sample.

An optional `line` vector groups plots into several independent lines (strips), so the
difference between the last plot of one line and the first plot of the next — which
aren't actually adjacent on the ground — is excluded from the variance estimate. The example is that of
chapter 7 of the book (Quadro 7.1): 75 plots of 0.1 ha in the same stand of 45 ha, in 15 lines of 5 plots
(each row of the matrix is a line, in the order of the plots). The book gets a mean of 22.11 m³, a variance of
the mean of 0.1394 m⁶ and a total of 9 950 m³:

```@example inv_systematic
using ForestInventory, ForestFoundations

rows = [
    8.0 9.4 7.3 10.1 8.7;
    8.6 12.7 9.8 11.6 8.8;
    13.1 7.6 9.5 8.5 10.4;
    16.6 19.0 18.6 10.8 22.9;
    21.6 15.1 15.8 18.0 16.2;
    23.6 23.7 23.3 23.0 24.0;
    19.7 18.4 22.8 27.1 19.5;
    30.6 29.4 27.8 28.7 30.5;
    20.4 24.6 19.2 28.2 19.9;
    28.0 38.6 28.3 34.9 22.9;
    30.1 27.6 23.1 28.4 21.4;
    22.6 24.7 26.7 31.7 28.7;
    26.7 24.8 26.9 31.8 21.8;
    31.8 28.7 33.9 28.3 30.7;
    25.5 33.1 35.1 30.1 28.3
]
value = vec(permutedims(rows)) * u"m^3"    # the plots of each line, in order
line = repeat(1:15, inner=5)

systematicsampling(value, 0.1u"ha", 45u"ha"; line=line)
```

The admissible error `e` is only used to report whether it was met (`Er <= e`), since the precision of a
systematic sample depends on how the variable is arranged along the line.

## Systematic Sampling with Multiple Random Starts

[`multistartsystematicsampling`](@ref) lays out several independent systematic lines, each from its own random
start, so every start is a cluster of plots and the calculations are those of cluster sampling. The example is
that of chapter 10 of the book (Quadro 10.1): 9 starts of 8 plots of 0.25 ha in a stand of 13 000 ha, so 6 500
possible starts. The book gets a mean of 49.769 m³, ``s_e^2 = 27.487`` m⁶, ``r = 0.088``, a variance of the mean
of 6.987 m⁶, a total of 2 587 988 m³ and 13 starts:

```@example inv_multistart
using ForestInventory, ForestFoundations

rows = [
    25.97 34.62 25.97 60.27 49.75 81.64 85.90 72.43;
    54.02 56.77 63.78 61.94 61.28 38.93 44.68 48.30;
    66.43 50.46 43.81 39.63 42.47 38.84 37.13 35.93;
    59.33 56.18 39.01 32.93 95.67 32.82 65.09 46.29;
    62.69 68.88 65.84 85.93 48.13 55.98 70.91 37.90;
    47.85 69.66 78.11 36.71 28.09 53.57 40.34 56.22;
    36.64 39.23 42.81 35.51 27.54 54.78 42.56 88.45;
    34.72 41.21 41.79 34.97 41.47 29.11 25.60 21.99;
    43.13 59.13 92.98 49.14 46.68 42.60 26.17 30.14
]
value = vec(permutedims(rows)) * u"m^3"    # the plots of each start, in order
start = repeat(1:9, inner=8)

report = multistartsystematicsampling(start, value, 0.25u"ha", 13000u"ha")
resulttable(report)
```

## One-Stage Cluster Sampling

[`clustersampling`](@ref) treats a group of `M` neighboring plots as the sampling unit. The example is that of
chapter 9 of the book (Quadro 9.1): 25 clusters of 4 plots of 0.25 ha in a stand of 4 000 ha, so 4 000 possible
clusters. The book gets a mean of 64.548 m³, ``s_d^2 = 432.58`` m⁶, ``s_e^2 = 46.25`` m⁶, ``r = 0.0966``, a
variance of the mean of 6.176 m⁶, a total of 1 032 768 m³ and 17 clusters:

```@example inv_cluster
using ForestInventory, ForestFoundations

rows = [
    84.7 78.2 63.0 88.5;
    70.1 65.4 59.6 58.1;
    75.6 77.7 66.1 43.5;
    81.3 38.7 62.8 31.1;
    99.5 47.0 82.1 77.4;
    97.3 66.1 57.4 113.9;
    70.2 43.5 48.3 69.3;
    81.7 68.8 22.6 78.2;
    60.9 71.4 75.2 81.1;
    67.3 90.5 72.6 67.1;
    75.4 88.5 104.6 105.5;
    71.0 57.0 79.3 71.8;
    61.8 100.9 69.4 61.0;
    116.9 91.0 36.9 57.3;
    56.0 22.9 94.8 130.2;
    68.2 32.5 68.2 58.6;
    55.4 61.3 76.3 40.0;
    55.0 63.0 56.8 36.8;
    78.8 59.6 68.1 46.8;
    33.6 44.4 34.1 79.6;
    48.8 26.6 67.5 51.0;
    82.4 61.8 23.4 60.2;
    69.8 49.5 53.9 0.0;
    58.4 61.5 45.5 50.2;
    27.4 63.6 69.4 60.8
]
value = vec(permutedims(rows)) * u"m^3"    # the plots of each cluster, in order
cluster = repeat(1:25, inner=4)

report = clustersampling(cluster, value, 0.25u"ha", 4000u"ha")
resulttable(report)
```

```@example inv_cluster
clustertable(report)
```

## Horizontal Point Sampling (Bitterlich)

[`horizontalpointsampling`](@ref) estimates the value and basal area per hectare from
angle-count ("Bitterlich") data: at each point, every tree that is "in" for the chosen
basal area factor (`baf`) represents exactly `baf` m²/ha of basal area, regardless of its
own size — no fixed plot radius is needed at all. Point 6 below was visited but had no
"in" trees; it is not a row of `point`/`diameter`/`value`, but is still counted through
`npoints=6` so the zero observation is not silently dropped from the mean. This design is not in the book, so the data below is illustrative.

```@example inv_hps
using ForestInventory, ForestFoundations

point    = [1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 5, 5, 5]
diameter = [25.0, 30.0, 20.0, 28.0, 22.0, 35.0, 30.0, 25.0, 20.0, 22.0, 30.0, 28.0, 26.0]u"cm"
value    = [0.35, 0.55, 0.22, 0.48, 0.28, 0.85, 0.55, 0.35, 0.22, 0.28, 0.55, 0.48, 0.40]u"m^3"

# baf=2 m²/ha, 6 points visited (5 counted trees + 1 empty), 0.1 ha per point, 10 ha stand
report = horizontalpointsampling(point, diameter, value, 2.0u"m^2/ha", 6, 0.1u"ha", 10u"ha")
resulttable(report)
```

```@example inv_hps
pointtable(report)
```


## Two-Stage Sampling

[`twostagesampling`](@ref) draws `n` primary units (e.g. stands) from a population of `N`, and measures `m` of the `M`
plots of each. The example is that of chapter 8 of the book: a forest of 1 000 000 ha, with 100 primary units of
10 000 ha, each of 40 000 plots of 0.25 ha. The pilot inventory (Quadro 8.2, the first 18 primary units of 6 plots)
requires 36 primary units, which the definitive inventory (Quadro 8.3) measures. The book gets a mean of 55.877 m³,
``s_d^2 = 535.637`` m⁶ and a total of 223 508 000 m³ (its ``s_e^2`` and variance of the mean are slightly off, see the
CHANGELOG):

```@example inv_twostage
using ForestInventory, ForestFoundations

rows = [
    29.07 72.37 96.03 55.90 106.91 54.51;
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
    40.25 33.02 32.24 23.87 64.41 18.85
]
value = vec(permutedims(rows)) * u"m^3"    # the plots of each primary unit, in order
primary = repeat(1:36, inner=6)

# the pilot: the first 18 primary units
pilot = resulttable(twostagesampling(primary[1:108], value[1:108], 0.25u"ha", 100, 40000))
pilot.nreq
```

```@example inv_twostage
report = twostagesampling(primary, value, 0.25u"ha", 100, 40000)
resulttable(report)
```

## Horizontal Line Sampling (Strand)

[`strandsampling`](@ref) walks a line of fixed length and counts, with the same angle gauge,
every tree whose diameter subtends the critical angle on either side of the line. Unlike
point sampling, the factor of a counted tree depends on its own diameter (trees per hectare
``TF_i = 500\,000\,k/(D_i L)``): it is the sum of the diameters, not of the basal areas,
that the method estimates directly. Line 4 below was walked but had no counted trees, and is
counted through `nlines = 4`:

```@example inv_strand
using ForestInventory, ForestFoundations

line     = [1, 1, 1, 2, 2, 3]
diameter = [28.1, 33.2, 31.0, 40.1, 27.9, 34.6]u"cm"
value    = [0.65, 1.01, 0.88, 1.78, 0.60, 1.19]u"m^3"

report = strandsampling(line, diameter, value, 2.0u"m^2/ha", 20.0u"m", 4, 10u"ha")
resulttable(report)
```

```@example inv_strand
pointtable(report)
```


## Prodan's k-Tree Sampling

[`prodansampling`](@ref) measures, at each point, the `k` nearest trees and the distance to
the farthest. That distance defines a circular plot whose boundary tree counts half, so the
tree density is ``N_j = (k - 0.5)/(\pi L_j^2)``. Since the estimator is exactly unbiased only
for a random pattern of the trees, the result is the mean and the standard deviation among
points (`sN`, `sG`, `sxha`), with no confidence interval or sampling intensity:

```@example inv_prodan
using ForestInventory, ForestFoundations

point    = [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
diameter = [30.56, 30.40, 36.29, 58.25, 29.60, 36.29, 28.10, 33.20, 31.05, 40.12, 27.88, 34.60]u"cm"
value    = [0.8108, 0.7982, 1.3397, 4.8380, 0.7367, 1.3397, 0.6534, 1.0102, 0.8811, 1.7765, 0.6011, 1.1893]u"m^3"
distance = [2.1, 2.8, 3.4, 3.9, 4.3, 5.2, 1.9, 2.6, 3.1, 3.7, 4.1, 4.9]u"m"

report = prodansampling(point, diameter, value, distance, 10u"ha")
resulttable(report)
```

```@example inv_prodan
pointtable(report)
```


## 3P Sampling

[`threepsampling`](@ref) takes a 100% cruise in which every tree has a quick visual
prediction, and only the trees selected with probability proportional to the prediction have
an exact measurement (`missing` for the others). The total is the sum of all the predictions
times the mean ratio between measurement and prediction, ``\hat{X} = K\bar{R}``:

```@example inv_threep
using ForestInventory, ForestFoundations

predicted = [0.20, 0.45, 0.30, 0.55, 0.15, 0.60, 0.25, 0.40, 0.35, 0.50]u"m^3"
measured  = [missing, 0.48u"m^3", missing, 0.58u"m^3", missing,
             0.63u"m^3", missing, missing, missing, 0.52u"m^3"]

threepsampling(predicted, measured, 5u"ha")
```


## Successive Occasions

The four designs for continuous forest inventory estimate the mean and the total of each occasion and the change
(growth) between them: `independentoccasionssampling`, `completereplacementsampling`, `doublesampling` and
`partialreplacementsampling`. Each returns a [`SamplingReport`](@ref) with the tables `occasion1`, `occasion2` and
`change`. The examples are those of chapter 11 of the book; its volumes are per hectare, so they are given as the
volumes of units of 1 ha.

### Independent Samples

[`independentoccasionssampling`](@ref) measures two independent samples, so the variance of the change is the sum of
the variances of the two means. The example is Quadro 11.1: *Pinus elliottii* in Guarapuava, PR, with 44 plots in 1975 and
44 in 1976 of a population of 1 500 units. The book gets the means 199.98 and 248.19 m³, a change of 48.21 m³ with a
standard error of 10.22 m³ and a total change of 72 315 m³:

```@example inv_independent
using ForestInventory, ForestFoundations

v1 = [
    146.7, 182.1, 258.6, 198.8, 196.7, 155.6, 201.0, 183.8, 212.4, 202.5, 217.5, 138.9, 242.7,
    200.0, 208.3, 193.0, 279.3, 181.0, 223.4, 210.4, 206.1, 171.2, 264.4, 278.6, 253.3, 41.0,
    242.9, 241.6, 251.2, 144.5, 225.9, 178.5, 182.0, 176.2, 199.5, 156.8, 136.2, 225.4,
    184.3, 245.8, 164.9, 198.8, 168.6, 228.7
]u"m^3"
v2 = [
    330.4, 222.3, 300.2, 293.8, 270.2, 230.4, 312.3, 341.6, 326.5, 46.3, 312.8, 273.8, 272.3,
    225.2, 282.9, 233.3, 213.4, 224.2, 229.1, 188.0, 176.7, 221.0, 247.5, 249.5, 230.3,
    136.4, 285.4, 265.9, 248.1, 250.2, 238.0, 265.9, 210.6, 280.0, 198.2, 290.3, 248.2,
    235.1, 265.2, 236.7, 280.4, 285.4, 198.3, 248.0
]u"m^3"

report = independentoccasionssampling(v1, v2, 1u"ha", 1500, 1500)
occasion1(report)
```

```@example inv_independent
occasion2(report)
```

```@example inv_independent
change(report)
```

### Complete Replacement

[`completereplacementsampling`](@ref) remeasures the same plots, so the variance of the change subtracts the covariance
between the occasions. The example is Quadro 11.2: *Eucalyptus* sp. in Lençóis Paulista, SP, with 45 permanent plots measured in
1976 and 1977 of a population of 1 500 units. The book gets the means 121.15 and 160.34 m³, a change of 39.19 m³ with a
standard error of 1.884 m³ and a total change of 58 784 m³:

```@example inv_complete
using ForestInventory, ForestFoundations

v1 = [
    141.89, 147.74, 126.26, 137.08, 39.30, 133.19, 110.21, 126.79, 132.09, 73.79, 76.22, 63.42,
    77.47, 37.24, 98.34, 87.04, 116.52, 83.75, 110.99, 77.87, 139.73, 182.11, 94.68, 181.69,
    161.44, 172.58, 191.29, 137.93, 167.08, 110.10, 123.59, 107.45, 137.35, 137.77, 79.56,
    80.79, 118.52, 133.61, 174.91, 140.90, 176.79, 185.32, 135.84, 109.97, 73.65
]u"m^3"
v2 = [
    174.99, 187.58, 153.81, 176.28, 61.76, 160.87, 150.23, 164.15, 170.99, 101.36, 114.85, 87.89,
    102.06, 51.37, 152.57, 119.94, 148.93, 129.42, 152.01, 109.04, 165.91, 208.33, 143.74,
    197.50, 191.62, 266.13, 280.29, 178.37, 218.88, 149.68, 161.67, 174.45, 172.66, 179.97,
    105.36, 126.26, 150.55, 191.47, 248.38, 183.30, 226.76, 224.25, 154.10, 142.97, 102.68
]u"m^3"

report = completereplacementsampling(v1, v2, 1u"ha", 1500, 1500)
change(report)
```

### Partial Replacement

[`partialreplacementsampling`](@ref) remeasures a part of the plots and adds new ones, with the estimators of Ware & Cunia
(1962). `missing` marks the occasion in which a plot was not measured. The example is Quadro 11.6: *Eucalyptus* sp. in
Lençóis Paulista, SP, with 17 plots measured only in 1976, 13 permanent plots and 13 plots measured only in 1977, of a
population of 1 582 plots. The book gets a current mean of 167.92 m³ with a standard error of 7.568 m³, and a change of 44.32 m³
with a standard error of 4.689 m³:

```@example inv_partial
using ForestInventory, ForestFoundations

xu = [141.9, 147.7, 126.3, 137.1, 133.2, 110.2, 126.8, 132.1, 73.8, 76.2, 63.4, 77.4, 98.3, 87.4,
    116.5, 83.8, 111.0]
xm = [77.9, 139.7, 182.1, 94.7, 181.7, 161.4, 172.6, 191.3, 137.9, 167.1, 110.1, 123.6, 107.5]
ym = [109.0, 165.9, 208.3, 143.7, 197.5, 191.6, 266.1, 280.3, 178.4, 218.9, 149.7, 161.7, 174.5]
yn = [172.7, 180.0, 105.4, 126.3, 150.6, 191.5, 248.4, 183.3, 226.8, 224.3, 154.1, 142.9, 102.7]

value1 = vcat(xu, xm, fill(missing, 13)) * u"m^3"
value2 = vcat(fill(missing, 17), ym, yn) * u"m^3"

report = partialreplacementsampling(value1, value2, 1u"ha", 1582)
occasion2(report)
```

```@example inv_partial
change(report)
```

### Double Sampling

[`doublesampling`](@ref) estimates the second occasion by regression, from a large sample of the first occasion and a subsample
that is remeasured. The example is Quadro 11.3: a stand of 500 ha of *Pinus elliottii* (500 plots of 1 000 m²), with 56 plots in 1975, 30
temporary and 26 permanent, and the 26 permanent plots remeasured in 1976. The book gets ``b = 1.097``, a mean of 169.952 m³
with a standard error of 6.664 m³, and a change of 49.86 m³ with a standard error of 2.790 m³:

```@example inv_double
using ForestInventory, ForestFoundations

xu = [132.3, 122.9, 133.3, 124.8, 127.5, 116.4, 120.0, 140.6, 116.8, 113.9, 138.7, 144.5, 104.0,
    157.4, 116.2, 100.1, 164.4, 195.7, 129.8, 102.1, 140.7, 187.0, 64.6, 129.6, 147.3, 151.9,
    110.7, 147.1, 101.8, 146.5]
xm = [139.1, 190.3, 100.7, 102.7, 98.4, 68.2, 42.2, 91.6, 91.6, 77.4, 91.4, 91.6, 123.5, 108.7,
    108.7, 93.6, 42.0, 53.8, 138.0, 173.1, 121.1, 81.8, 64.0, 145.9, 168.9, 186.8]
ym = [177.6, 232.9, 165.5, 180.8, 164.3, 103.6, 68.0, 145.5, 137.6, 123.7, 146.6, 141.3, 174.4,
    164.9, 164.9, 134.0, 57.4, 84.4, 202.0, 207.9, 171.8, 154.0, 101.2, 189.2, 226.9, 239.3]

value1 = vcat(xu, xm) * u"m^3"
value2 = vcat(fill(missing, 30), ym) * u"m^3"

report = doublesampling(value1, value2, 1u"ha", 500)
occasion2(report)
```

```@example inv_double
change(report)
```
