# Per-point Prodan estimate: sorts each point's trees by ascending distance (the
# convention the method's own field protocol already produces -- trees are found and
# measured outward from the point center, so the farthest one is discovered last), then
# applies Prodan's boundary weighting (the k-th tree contributes half) to basal area and
# `value`, and (k - 0.5) trees to density -- all divided by the circular plot area the
# k-th tree's distance defines.
function _prodanpointtable(point::AbstractVector, g::AbstractVector{<:Area}, value::AbstractVector, distance::AbstractVector{<:Length})
  working = DataFrame(point=point, g=g, value=value, distance=distance)
  grouped = groupby(working, :point, sort=true)
  counts = combine(grouped, nrow => :k).k
  k = counts[1]
  k >= 2 || throw(ArgumentError("at least 2 trees per point are required."))
  all(==(k), counts) || throw(ArgumentError("every point must contribute the same number of trees k (got counts $counts)."))

  return combine(grouped) do df
    ordered = sort(df, :distance)
    L = ordered.distance[end]
    # the expansion factor of the plot's circular area (in ha^-1 or ac^-1, picked from L's
    # own unit system) turns every quantity below into a per-hectare (or per-acre) rate
    ef = expansionfactor(π * abs2(L))
    weight = fill(1.0, k)
    weight[end] = 0.5
    (
      k=k, L=L,
      N=(k - 0.5) * ef,
      G=sum(weight .* ordered.g) * ef,
      xha=sum(weight .* ordered.value) * ef,
    )
  end
end

"""
    prodansampling(point::AbstractVector, diameter::AbstractVector{<:Length}, value::AbstractVector,
                    distance::AbstractVector{<:Length}, totalarea::Area)

Estimates the density, the basal area and the mean per hectare of a forest by Prodan's k-tree distance sampling, a nearest-neighbor method without plot.

# Description

At each sample point the observer measures outward from the center the ``k`` nearest trees (Prodan, 1968, used ``k = 6``, hence "Prodan's 6-tree method"; here any ``k \\ge 2``, inferred from the trees supplied per `point`) and records the horizontal distance ``L`` to the farthest, the ``k``-th. That distance defines a circle of radius ``L`` around the point, except that the ``k``-th tree is on its boundary, so the circle holds ``k - 0.5`` trees and not ``k`` (the correction of Prodan, a different convention from the factor ``(k - 1)/k`` of Moore, 1954; see Kershaw et al., 2017, Sec. 9.5.1). Basal area and any measured quantity take the same weight on the boundary tree before they are divided by the area of the circle.

The true inclusion zones of a nearest-neighbor method are irregular Voronoi polygons, which a single distance does not recover, so every density estimator of the method, Prodan's included, is exactly unbiased only for a random (Poisson) pattern of the trees (Payandeh & Ek, 1986). For this reason this function reports the mean and the standard deviation among points but, unlike the other designs of the package, **no confidence interval and no required sample size**: those calculations assume an unbiased estimator.

# Arguments

- `point::AbstractVector`: the sample point of every measured tree, in the order of `diameter`/`value`/`distance`. Every point has the same number of trees ``k`` (at least 2).
- `diameter::AbstractVector`: the diameter of every tree, `Length` quantities (any unit, imperial included).
- `value::AbstractVector`: the measurement of every tree (typically the volume). A `Unitful` quantity keeps its unit; a plain number stays unitless.
- `distance::AbstractVector`: the horizontal distance from the point to every tree, `Length` quantities. The trees need not be sorted: each point is sorted by distance and the farthest fixes the radius.
- `totalarea`: the area of the forest, an `Area`.

# Returns

A [`SamplingReport`](@ref) with `pointtable`, one row per point with `k` (trees measured), `L` (distance of the ``k``-th tree), `N` (trees per hectare), `G` (basal area per hectare) and `xha` (value per hectare), and `resulttable`, one row with:

- `N`, `G`, `xha`: the mean among points of each per-point estimate.
- `sN`, `sG`, `sxha`: the standard deviation among points of each: only a descriptive spread, not a confidence interval.
- `Gt`, `X`: ``G`` and `xha` times the area of the forest: the total basal area and the total of `value`.
- `k`: the trees measured per point. `n`: the number of points. `A`: the area of the forest.

# Mathematical basis

For point ``j``, with the basal areas ``g_{j,1}, \\dots, g_{j,k}`` and the values ``x_{j,1}, \\dots, x_{j,k}`` ordered from the nearest to the farthest tree and the distance of the boundary ``L_j``:

```math
N_j = \\frac{k - 0.5}{\\pi L_j^2} \\qquad
G_j = \\frac{\\sum_{i=1}^{k-1} g_{j,i} + g_{j,k}/2}{\\pi L_j^2} \\qquad
x_{ha,j} = \\frac{\\sum_{i=1}^{k-1} x_{j,i} + x_{j,k}/2}{\\pi L_j^2}
```

The estimates are the arithmetic means among the ``n`` points, and the spreads their standard deviations (``n - 1`` in the denominator):

```math
N = \\frac{\\sum_{j=1}^{n} N_j}{n} \\qquad
s_N = \\sqrt{\\frac{\\sum_{j=1}^{n} (N_j - N)^2}{n - 1}} \\qquad
G_t = G A \\qquad X = x_{ha} A
```

and the same for ``G`` and ``x_{ha}``. The ``k``-th tree contributes half, as it would on average at the edge of a fixed-area plot.

# Technical description

Due to Prodan (1968); see Kershaw et al. (2017, ch. 9, Sec. 9.5.1) for the density estimator and the discussion of its bias, and Péllico Netto & Brena (1997) for the presentation of the sampling statistics that the rest of this package follows. ``k = 6`` is the recommendation of Prodan; Lessard et al. (1994) and Lynch & Rusydi (1999) report acceptable results for ``k \\ge 5``.

# Examples

```julia-repl
julia> point    = [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2];

julia> diameter = [30.56, 30.40, 36.29, 58.25, 29.60, 36.29,
                    28.10, 33.20, 31.05, 40.12, 27.88, 34.60]u"cm";

julia> value    = [0.8108, 0.7982, 1.3397, 4.8380, 0.7367, 1.3397,
                    0.6534, 1.0102, 0.8811, 1.7765, 0.6011, 1.1893]u"m^3";

julia> distance = [2.1, 2.8, 3.4, 3.9, 4.3, 5.2, 1.9, 2.6, 3.1, 3.7, 4.1, 4.9]u"m";

julia> report = prodansampling(point, diameter, value, distance, 10u"ha");

julia> resulttable(report).N[1]
688.3030915777651 ha^-1

julia> resulttable(report).G[1]
67.86870542693703 m^2 ha^-1

julia> pointtable(report)
2×6 DataFrame
 Row │ point  k      L          N              G                  xha
     │ Int64  Int64  Quantity…  Quantity…      Quantity…          Quantity…
─────┼──────────────────────────────────────────────────────────────────────────────
   1 │     1      6      5.2 m   647.45 ha^-1  74.9144 m^2 ha^-1  1082.21 m^3 ha^-1
   2 │     2      6      4.9 m  729.156 ha^-1   60.823 m^2 ha^-1  731.403 m^3 ha^-1
```

The trees need not be sorted by distance: the six rows of point 1 in any order, with `diameter`/`value` permuted to match, give the same result.

# References
Prodan, M. (1968). Punkt-Stichproben für die Forsteinrichtung. *Der Forst- und Holzwirt*, 23(11), 225-226.
Kershaw, J. A., Ducey, M. J., Beers, T. W. & Husch, B. (2017). *Forest Mensuration*, 5th ed. Wiley.
"""
function prodansampling(point::AbstractVector, diameter::AbstractVector{<:Length}, value::AbstractVector, distance::AbstractVector{<:Length}, totalarea::Area)
  checklengths(:point => point, :diameter => diameter, :value => value, :distance => distance)

  d = diameter
  L = distance
  g = sectionalarea.(d)

  pointtable = _prodanpointtable(point, g, value, L)
  n = nrow(pointtable)
  n >= 1 || throw(ArgumentError("at least one sample point is required."))
  k = pointtable.k[1]

  n̄ = mean(pointtable.N)
  ḡ = mean(pointtable.G)
  v̄ = mean(pointtable.xha)
  nsd = n > 1 ? std(pointtable.N) : zero(n̄)
  gsd = n > 1 ? std(pointtable.G) : zero(ḡ)
  vsd = n > 1 ? std(pointtable.xha) : zero(v̄)

  resulttable = DataFrame(
    N=n̄, sN=nsd,
    G=ḡ, sG=gsd, Gt=ḡ * totalarea,
    xha=v̄, sxha=vsd, X=v̄ * totalarea,
    k=k, n=n, A=totalarea,
  )

  return SamplingReport((; pointtable=pointtable, resulttable=resulttable))
end
