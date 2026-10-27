"""
    SamplingReport

Wraps the set of tables produced by a multi-table forest-inventory sampling design (for
example a per-stratum auxiliary table alongside the final result table) so they print and
travel together instead of as a bare `Vector`/`Tuple`. Designs that only ever produce one
table (e.g. [`simplecasualsampling`](@ref)) return a plain `DataFrame` instead —
wrapping a single table would add nothing.

# Fields
- `tables::NamedTuple`: the report's tables, keyed by name, e.g.
  `(anova=..., auxiliarytable=..., resulttable=...)` for a stratified design, or
  `(occasion1=..., occasion2=..., change=...)` for a successive-occasions design.

# Examples

Individual tables are reachable both through `.tables` and directly as properties:

```julia-repl
julia> report = stratifiedsampling(stratum, value, plotarea, strataarea);
julia> report.resulttable   # shorthand for report.tables.resulttable
julia> report.auxiliarytable
```
"""
struct SamplingReport
  tables::NamedTuple
end

Base.propertynames(r::SamplingReport) = (:tables, propertynames(getfield(r, :tables))...)

function Base.getproperty(r::SamplingReport, name::Symbol)
  name === :tables && return getfield(r, :tables)
  return getfield(r, :tables)[name]
end

function Base.show(io::IO, r::SamplingReport)
  tables = getfield(r, :tables)
  names = propertynames(tables)
  for (i, name) in enumerate(names)
    println(io, titlecase(replace(string(name), r"table$" => " table")))
    show(io, tables[name])
    i < length(names) && (println(io); println(io))
  end
end

"""
    resulttable(report::SamplingReport)

The final one-row, one-column-per-statistic result table — same as `report.resulttable`.
Produced by every multi-table design except the four successive-occasions designs
(see [`occasion1`](@ref)/[`occasion2`](@ref)/[`change`](@ref) instead). The function-call
equivalent of dot access, matching how `AllometricModel`'s own results are reached
through functions (`coef`, `confint`, ...) rather than raw field access elsewhere in
this ecosystem; calling it on a report that has no such table errors the same way
`report.resulttable` already would.

# Examples

```julia-repl
julia> report = clustersampling(cluster, value, plotarea, totalarea);

julia> resulttable(report)   # same as report.resulttable
```
"""
resulttable(r::SamplingReport) = r.resulttable

"""
    clustertable(report::SamplingReport)

Per-cluster descriptive statistics from [`clustersampling`](@ref): `n` (plots), `xmi` (mean) and `s2i` (variance) of each cluster — same as `report.clustertable`.
"""
clustertable(r::SamplingReport) = r.clustertable

"""
    auxiliarytable(report::SamplingReport)

Per-stratum statistics from [`stratifiedsampling`](@ref): `nh` (plots), `xmh` (mean), `s2h` (variance), `sh` (standard deviation) and the weights `Wh`, `Whsh` and `Whs2h` — same as `report.auxiliarytable`.
"""
auxiliarytable(r::SamplingReport) = r.auxiliarytable

"""
    pointtable(report::SamplingReport)

Per-point statistics (the line of [`strandsampling`](@ref), or the point of [`prodansampling`](@ref), in the same shape): `n` (trees counted), `G` (basal area per hectare) and `xha` (value per hectare) from
[`horizontalpointsampling`](@ref) (only points with at least one counted tree) — same
as `report.pointtable`.
"""
pointtable(r::SamplingReport) = r.pointtable

"""
    starttable(report::SamplingReport)

Per-start descriptive statistics from [`multistartsystematicsampling`](@ref): `n`, `xmi` and `s2i` of each start — same as
`report.starttable`.
"""
starttable(r::SamplingReport) = r.starttable

"""
    primarytable(report::SamplingReport)

Per-primary-unit descriptive statistics from [`twostagesampling`](@ref): `n`, `xmi` and `s2i` of each primary unit — same as
`report.primarytable`.
"""
primarytable(r::SamplingReport) = r.primarytable

"""
    anova(report::SamplingReport)

The between/within-strata analysis-of-variance table from [`stratifiedsampling`](@ref), with `source`, `gl`, `SQ`, `QM`, `F` and `p` —
same as `report.anova`.
"""
anova(r::SamplingReport) = r.anova

"""
    occasion1(report::SamplingReport)

The first-occasion simple-random-sampling result table, from any of the four
sampling-on-successive-occasions designs — same as `report.occasion1`.
"""
occasion1(r::SamplingReport) = r.occasion1

"""
    occasion2(report::SamplingReport)

The second-occasion result table, from any of the four sampling-on-successive-occasions
designs — same as `report.occasion2`.
"""
occasion2(r::SamplingReport) = r.occasion2

"""
    change(report::SamplingReport)

The estimated growth (mean and total change) between occasions, from any of the four
sampling-on-successive-occasions designs, with `dm`, `sdm`, `Ea`, `Er`, `D`, `Dl` and `Du` — same as `report.change`.
"""
change(r::SamplingReport) = r.change
