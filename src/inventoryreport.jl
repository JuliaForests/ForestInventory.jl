"""
    SamplingReport

The result of every forest-inventory sampling design: its tables (the final result table
and, when the design produces them, the auxiliary tables, for example one per stratum) print
and travel together, and every design is read the same way, `report.resulttable`.

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
    longtable(report::SamplingReport, table::Symbol=:resulttable)

A table of the report in the vertical layout of the reports, a row per variable with its
name, its unit and its value ([`longtable`](@ref) of ForestFoundations). By default the
result table; `table` names another one, e.g. `:auxiliarytable`, or `:occasion1`,
`:occasion2` and `:change` in the successive-occasions designs.

# Examples

```julia-repl
julia> report = simplecasualsampling(v, 0.1u"ha", 45u"ha");

julia> longtable(report)              # the result table, a row per variable

julia> longtable(report, :change)     # another table of the report
```
"""
function ForestFoundations.longtable(r::SamplingReport, table::Symbol=:resulttable)
  tables = getfield(r, :tables)
  haskey(tables, table) ||
    throw(ArgumentError("The report has no table $table; its tables are $(join(keys(tables), ", "))."))
  return longtable(tables[table])
end

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
