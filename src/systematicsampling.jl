# Successive differences between consecutive plots along a systematic line, skipping the
# pair that straddles two different lines when a `line` grouping is given -- generic to
# any number of lines (transects), reducing to a single running difference when there is
# only one.
function _successivedifferences(value::AbstractVector, line::Union{AbstractVector,Nothing})
  n = length(value)
  keep = isnothing(line) ? trues(n - 1) : [line[i] == line[i+1] for i in 1:(n-1)]
  return [value[i] - value[i+1] for i in 1:(n-1) if keep[i]]
end

"""
    systematicsampling(value::AbstractVector, plotarea::Area, totalarea::Area;
                        line::Union{AbstractVector,Nothing}=nothing, e::Real=10, α::Real=0.95)

Estimates the mean, the total and their confidence intervals of a forest by systematic sampling in a single stage, with the variance of the mean estimated from the differences between successive plots.

# Description

In systematic sampling the plots are laid out at a fixed interval along one or more lines (strips or transects) instead of at random, which is faster and cheaper in the field. Consecutive plots tend to be more alike than two random ones, so the variance of simple random sampling does not estimate the precision of a systematic sample. With a single random start there is no replication to estimate the error, and the method of the successive differences (Loetsch & Haller, 1973) is the usual approximation. The notation and the formulas follow Péllico Netto & Brena (1997), chapter 7.

# Arguments

- `value`: the measurement of every plot, **in the order in which the plots were measured along the line(s)**. A `Unitful` quantity keeps its unit; a plain number stays unitless.
- `plotarea`: the area of each plot, an `Area`.
- `totalarea`: the area of the forest, an `Area`.
- `line::Union{AbstractVector,Nothing}=nothing`: the line (strip) of every plot, in the order of `value`. The difference between the last plot of a line and the first of the next is left out, because they are not neighbours on the ground. With `nothing` the plots are a single line.
- `e::Real=10`: the admissible sampling error, in percent of the mean. This design reports whether it was met, not a required number of plots.
- `α::Real=0.95`: the confidence level.

# Returns

A [`SamplingReport`](@ref) with `resulttable`, a one-row `DataFrame` with the columns of [`simplecasualsampling`](@ref) (`xm`, `cv`, `s2xm`, `sxm`, `Ea`, `Er`, `xha`, `X`, `Xl`, `Xu`, `population`, `f`, `n`, `N`, `A`), with these differences:

- `cv` is the coefficient of variation of the plots, ``s/\\bar{x} \\times 100``, from their sample variance, since the method of the differences estimates only the variance of the mean;
- `e`: the admissible error, repeated for the comparison with `Er`;
- `m`: the number of lines.

There are no `nreq` and `nadd`: the precision depends on how the variable is arranged along the line, not only on ``n``.

# Mathematical basis

For ``m`` lines with ``n_j`` plots each and ``n = \\sum n_j`` plots, the mean is ``\\bar{x} = \\sum x_i / n`` (Eq. 7.5), and the variance of the mean comes from the squared differences between successive plots of the same line (Eqs. 7.11 and 7.14), with the finite-population correction ``(1 - f)``, ``f = n/N``:

```math
s_{\\bar{x}}^2 = \\frac{\\sum_{j=1}^{m} \\sum_{i=1}^{n_j - 1} (x_{ij} - x_{(i+1)j})^2}{2 n \\sum_{j=1}^{m} (n_j - 1)} (1 - f)
= \\frac{\\sum (x_i - x_{i+1})^2}{2 n (n - m)} (1 - f)
```

With a single line (``m = 1``) it is the classic ``\\sum (x_i - x_{i+1})^2 / [2n(n-1)]`` of Eq. 7.11. The correction is dropped when ``1 - f \\ge 0.98``. The standard error, the sampling errors, the total and its interval follow the same equations as simple random sampling, with ``t`` of ``n - 1`` degrees of freedom (Eqs. 7.16 to 7.21):

```math
s_{\\bar{x}} = \\sqrt{s_{\\bar{x}}^2} \\qquad
E_a = t \\, s_{\\bar{x}} \\qquad E_r = \\frac{t \\, s_{\\bar{x}}}{\\bar{x}} \\times 100 \\qquad
\\hat{X} = N \\bar{x} \\qquad
IC\\left[\\hat{X} - N t s_{\\bar{x}} \\le X \\le \\hat{X} + N t s_{\\bar{x}}\\right] = P
```

# Technical description

Systematic sampling has no closed formula for the sampling intensity like simple or stratified sampling, since its precision depends on the arrangement of the variable and not only on ``n``. This function reports whether the achieved ``E_r`` meets the admissible error `e`. The method of the differences needs the plots in the order of the line, and it is an approximation: it is exact for a random arrangement and conservative when there is a trend. The lines are the strips of the book: each one is a sequence of adjacent plots.

# Examples

The example of the book (chapter 7, Quadro 7.1): 75 plots of 0.1 ha in a stand of 45 ha (``N = 450``), in 15 lines of 5 plots, in m³ per plot. Each row of the matrix is a line, in the order of the plots.

```julia-repl
julia> rows = [8.0 9.4 7.3 10.1 8.7;
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
        25.5 33.1 35.1 30.1 28.3];
julia> v = vec(permutedims(rows)) * u"m^3";
julia> line = repeat(1:15, inner=5);

julia> report = systematicsampling(v, 0.1u"ha", 45u"ha"; line=line)
julia> report.resulttable.xm[1]
22.110666666666667 m^3

julia> report.resulttable.s2xm[1]
0.13944814814814818 m^6

julia> report.resulttable.Er[1]
3.37
```

The book gets ``\\bar{x} = 22.11`` m³, ``s_{\\bar{x}}^2 = 0.1394`` m⁶, ``E_r = 3.38`` % (with the tabulated ``t = 2.0``, the package uses the exact ``t = 1.993``) and ``\\hat{X} = 9\\,950`` m³, with ``\\sum (x_i - x_{i+1})^2 = 1\\,506.04``.

# References
Péllico Netto, S. & Brena, D. A. (1997). *Inventário Florestal*. Curitiba, chapter 7.
Loetsch, F. & Haller, K. E. (1973). *Forest Inventory*, vol. 1. München: BLV.
"""
function systematicsampling(value::AbstractVector, plotarea::Area, totalarea::Area;
  line::Union{AbstractVector,Nothing}=nothing, e::Real=10, α::Real=0.95)

  n = length(value)
  n > 1 || throw(ArgumentError("At least two sampled plots are required."))
  checklengths(:line => line, :value => value)

  N = round(Int, NoUnits(totalarea / plotarea))
  k = isnothing(line) ? 1 : length(unique(line))
  x̅ = mean(value)
  diffs = _successivedifferences(value, line)

  infinite = _isinfinitepopulation(n, N)
  f = infinite ? 1.0 : 1 - n / N
  population = infinite ? "infinite" : "finite"
  s²x̅ = sum(abs2, diffs) / (2n * (n - k)) * f
  sx̅ = sqrt(s²x̅)
  cv = std(value) / x̅ * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  absoluteerror = Ttab * sx̅
  relativeerror = absoluteerror / x̅ * 100
  hectarevalue = x̅ * expansionfactor(plotarea)
  totalvalue = x̅ * N
  ciupper = totalvalue + N * absoluteerror
  cilower = totalvalue - N * absoluteerror

  resulttable = DataFrame(
    xm=x̅, cv=round(cv, digits=2), s2xm=s²x̅, sxm=sx̅,
    Ea=absoluteerror, Er=round(relativeerror, digits=2), e=e,
    xha=hectarevalue, X=totalvalue, Xl=cilower, Xu=ciupper,
    population=population, f=round(n / N, digits=3), m=k, n=n, N=N, A=totalarea,
  )
  return SamplingReport((; resulttable))
end
