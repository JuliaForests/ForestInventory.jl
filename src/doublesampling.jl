"""
    doublesampling(value1::AbstractVector, value2::AbstractVector, plotarea::Area, npop::Integer; α::Real=0.95)

Estimates the growth (change) of a forest between two inventory occasions with a regression estimator, from a large sample on occasion 1 and a subsample remeasured on occasion 2, generic to any sample size.

# Description

``n_1 = u + m`` plots are measured on occasion 1 and only ``m`` of them (the permanent plots, the subsample) are remeasured on occasion 2 (``n_2 = m < n_1``); the other ``u`` are temporary plots of occasion 1 alone, which exist to fix the mean ``\\bar{x}`` of occasion 1 with more precision at a low cost. The mean of occasion 2 is estimated by the regression of the value of occasion 2 on that of occasion 1 over the permanent plots, applied to the mean of the large sample (Cochran, 1977, Sec. 12.9). The notation and the formulas follow Péllico Netto & Brena (1997), chapter 11, Sec. 11.3.

# Arguments

- `value1`: the measurement of occasion 1 of every plot (``n_1`` entries). A `Unitful` quantity keeps its unit; a plain number stays unitless.
- `value2`: the measurement of occasion 2, one entry per plot in the order of `value1`, and `missing` for the temporary plots that were not remeasured.
- `plotarea`: the area of each plot, an `Area`, used only for the mean per hectare.
- `npop::Integer`: the possible plots of the population.
- `α::Real=0.95`: the confidence level.

# Returns

A [`SamplingReport`](@ref) with `occasion1` (the simple random sample of the ``n_1`` plots, with the columns of [`independentoccasionssampling`](@ref)), `occasion2` (the regression estimate) and `change`. The columns of `occasion2`:

- `xm`: the mean per plot of occasion 2 estimated by regression, ``\\bar{y}_r``.
- `b`: the slope of the regression. `s2yx`: the residual variance of the regression, ``s_{yx}^2``.
- `sxm`, `Ea`, `Er`, `xha`, `X`, `Xl`, `Xu`: as in [`simplecasualsampling`](@ref), for ``\\bar{y}_r``.
- `m`: the permanent (remeasured) plots. `u`: the temporary plots, ``n_1 - m``. `N`: the possible plots.

`change` has the columns of [`independentoccasionssampling`](@ref).

# Mathematical basis

With ``\\bar{x}_1`` the mean of occasion 1 of the ``n_1`` plots, and ``\\bar{x}_m`` and ``\\bar{y}_m`` the means of occasion 1 and occasion 2 of the ``m`` permanent plots, the regression estimate is (Eq. 11.54)

```math
\\bar{y}_r = \\bar{y}_m + b \\left(\\bar{x}_1 - \\bar{x}_m\\right) \\qquad
b = \\frac{\\sum_{j=1}^{m} (x_j - \\bar{x}_m)(y_j - \\bar{y}_m)}{\\sum_{j=1}^{m} (x_j - \\bar{x}_m)^2}
```

where the slope is solved by the normal equations ``(X^\\top X)\\beta = X^\\top y`` over the permanent plots. Its variance combines the residual variance of the regression, ``s_{yx}^2``, from the ``m`` permanent plots, with the extra precision of the larger sample (Eqs. 11.56 and 11.60):

```math
s_{yx}^2 = \\frac{\\sum_{j=1}^{m} (y_j - \\bar{y}_m)^2 - b \\sum_{j=1}^{m} (x_j - \\bar{x}_m)(y_j - \\bar{y}_m)}{m - 2} \\qquad
s_y^2 = \\frac{\\sum_{j=1}^{m} (y_j - \\bar{y}_m)^2}{m - 1} \\qquad
s_{\\bar{y}_r}^2 = \\frac{s_{yx}^2}{m} + \\frac{s_y^2 - s_{yx}^2}{n_1}
```

The errors use ``t`` of ``m - 1`` degrees of freedom, ``E_a = t \\, s_{\\bar{y}_r}``, ``E_r = E_a / \\bar{y}_r \\times 100`` and the total is ``\\hat{Y} = N \\bar{y}_r`` (Eqs. 11.66 to 11.70). The change and its variance (Eqs. 11.72 and 11.75), with ``r`` the correlation of the ``m`` permanent plots, use ``t`` of ``(n_1 - 1) + (n_2 - 1)`` degrees of freedom, and ``\\hat{D} = N \\bar{d}_r`` (Eqs. 11.77 to 11.81):

```math
\\bar{d}_r = \\bar{y}_r - \\bar{x}_1 \\qquad
s_{\\bar{d}_r}^2 = \\frac{(s_y^2 - s_{yx}^2)\\left[1 + r(r - 2)\\right]}{r^2 \\, n_1} + \\frac{s_{yx}^2}{n_2}
```

# Examples

The example of the book (chapter 11, Quadro 11.3): a stand of 500 ha of *Pinus elliottii* with plots of 1 000 m², so ``N = 500`` plots, and volumes per hectare (the volumes of plots of 1 ha): 56 plots in 1975, 30 temporary and 26 permanent, and the 26 permanent plots remeasured in 1976.

```julia-repl
julia> xu = [132.3, 122.9, 133.3, 124.8, 127.5, 116.4, 120.0, 140.6, 116.8, 113.9, 138.7, 144.5, 104.0,
        157.4, 116.2, 100.1, 164.4, 195.7, 129.8, 102.1, 140.7, 187.0, 64.6, 129.6, 147.3, 151.9,
        110.7, 147.1, 101.8, 146.5];
julia> xm = [139.1, 190.3, 100.7, 102.7, 98.4, 68.2, 42.2, 91.6, 91.6, 77.4, 91.4, 91.6, 123.5, 108.7, 108.7,
        93.6, 42.0, 53.8, 138.0, 173.1, 121.1, 81.8, 64.0, 145.9, 168.9, 186.8];
julia> ym = [177.6, 232.9, 165.5, 180.8, 164.3, 103.6, 68.0, 145.5, 137.6, 123.7, 146.6, 141.3, 174.4, 164.9,
        164.9, 134.0, 57.4, 84.4, 202.0, 207.9, 171.8, 154.0, 101.2, 189.2, 226.9, 239.3];

julia> value1 = vcat(xu, xm) * u"m^3";
julia> value2 = vcat(fill(missing, 30), ym) * u"m^3";

julia> report = doublesampling(value1, value2, 1u"ha", 500);

julia> occasion2(report).xm[1]
169.92370951520402 m^3

julia> occasion2(report).b[1]
1.0971

julia> change(report)
1×7 DataFrame
 Row │ dm           sdm         Ea           Er       D            Dl           Du
     │ Quantity…    Quantity…   Quantity…    Float64  Quantity…    Quantity…    Quantity…
─────┼──────────────────────────────────────────────────────────────────────────────────────
   1 │ 49.8576 m^3  2.7901 m^3  5.55247 m^3    11.14  24928.8 m^3  22152.6 m^3  27705.1 m^3
```

The book gets ``b = 1.097``, ``\\bar{y}_r = 169.952`` m³ with ``s_{yx}^2 = 200.418`` and ``s_{\\bar{y}_r} = 6.664``, and the change ``\\bar{d}_r = 49.86`` m³ with ``s_{\\bar{d}_r}^2 = 7.786``, ``E_a = 5.552`` m³ and ``\\hat{D} = 24\\,930`` m³. The unit 34 of the book is printed with 181.8 m³ in 1976, but its own means and sums need 180.8 m³, the value used here.

# References
Péllico Netto, S. & Brena, D. A. (1997). *Inventário Florestal*. Curitiba, chapter 11.
Cochran, W. G. (1977). *Sampling Techniques*, 3rd ed. New York: Wiley, Sec. 12.9.
"""
function doublesampling(value1::AbstractVector, value2::AbstractVector, plotarea::Area, npop::Integer; α::Real=0.95)

  checklengths(:value1 => value1, :value2 => value2)
  permanentidx = findall(!ismissing, value2)
  m = length(permanentidx)
  n = length(value1)
  m >= 3 || throw(ArgumentError("At least three permanent (remeasured) plots are required."))

  o1 = _occasionmean(value1, npop)
  occasion1 = _occasiontable(o1, npop, plotarea, α)
  vunit = unit(o1.x̅)

  x1perm = ustrip.(value1[permanentidx])
  y2perm = ustrip.(value2[permanentidx])
  x̄1perm, ȳ2perm = mean(x1perm), mean(y2perm)
  x̄1all = ustrip(o1.x̅)

  b = _olsslope(y2perm, x1perm)
  sstY = sum(abs2, y2perm .- ȳ2perm)
  ssxy = sum((x1perm .- x̄1perm) .* (y2perm .- ȳ2perm))
  syx² = (sstY - b * ssxy) / (m - 2)
  sy² = sstY / (m - 1)

  yreg = ȳ2perm + b * (x̄1all - x̄1perm)
  meanvar = syx² / m + (sy² - syx²) / n
  yregq, sx̅ = yreg * vunit, sqrt(meanvar) * vunit

  t = -quantile(TDist(m - 1), (1 - α) / 2)
  absoluteerror = t * sx̅
  relativeerror = absoluteerror / yregq * 100
  hectarevalue = yregq * expansionfactor(plotarea)
  totalvalue = yregq * npop

  occasion2 = DataFrame(
    xm=yregq, b=round(b, digits=4), s2yx=syx² * vunit^2, sxm=sx̅,
    Ea=absoluteerror, Er=round(relativeerror, digits=2), xha=hectarevalue,
    X=totalvalue, Xl=totalvalue - npop * absoluteerror, Xu=totalvalue + npop * absoluteerror,
    m=m, u=n - m, N=npop,
  )

  # the change and its variance (Eqs. 11.72 and 11.75), with r the correlation of the m permanent plots
  r = ssxy / sqrt(sum(abs2, x1perm .- x̄1perm) * sstY)
  r == 0 && throw(ArgumentError("The permanent plots have no correlation between the occasions (r = 0)."))
  growth = yregq - o1.x̅
  growthvar = (sy² - syx²) * (1 + r * (r - 2)) / (r^2 * n) + syx² / m
  change = _changetable(growth, growthvar, (n - 1) + (m - 1), npop, α)

  return SamplingReport((; occasion1, occasion2, change))
end
