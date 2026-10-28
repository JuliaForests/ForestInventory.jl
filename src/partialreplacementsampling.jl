"""
    partialreplacementsampling(value1::AbstractVector, value2::AbstractVector, plotarea::Area, npop::Integer; α::Real=0.95)

Estimates the mean, the total and the change (growth) of a forest between two inventory occasions by sampling with partial replacement of the plots, with the estimators of Ware & Cunia (1962), generic to any sample size.

# Description

A middle ground between [`completereplacementsampling`](@ref) (remeasure every plot) and [`doublesampling`](@ref) (remeasure a subsample for a regression). On occasion 1 there are ``u`` temporary plots and ``m`` permanent plots (``n_1 = u + m``); on occasion 2 the ``m`` permanent plots are remeasured and ``n`` new temporary plots are added (``n_2 = m + n``). Three groups of plots appear in the two vectors:

- measured on occasion 1 only (``u`` plots: `value1[i]` present, `value2[i]` missing);
- measured on both occasions (``m`` permanent plots: both present);
- measured on occasion 2 only (``n`` new temporary plots: `value1[i]` missing, `value2[i]` present).

The strong correlation between the volumes of the permanent plots on the two occasions makes a linear regression possible, which carries to the temporary plots the information of the occasion in which they were not measured, so every estimate uses all the ``u + m + n`` plots. With ``m = 0`` the design becomes independent samples, with ``u = n = 0`` complete replacement and with ``n = 0`` double sampling. The notation and the formulas follow Péllico Netto & Brena (1997), chapter 11, Sec. 11.4.

# Arguments

- `value1`, `value2`: the measurements of occasion 1 and occasion 2, of the same length, one entry per plot, `missing` where the plot was not measured on that occasion. A `Unitful` quantity keeps its unit; a plain number stays unitless.
- `plotarea`: the area of each plot, an `Area`, used only for the mean per hectare.
- `npop::Integer`: the possible plots of the population.
- `α::Real=0.95`: the confidence level.

# Returns

A [`SamplingReport`](@ref) with `occasion1` (the simple random sample of the ``n_1 = u + m`` plots of occasion 1, with the columns of [`independentoccasionssampling`](@ref), Eqs. 11.43 to 11.51), `occasion2` (the current mean) and `change`. The columns of `occasion2`:

- `xm`: the current mean per plot of occasion 2, ``\\bar{y}`` (Eq. 11.85).
- `sxm`, `Ea`, `Er`, `xha`, `X`, `Xl`, `Xu`: as in [`simplecasualsampling`](@ref), for ``\\bar{y}`` (``t`` of ``n_2 - 1`` degrees of freedom).
- `a`, `c`: the constants of the estimator. `r`: the correlation of the permanent plots.
- `m`, `u`, `n`: the permanent, the occasion-1-only and the new temporary plots (here `n` is the number of new plots, as in the book). `N`: the possible plots.

`change` has the columns of [`independentoccasionssampling`](@ref), for the direct estimate of the change ``\\bar{d}_p`` (Eq. 11.104), with ``t`` of ``(n_1 - 1) + (n_2 - 1)`` degrees of freedom.

# Mathematical basis

With ``P_u = u/n_1`` and ``P_m = m/n_1``, ``\\bar{x}_u`` and ``\\bar{x}_m`` the means of occasion 1 of the temporary and of the permanent plots, ``\\bar{y}_m`` and ``\\bar{y}_n`` the means of occasion 2 of the permanent and of the new plots, ``s_x^2`` the variance of the ``n_1`` plots of occasion 1 and ``s_y^2`` that of the ``n_2`` plots of occasion 2, and ``s_{x_m}``, ``s_{y_m}`` and ``r`` the standard deviations and the correlation of the ``m`` permanent plots (Eqs. 11.63, 11.64), the regression coefficient is ``b_{yx} = r s_{y_m}/s_{x_m}`` (Eq. 11.89), and the constants (Eqs. 11.87 and 11.88) are

```math
c = \\frac{m}{n_2 - P_u n r^2} \\qquad
a = \\frac{m P_u}{n_2 - P_u n r^2} \\, b_{yx}
```

The current mean (the best estimate of the mean of occasion 2, Eq. 11.85), its variance (Eq. 11.86), the standard error, the errors, the total and its interval (Eqs. 11.94 to 11.97) are

```math
\\bar{y} = a \\bar{x}_u - a \\bar{x}_m + c \\bar{y}_m + (1 - c) \\bar{y}_n \\qquad
s_{\\bar{y}}^2 = a^2 s_x^2 \\left(\\frac{1}{u} + \\frac{1}{m}\\right) + c^2 \\frac{s_y^2}{m} + (1 - c)^2 \\frac{s_y^2}{n} - 2 a c \\, r \\frac{s_{x_m} s_{y_m}}{m}
```

```math
E_a = t \\, s_{\\bar{y}} \\qquad E_r = \\frac{t \\, s_{\\bar{y}}}{\\bar{y}} \\times 100 \\qquad
\\hat{Y} = N \\bar{y} \\qquad
IC\\left[\\hat{Y} - N t s_{\\bar{y}} \\le Y \\le \\hat{Y} + N t s_{\\bar{y}}\\right] = P
```

The change is estimated directly (Eqs. 11.103 to 11.107), from the coefficients ``A`` and ``B`` that minimize its variance, with the degrees of freedom ``(n_1 - 1) + (n_2 - 1)`` for ``t``:

```math
A = \\frac{m}{n_2 - P_u n r^2} + \\frac{n P_m}{n_2 - P_u n r^2} \\, r \\frac{s_{x_m}}{s_{y_m}} \\qquad
B = - \\frac{m P_u}{n_2 - P_u n r^2} \\, r \\frac{s_{y_m}}{s_{x_m}} - \\frac{n_2 P_m}{n_2 - P_u n r^2}
```

```math
\\bar{d}_p = A \\bar{y}_m + (1 - A) \\bar{y}_n + B \\bar{x}_m - (1 + B) \\bar{x}_u
```

```math
s_{\\bar{d}_p}^2 = A^2 \\frac{s_y^2}{m} + (1 - A)^2 \\frac{s_y^2}{n} + B^2 \\frac{s_x^2}{m} + (1 + B)^2 \\frac{s_x^2}{u} + 2 A B \\, r \\frac{s_{x_m} s_{y_m}}{m}
```

The errors, the total change ``\\hat{D} = N \\bar{d}_p`` and its interval follow Eqs. 11.77 to 11.81.

# Technical description

The book also gives the best estimate of the mean of occasion 1, ``\\bar{x}`` (Eqs. 11.100 to 11.102), from which the change ``\\bar{y} - \\bar{x}`` follows; it is the same number as ``\\bar{d}_p``. The required sampling intensity of the book's Sec. 11.4 is not computed here.

# Examples

The example of the book (chapter 11, Quadro 11.6): a stand of *Eucalyptus* sp. of 1 582 plots, volumes in 1976 and 1977, with 17 plots measured only in 1976, 13 permanent plots and 13 new plots measured only in 1977. The volumes of the book are per hectare, so they are given as the volumes of plots of 1 ha.

```julia-repl
julia> xu = [141.9, 147.7, 126.3, 137.1, 133.2, 110.2, 126.8, 132.1, 73.8, 76.2, 63.4, 77.4, 98.3, 87.4,
        116.5, 83.8, 111.0];
julia> xm = [77.9, 139.7, 182.1, 94.7, 181.7, 161.4, 172.6, 191.3, 137.9, 167.1, 110.1, 123.6, 107.5];
julia> ym = [109.0, 165.9, 208.3, 143.7, 197.5, 191.6, 266.1, 280.3, 178.4, 218.9, 149.7, 161.7, 174.5];
julia> yn = [172.7, 180.0, 105.4, 126.3, 150.6, 191.5, 248.4, 183.3, 226.8, 224.3, 154.1, 142.9, 102.7];

julia> value1 = vcat(xu, xm, fill(missing, 13)) * u"m^3";
julia> value2 = vcat(fill(missing, 17), ym, yn) * u"m^3";

julia> report = partialreplacementsampling(value1, value2, 1u"ha", 1582);

julia> occasion2(report).xm[1]
167.93032593719232 m^3

julia> change(report).dm[1]
44.32321388247928 m^3
```

The book gets ``\\bar{y} = 167.92`` m³ (``s_{\\bar{y}} = 7.568``), ``\\bar{d}_p = 44.32`` m³ (``s_{\\bar{d}_p} = 4.689``) and ``\\hat{D} = 70\\,114`` m³; the differences are the rounding of the intermediate values.

# References
Péllico Netto, S. & Brena, D. A. (1997). *Inventário Florestal*. Curitiba, chapter 11.
Ware, K. D. & Cunia, T. (1962). *Continuous forest inventory with partial replacement of samples*. Forest Science Monograph 3.
"""
function partialreplacementsampling(value1::AbstractVector, value2::AbstractVector, plotarea::Area, npop::Integer; α::Real=0.95)

  checklengths(:value1 => value1, :value2 => value2)

  idx = eachindex(value1)
  matched = findall(i -> !ismissing(value1[i]) && !ismissing(value2[i]), idx)
  unmatched = findall(i -> !ismissing(value1[i]) && ismissing(value2[i]), idx)
  newtemp = findall(i -> ismissing(value1[i]) && !ismissing(value2[i]), idx)
  m, u, n = length(matched), length(unmatched), length(newtemp)
  m >= 3 || throw(ArgumentError("At least three matched (remeasured on both occasions) plots are required."))
  u >= 2 || throw(ArgumentError("At least two occasion-1-only temporary plots are required."))
  n >= 2 || throw(ArgumentError("At least two new temporary (occasion-2-only) plots are required."))

  n1, n2 = u + m, m + n
  o1 = _occasionmean(collect(skipmissing(value1)), npop)
  occasion1 = _occasiontable(o1, npop, plotarea, α)
  vunit = unit(o1.x̅)

  # the four subsamples, without unit, and the proportions of the permanent and temporary plots
  xu, xm = ustrip.(value1[unmatched]), ustrip.(value1[matched])
  ym, yn = ustrip.(value2[matched]), ustrip.(value2[newtemp])
  Pu, Pm = u / n1, m / n1
  # the variances of the plots of each occasion, and the standard deviations and the correlation of the permanent plots
  s²x, s²y = var(vcat(xu, xm)), var(vcat(ym, yn))
  sx, sy, r = std(xm), std(ym), cor(xm, ym)

  # the current mean (Eqs. 11.85, 11.87, 11.88 and 11.89) and its variance (Eq. 11.86)
  den = n2 - Pu * n * r^2
  byx = r * sy / sx
  c = m / den
  a = m * Pu / den * byx
  ȳ = a * mean(xu) - a * mean(xm) + c * mean(ym) + (1 - c) * mean(yn)
  s²ȳ = a^2 * s²x * (1 / u + 1 / m) + c^2 * s²y / m + (1 - c)^2 * s²y / n - 2 * a * c * r * sx * sy / m

  ȳq, sȳ = ȳ * vunit, sqrt(s²ȳ) * vunit
  t = -quantile(TDist(n2 - 1), (1 - α) / 2)
  absoluteerror = t * sȳ
  relativeerror = absoluteerror / ȳq * 100
  hectarevalue = ȳq * expansionfactor(plotarea)
  totalvalue = ȳq * npop

  occasion2 = DataFrame(
    xm=ȳq, sxm=sȳ, Ea=absoluteerror, Er=round(relativeerror, digits=2), xha=hectarevalue,
    X=totalvalue, Xl=totalvalue - npop * absoluteerror, Xu=totalvalue + npop * absoluteerror,
    a=round(a, digits=4), c=round(c, digits=4), r=round(r, digits=4), m=m, u=u, n=n, N=npop,
  )

  # the direct estimate of the change and its variance (Eqs. 11.104 to 11.107)
  A = m / den + n * Pm / den * r * sx / sy
  B = -m * Pu / den * r * sy / sx - n2 * Pm / den
  d = A * mean(ym) + (1 - A) * mean(yn) + B * mean(xm) - (1 + B) * mean(xu)
  s²d = A^2 * s²y / m + (1 - A)^2 * s²y / n + B^2 * s²x / m + (1 + B)^2 * s²x / u + 2 * A * B * r * sx * sy / m
  change = _changetable(d * vunit, s²d, (n1 - 1) + (n2 - 1), npop, α)

  return SamplingReport((; occasion1, occasion2, change))
end
