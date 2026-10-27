@testset "doublesampling" begin
  # The example of the book (chapter 11, Quadro 11.3): a stand of 500 ha of Pinus elliottii, with
  # circular plots of 1 000 m², so N = 500 plots of 1 ha (the volumes are per hectare): n1 = 56
  # plots in 1975 (u = 30 temporary and m = 26 permanent) and the 26 permanent plots remeasured in
  # 1976. The book takes t = 2.06 for the regression and t = 1.99 for the change.
  value1 = vcat(da_xu, da_xm) * u"m^3"
  value2 = vcat(fill(missing, 30), da_ym) * u"m^3"

  report = doublesampling(value1, value2, 1u"ha", 500)
  @test report isa SamplingReport

  # the first occasion is a simple random sample of the 56 plots (Eqs. 11.43 to 11.51)
  o1 = occasion1(report)
  @test ustrip(o1.xm[1]) ≈ 120.07 atol = 1e-2
  @test ustrip(o1.sxm[1])^2 ≈ 20.372 atol = 1e-3
  @test ustrip(o1.sxm[1]) ≈ 4.514 atol = 1e-3
  @test ustrip(o1.X[1]) ≈ 60_035 atol = 5

  # the second occasion by regression (Eqs. 11.54 to 11.60): the book b = 1.097, s2yx = 200.418
  o2 = occasion2(report)
  @test ustrip(o2.xm[1]) ≈ 169.95 atol = 5e-2
  @test o2.b[1] ≈ 1.0971 atol = 1e-3
  @test ustrip(o2.s2yx[1]) ≈ 200.418 atol = 1e-2
  @test ustrip(o2.sxm[1])^2 ≈ 44.417 atol = 5e-2
  @test ustrip(o2.sxm[1]) ≈ 6.664 atol = 5e-3
  @test ustrip(o2.Ea[1]) ≈ tq(25) * ustrip(o2.sxm[1])
  @test ustrip(o2.Ea[1]) ≈ 13.728 atol = 5e-2
  @test o2.Er[1] ≈ 8.08 atol = 2e-2
  @test ustrip(o2.X[1]) ≈ 84_976 atol = 20
  @test o2.m[1] == 26
  @test o2.u[1] == 30
  @test o2.N[1] == 500

  # the change (Eqs. 11.72 and 11.75): the book d = 49.86, s2d = 7.786, with (n1 - 1) + (n2 - 1) = 80
  chg = change(report)
  @test ustrip(chg.dm[1]) ≈ 49.86 atol = 5e-2
  @test ustrip(chg.sdm[1])^2 ≈ 7.786 atol = 2e-2
  @test ustrip(chg.sdm[1]) ≈ 2.790 atol = 2e-3
  @test ustrip(chg.Ea[1]) ≈ tq(80) * ustrip(chg.sdm[1])
  @test ustrip(chg.Ea[1]) ≈ 5.552 atol = 2e-2
  @test chg.Er[1] ≈ 11.14 atol = 2e-2
  @test ustrip(chg.D[1]) ≈ 24_930 atol = 20
  @test ustrip(chg.Dl[1]) ≈ 22_153 atol = 20
  @test ustrip(chg.Du[1]) ≈ 27_706 atol = 20

  @test_throws ArgumentError doublesampling(value1, [fill(missing, 55)..., 226.9]u"m^3", 1u"ha", 500)
  @test_throws DimensionMismatch doublesampling(value1, value2[1:(end-1)], 1u"ha", 500)

  @testset "units" begin
    reportU = doublesampling(value1, value2, 1u"ha", 500)
    @test occasion2(reportU) == o2
  end

  @testset "no unit assumed for a bare-number input" begin
    reportBare = doublesampling(ustrip.(value1), ustrip.(value2), 1u"ha", 500)
    @test occasion2(reportBare).xm[1] isa Real
    @test !(occasion2(reportBare).xm[1] isa Unitful.Quantity)
    @test occasion2(reportBare).xm[1] ≈ 169.95 atol = 5e-2
  end
end
