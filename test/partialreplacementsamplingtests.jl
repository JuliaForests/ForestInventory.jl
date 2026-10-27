@testset "partialreplacementsampling" begin
  # The example of the book (chapter 11, Quadro 11.6): a stand of Eucalyptus sp. in Lençóis Paulista,
  # SP, with N = 1 582 plots (volumes per hectare, so plots of 1 ha): u = 17 plots measured only in
  # 1976, m = 13 permanent plots and n = 13 new plots measured only in 1977. The estimators are those
  # of Ware & Cunia (1962), Eqs. 11.85 to 11.109 of the book.
  value1 = vcat(arp_xu, arp_xm, fill(missing, 13)) * u"m^3"
  value2 = vcat(fill(missing, 17), arp_ym, arp_yn) * u"m^3"

  report = partialreplacementsampling(value1, value2, 1u"ha", 1582)
  @test report isa SamplingReport

  # the first occasion: the 30 plots, with t = 2.045 (29 degrees of freedom)
  o1 = occasion1(report)
  @test ustrip(o1.xm[1]) ≈ 123.02 atol = 1e-2
  @test ustrip(o1.sxm[1])^2 ≈ 41.28 atol = 1e-2
  @test ustrip(o1.sxm[1]) ≈ 6.425 atol = 1e-3
  @test ustrip(o1.Ea[1]) ≈ 13.139 atol = 5e-3
  @test o1.Er[1] ≈ 10.68 atol = 1e-2
  @test ustrip(o1.X[1]) ≈ 194_617 atol = 15

  # the second occasion: the current mean, y = 167.92, with s = 7.568 and t = 2.06 (25 degrees of freedom)
  o2 = occasion2(report)
  @test o2.r[1] ≈ 0.869 atol = 1e-3
  @test o2.c[1] ≈ 0.6362 atol = 5e-4
  @test o2.a[1] ≈ 0.4029 atol = 5e-4
  @test ustrip(o2.xm[1]) ≈ 167.92 atol = 2e-2
  @test ustrip(o2.sxm[1]) ≈ 7.568 atol = 1e-2
  @test ustrip(o2.Ea[1]) ≈ tq(25) * ustrip(o2.sxm[1])
  @test ustrip(o2.Ea[1]) ≈ 15.590 atol = 2e-2
  @test o2.Er[1] ≈ 9.28 atol = 2e-2
  @test ustrip(o2.X[1]) ≈ 265_649 atol = 30
  @test o2.m[1] == 13
  @test o2.u[1] == 17
  @test o2.n[1] == 13
  @test o2.N[1] == 1582

  # the change, estimated directly (Eqs. 11.104 to 11.107): d = 44.32, s2d = 21.99, with
  # (n1 - 1) + (n2 - 1) = 54 degrees of freedom
  chg = change(report)
  @test ustrip(chg.dm[1]) ≈ 44.32 atol = 1e-2
  @test ustrip(chg.sdm[1])^2 ≈ 21.99 atol = 1e-2
  @test ustrip(chg.sdm[1]) ≈ 4.689 atol = 2e-3
  @test ustrip(chg.Ea[1]) ≈ tq(54) * ustrip(chg.sdm[1])
  @test ustrip(chg.Ea[1]) ≈ 9.378 atol = 3e-2
  @test chg.Er[1] ≈ 21.16 atol = 8e-2
  @test ustrip(chg.D[1]) ≈ 70_114 atol = 10

  @test_throws DimensionMismatch partialreplacementsampling(value1, value2[1:(end-1)], 1u"ha", 1582)
  @test_throws ArgumentError partialreplacementsampling(
    [1.0, 2.0, missing, missing]u"m^3", [missing, missing, 3.0, 4.0]u"m^3", 1u"ha", 200)

  @testset "units" begin
    reportU = partialreplacementsampling(value1, value2, 1u"ha", 1582)
    @test occasion2(reportU) == o2
  end

  @testset "no unit assumed for a bare-number input" begin
    reportBare = partialreplacementsampling(ustrip.(value1), ustrip.(value2), 1u"ha", 1582)
    @test occasion2(reportBare).xm[1] isa Real
    @test !(occasion2(reportBare).xm[1] isa Unitful.Quantity)
    @test occasion2(reportBare).xm[1] ≈ 167.92 atol = 2e-2
  end
end
