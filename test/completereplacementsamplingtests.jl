@testset "completereplacementsampling" begin
  # The example of the book (chapter 11, Quadro 11.2): Eucalyptus sp. in Lençóis Paulista, SP,
  # with 45 permanent plots measured in 1976 and 1977, in a population of N = 1 500 units (the
  # volumes of the book are per hectare, so they are the volumes of units of 1 ha). The book takes
  # t = 2.021 for each occasion and t = 1.99 for the change, the package the exact ones.
  v1 = art_x * u"m^3"
  v2 = art_y * u"m^3"

  report = completereplacementsampling(v1, v2, 1u"ha", 1500, 1500)
  @test report isa SamplingReport

  o1 = occasion1(report)
  @test ustrip(o1.xm[1]) ≈ 121.15 atol = 1e-2
  @test ustrip(o1.sxm[1])^2 ≈ 33.14 atol = 1e-2
  o2 = occasion2(report)
  @test ustrip(o2.xm[1]) ≈ 160.34 atol = 1e-2
  @test ustrip(o2.sxm[1]) ≈ 7.124 atol = 1e-3
  @test ustrip(o2.Ea[1]) ≈ 14.397 atol = 5e-2
  @test o2.Er[1] ≈ 8.98 atol = 5e-2
  @test ustrip(o2.X[1]) ≈ 240_512 atol = 20

  # the variance of the change subtracts the covariance (Eq. 11.35): 33.14 + 50.75 - 2(1 807.46)/45
  chg = change(report)
  @test ustrip(chg.dm[1]) ≈ 39.19 atol = 1e-2
  @test ustrip(chg.sdm[1])^2 ≈ 3.548 atol = 1e-3
  @test ustrip(chg.sdm[1]) ≈ 1.884 atol = 1e-3
  @test ustrip(chg.Ea[1]) ≈ tq(88) * ustrip(chg.sdm[1])    # (n1 - 1) + (n2 - 1) = 88 degrees of freedom
  @test ustrip(chg.Ea[1]) ≈ 3.749 atol = 1e-2
  @test chg.Er[1] ≈ 9.56 atol = 2e-2
  @test ustrip(chg.D[1]) ≈ 58_784 atol = 2
  @test ustrip(chg.Dl[1]) ≈ 53_161 atol = 15
  @test ustrip(chg.Du[1]) ≈ 64_407 atol = 15

  @test_throws DimensionMismatch completereplacementsampling(v1, [1.0, 2.0]u"m^3", 1u"ha", 1500, 1500)

  @testset "units" begin
    reportU = completereplacementsampling(v1, v2, 1u"ha", 1500, 1500)
    @test change(reportU) == chg
  end
end
