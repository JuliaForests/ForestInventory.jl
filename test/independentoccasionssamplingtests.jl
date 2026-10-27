@testset "independentoccasionssampling" begin
  # The example of the book (chapter 11, Quadro 11.1): Pinus elliottii in Guarapuava, PR, with 44
  # independent plots in 1975 and 44 in 1976, in a population of N = 1 500 units. The volumes of the
  # book are per hectare, so they are the volumes of units of 1 ha. The book takes t = 2.021 for
  # each occasion and t = 1.99 for the change, the package the exact ones (2.017 and 1.988).
  v1 = ai_x * u"m^3"
  v2 = ai_y * u"m^3"

  report = independentoccasionssampling(v1, v2, 1u"ha", 1500, 1500)
  @test report isa SamplingReport

  o1 = occasion1(report)
  @test ustrip(o1.xm[1]) ≈ 199.98 atol = 1e-2
  @test ustrip(o1.sxm[1]) ≈ 6.598 atol = 1e-3
  @test ustrip(o1.Ea[1]) ≈ tq(43) * ustrip(o1.sxm[1])
  @test ustrip(o1.Ea[1]) ≈ 13.335 atol = 5e-2
  @test o1.Er[1] ≈ 6.67 atol = 5e-2
  @test ustrip(o1.X[1]) ≈ 299_970 atol = 20
  @test ustrip(o1.Xl[1]) ≈ 279_968 atol = 50
  @test o1.n[1] == 44
  @test o1.N[1] == 1500

  o2 = occasion2(report)
  @test ustrip(o2.xm[1]) ≈ 248.19 atol = 1e-2
  @test ustrip(o2.sxm[1]) ≈ 7.809 atol = 1e-3
  @test ustrip(o2.Ea[1]) ≈ 15.782 atol = 5e-2
  @test o2.Er[1] ≈ 6.36 atol = 5e-2
  @test ustrip(o2.X[1]) ≈ 372_285 atol = 20
  @test o2.n[1] == 44

  # the change: the variance is the sum of the variances of the means (Eq. 11.23), and the
  # degrees of freedom of t are (n1 - 1) + (n2 - 1) = 86 (Eq. 11.27)
  chg = change(report)
  @test ustrip(chg.dm[1]) ≈ 48.21 atol = 1e-2
  @test ustrip(chg.sdm[1]) ≈ 10.22 atol = 1e-2
  @test ustrip(chg.sdm[1])^2 ≈ 104.53 atol = 1e-2
  @test ustrip(chg.Ea[1]) ≈ tq(86) * ustrip(chg.sdm[1])
  @test ustrip(chg.Ea[1]) ≈ 20.34 atol = 5e-2
  @test chg.Er[1] ≈ 42.19 atol = 5e-2
  @test ustrip(chg.D[1]) ≈ 72_315 atol = 5
  @test ustrip(chg.Dl[1]) ≈ 41_808 atol = 50
  @test ustrip(chg.Du[1]) ≈ 102_821 atol = 50

  @testset "units" begin
    reportU = independentoccasionssampling(v1, v2, 1u"ha", 1500, 1500)
    @test change(reportU) == chg
  end

  @testset "no unit assumed for a bare-number input" begin
    reportBare = independentoccasionssampling(ai_x, ai_y, 1u"ha", 1500, 1500)
    @test change(reportBare).dm[1] isa Real
    @test !(change(reportBare).dm[1] isa Unitful.Quantity)
    @test change(reportBare).dm[1] ≈ 48.21 atol = 1e-2
  end
end
