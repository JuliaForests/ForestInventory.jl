@testset "doublesampling" begin
  # Double sampling exercise from Nardini (2023, PPGAAA/UFSM-FW dissertation underlying
  # this package -- see the package docs "About the Authors"): 56 plots measured at the
  # first occasion, only the last 26 ("permanent" subset) remeasured at the second (0.05
  # ha plots, 10000 possible plots). Raw data from the original desktop tool's
  # repository (github.com/claitonnardini/Inventario-Florestal.jl, Dados/ad.csv). Our
  # mean growth matches the dissertation's reported result (49.89 m^3) exactly.
  value1 = [
    132.3, 122.9, 133.3, 124.8, 127.5, 116.4, 120.0, 140.6, 116.8, 113.9, 138.7, 144.5, 104.0,
    157.4, 116.2, 100.1, 164.4, 195.7, 129.8, 102.1, 140.7, 187.0, 64.6, 129.6, 147.3, 151.9,
    110.7, 147.1, 101.8, 146.5,
    139.1, 190.3, 100.7, 102.7, 98.4, 68.2, 42.2, 91.6, 91.6, 77.4, 91.4, 91.6, 123.5, 108.7,
    108.7, 93.6, 42.0, 53.8, 138.0, 173.1, 121.1, 81.8, 64.0, 145.9, 168.9, 186.8,
  ]u"m^3"
  value2 = [
    fill(missing, 30)...,
    177.6, 232.9, 165.5, 181.8, 164.3, 103.6, 68.0, 145.5, 137.6, 123.7, 146.6, 141.3, 174.4,
    164.9, 164.9, 134.0, 57.4, 84.4, 202.0, 207.9, 171.8, 154.0, 101.2, 189.2, 226.9, 239.3,
  ]u"m^3"

  report = doubleSampling(value1, value2, 0.05, 10000)
  @test report isa SamplingReport
  o2 = occasion2(report)
  chg = change(report)

  @test ustrip(o2.mean[1]) ≈ 169.961 atol = 1e-3
  @test o2.regressionSlope[1] ≈ 1.0969 atol = 1e-3
  @test o2.m[1] == 26
  @test o2.temporaryCount[1] == 30
  @test ustrip(chg.meanGrowth[1]) ≈ 49.8947 atol = 1e-3

  @test_throws ArgumentError doubleSampling(value1,
    [fill(missing, 55)..., 226.9]u"m^3", 0.05, 10000)
  @test_throws DimensionMismatch doubleSampling(value1, value2[1:end-1], 0.05, 10000)

  @testset "units" begin
    reportU = doubleSampling(value1, value2, 0.05u"ha", 10000)
    @test occasion2(reportU) == o2
  end

  @testset "no unit assumed for a bare-number input" begin
    reportBare = doubleSampling(ustrip.(value1), ustrip.(value2), 0.05, 10000)
    @test occasion2(reportBare).mean[1] isa Real
    @test !(occasion2(reportBare).mean[1] isa Unitful.Quantity)
    @test occasion2(reportBare).mean[1] ≈ 169.961 atol = 1e-3
  end
end
