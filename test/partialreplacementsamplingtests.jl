@testset "partialreplacementsampling" begin
  # Partial-replacement successive-occasions exercise from Nardini (2023, PPGAAA/UFSM-FW
  # dissertation underlying this package -- see the package docs "About the Authors"):
  # 17 plots measured only at the first occasion, 13 matched plots remeasured at both,
  # and 14 new plots measured only at the second (0.05 ha plots, 10000 possible plots).
  # Raw data from the original desktop tool's repository
  # (github.com/claitonnardini/Inventario-Florestal.jl, Dados/arp.csv). Our mean growth
  # (44.99 m^3) is close to, but not identical to, the dissertation's reported result
  # (44.32 m^3/ha) -- expected, since this package solves the regression-weighted
  # estimator with its own iterative sample-size loop.
  value1 = [
    141.9, 147.7, 126.3, 137.1, 133.2, 110.2, 126.8, 132.1, 73.8, 76.2, 63.4, 77.4, 98.3, 87.4,
    116.5, 83.8, 111.0,
    77.9, 139.7, 182.1, 94.7, 181.7, 161.4, 172.6, 191.3, 137.9, 167.1, 110.1, 123.6, 107.5,
    fill(missing, 13)...,
  ]u"m^3"
  value2 = [
    fill(missing, 17)...,
    109.0, 165.9, 208.3, 143.7, 197.5, 191.6, 266.1, 280.3, 178.4, 218.9, 149.7, 161.7, 174.5,
    172.7, 180.0, 105.4, 126.3, 150.6, 191.5, 248.4, 183.3, 226.8, 224.3, 154.1, 142.9, 102.7,
  ]u"m^3"

  report = partialReplacementSampling(value1, value2, 0.05, 10000)
  @test report isa SamplingReport
  o2 = occasion2(report)
  chg = change(report)

  @test ustrip(o2.mean[1]) ≈ 168.008 atol = 1e-3
  @test o2.m[1] == 13
  @test o2.u[1] == 17
  @test o2.v[1] == 13
  @test 0 <= o2.regressionWeight[1] <= 1
  @test ustrip(chg.meanGrowth[1]) ≈ 44.9851 atol = 1e-3

  @test_throws DimensionMismatch partialReplacementSampling(value1, value2[1:end-1], 0.05, 10000)
  @test_throws ArgumentError partialReplacementSampling(
    [1.0, 2.0, missing, missing]u"m^3", [missing, missing, 3.0, 4.0]u"m^3", 0.05, 200)

  @testset "units" begin
    reportU = partialReplacementSampling(value1, value2, 0.05u"ha", 10000)
    @test occasion2(reportU) == o2
  end

  @testset "no unit assumed for a bare-number input" begin
    reportBare = partialReplacementSampling(ustrip.(value1), ustrip.(value2), 0.05, 10000)
    @test occasion2(reportBare).mean[1] isa Real
    @test !(occasion2(reportBare).mean[1] isa Unitful.Quantity)
    @test occasion2(reportBare).mean[1] ≈ 168.008 atol = 1e-3
  end
end
