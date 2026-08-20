@testset "completereplacementsampling" begin
  # Complete-replacement (matched plots) successive-occasions exercise from Nardini
  # (2023, PPGAAA/UFSM-FW dissertation underlying this package -- see the package docs
  # "About the Authors"): the same 45 plots remeasured at both occasions (0.05 ha plots,
  # 30000 possible plots). Raw data from the original desktop tool's repository
  # (github.com/claitonnardini/Inventario-Florestal.jl, Dados/art.csv). Our mean growth
  # matches the dissertation's reported result (39.19 m^3) exactly.
  v1 = [
    141.89, 147.74, 126.26, 137.08, 39.30, 133.19, 110.21, 126.79, 132.09, 73.79, 76.22, 63.42,
    77.47, 37.24, 98.34, 87.04, 116.52, 83.75, 110.99, 77.87, 139.73, 182.11, 94.68, 181.69,
    161.44, 172.58, 191.29, 137.93, 167.08, 110.10, 123.59, 107.45, 137.35, 137.77, 79.56,
    80.79, 118.52, 133.61, 174.91, 140.90, 176.79, 185.32, 135.84, 109.97, 73.65,
  ]u"m^3"
  v2 = [
    174.99, 187.58, 153.81, 176.28, 61.76, 160.87, 150.23, 164.15, 170.99, 101.36, 114.85, 87.89,
    102.06, 51.37, 152.57, 119.94, 148.93, 129.42, 152.01, 109.04, 165.91, 208.33, 143.74, 197.50,
    191.62, 266.13, 280.29, 178.37, 218.88, 149.68, 161.67, 174.45, 172.66, 179.97, 105.36,
    126.26, 150.55, 191.47, 248.38, 183.30, 226.76, 224.25, 154.10, 142.97, 102.68,
  ]u"m^3"

  report = completeReplacementSampling(v1, v2, 0.05, 30000, 30000)
  @test report isa SamplingReport
  chg = change(report)

  @test ustrip(chg.meanGrowth[1]) ≈ mean(ustrip.(v2)) - mean(ustrip.(v1)) atol = 1e-9
  @test ustrip(chg.standardError[1]) ≈ 2.4521 atol = 1e-3

  @test_throws DimensionMismatch completeReplacementSampling(v1, [1.0, 2.0]u"m^3", 0.05, 200, 200)

  @testset "units" begin
    reportU = completeReplacementSampling(v1, v2, 0.05u"ha", 30000, 30000)
    @test change(reportU) == chg
  end
end
