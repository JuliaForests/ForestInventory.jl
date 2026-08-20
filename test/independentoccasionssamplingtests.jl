@testset "independentoccasionssampling" begin
  # Independent-samples successive-occasions exercise from Nardini (2023, PPGAAA/UFSM-FW
  # dissertation underlying this package -- see the package docs "About the Authors"): 44
  # plots at each occasion (0.05 ha plots, 30000 possible plots per occasion). Raw data
  # from the original desktop tool's repository
  # (github.com/claitonnardini/Inventario-Florestal.jl, Dados/ind.csv). Our mean growth
  # matches the dissertation's reported result (48.21 m^3) exactly; the error/CI figures
  # can differ since the dissertation's population sizes for each occasion aren't
  # recorded in its result sheet.
  v1 = [
    146.7, 182.1, 258.6, 198.8, 196.7, 155.6, 201.0, 183.8, 212.4, 202.5, 217.5, 138.9, 242.7,
    200.0, 208.3, 193.0, 279.3, 181.0, 223.4, 210.4, 206.1, 171.2, 264.4, 278.6, 253.3, 41.0,
    242.9, 241.6, 251.2, 144.5, 225.9, 178.5, 182.0, 176.2, 199.5, 156.8, 136.2, 225.4, 184.3,
    245.8, 164.9, 198.8, 168.6, 228.7,
  ]u"m^3"
  v2 = [
    330.4, 222.3, 300.2, 293.8, 270.2, 230.4, 312.3, 341.6, 326.5, 46.3, 312.8, 273.8, 272.3,
    225.2, 282.9, 233.3, 213.4, 224.2, 229.1, 188.0, 176.7, 221.0, 247.5, 249.5, 230.3, 136.4,
    285.4, 265.9, 248.1, 250.2, 238.0, 265.9, 210.6, 280.0, 198.2, 290.3, 248.2, 235.1, 265.2,
    236.7, 280.4, 285.4, 198.3, 248.0,
  ]u"m^3"

  report = independentOccasionsSampling(v1, v2, 0.05, 30000, 30000)
  @test report isa SamplingReport
  chg = change(report)

  @test ustrip(chg.meanGrowth[1]) ≈ 48.2091 atol = 1e-3
  @test ustrip(chg.standardError[1]) ≈ 10.3695 atol = 1e-3
  @test ustrip(chg.totalGrowth[1]) ≈ 1.44627e6 atol = 1e1

  @test ustrip(occasion1(report).mean[1]) ≈ mean(ustrip.(v1)) atol = 1e-9
  @test ustrip(occasion2(report).mean[1]) ≈ mean(ustrip.(v2)) atol = 1e-9
  @test occasion1(report).n[1] == 44
  @test occasion2(report).n[1] == 44

  @testset "units" begin
    reportU = independentOccasionsSampling(v1, v2, 0.05u"ha", 30000, 30000)
    @test change(reportU) == chg
  end

  @testset "no unit assumed for a bare-number input" begin
    reportBare = independentOccasionsSampling(ustrip.(v1), ustrip.(v2), 0.05, 30000, 30000)
    @test change(reportBare).meanGrowth[1] isa Real
    @test !(change(reportBare).meanGrowth[1] isa Unitful.Quantity)
  end
end
