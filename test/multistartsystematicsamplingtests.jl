@testset "multistartsystematicsampling" begin
  # Systematic sampling with multiple random starts exercise from Nardini (2023,
  # PPGAAA/UFSM-FW dissertation underlying this package -- see the package docs "About
  # the Authors"): 9 random starts of 8 plots each (0.25 ha, 13000 ha stand). Raw data
  # from the original desktop tool's repository
  # (github.com/claitonnardini/Inventario-Florestal.jl, Dados/mult.csv). Our mean/total
  # match the dissertation's reported results (199.08 m^3/ha, 2588062 m^3) exactly.
  start = repeat(1:9, inner=8)
  value = [
    25.97, 34.62, 25.97, 60.27, 49.75, 81.64, 85.90, 72.43,
    54.02, 56.77, 63.78, 61.94, 61.28, 38.93, 44.68, 48.30,
    66.43, 50.46, 43.81, 39.63, 42.47, 38.84, 37.13, 35.93,
    59.37, 56.18, 39.01, 32.93, 95.67, 32.82, 65.09, 46.29,
    62.69, 68.88, 65.84, 85.93, 48.13, 55.98, 70.91, 37.90,
    47.85, 69.66, 78.11, 36.71, 28.09, 53.57, 40.34, 56.22,
    36.64, 39.23, 42.81, 35.51, 27.54, 54.78, 42.56, 88.45,
    34.72, 41.21, 41.79, 34.97, 41.47, 29.11, 25.60, 21.99,
    43.13, 59.13, 92.98, 49.14, 46.68, 42.60, 26.17, 30.14,
  ]u"m^3"

  report = multistartSystematicSampling(start, value, 0.25, 13000)
  @test report isa SamplingReport
  rt = resultTable(report)

  @test ustrip(rt.mean[1]) ≈ 49.7704 atol = 1e-3
  @test rt.M[1] == 8
  @test rt.n[1] == 9
  @test rt.N[1] == 6500

  @testset "units" begin
    reportU = multistartSystematicSampling(start, value, 0.25u"ha", 13000u"ha")
    @test resultTable(reportU) == rt
  end
end
