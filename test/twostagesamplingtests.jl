@testset "twostagesampling" begin
  # Two-stage sampling exercise from Nardini (2023, PPGAAA/UFSM-FW dissertation
  # underlying this package -- see the package docs "About the Authors"): 36 primary
  # units sub-sampled with 6 secondary units each (0.25 ha plots, out of N=100 possible
  # primaries and M=40000 possible secondary units per primary). Raw data from the
  # original desktop tool's repository
  # (github.com/claitonnardini/Inventario-Florestal.jl, Dados/de.csv). Our mean/total
  # match the dissertation's reported results (223.51 m^3/ha, 223509259 m^3) exactly.
  primary = repeat(1:36, inner=6)
  value = [
    29.07, 72.37, 96.03, 55.90, 106.91, 54.51,
    113.55, 53.14, 63.08, 125.64, 113.10, 27.21,
    37.26, 58.57, 40.27, 58.81, 57.30, 72.97,
    52.66, 58.97, 37.86, 74.61, 69.56, 70.67,
    60.81, 37.76, 63.43, 75.40, 46.24, 44.19,
    72.98, 65.00, 94.62, 55.52, 57.81, 48.97,
    120.71, 73.83, 108.75, 47.76, 101.33, 98.52,
    62.48, 70.30, 41.59, 49.85, 65.64, 65.05,
    73.35, 38.46, 138.26, 106.78, 154.13, 145.69,
    56.35, 19.51, 20.25, 23.54, 38.34, 25.06,
    73.65, 51.87, 62.15, 43.96, 30.82, 41.77,
    65.58, 36.48, 38.76, 26.37, 83.33, 47.14,
    34.31, 29.64, 30.69, 23.87, 23.63, 37.50,
    45.98, 50.52, 44.55, 72.07, 42.39, 86.99,
    55.58, 67.72, 44.15, 107.90, 51.05, 101.61,
    26.79, 27.88, 69.14, 36.53, 41.43, 55.25,
    52.11, 47.06, 80.57, 20.37, 166.67, 118.35,
    44.40, 117.04, 82.26, 130.69, 123.53, 131.59,
    41.48, 90.94, 60.43, 47.60, 58.71, 44.05,
    59.69, 65.41, 72.89, 55.52, 41.42, 76.90,
    50.27, 52.25, 31.06, 79.50, 88.10, 60.19,
    62.24, 27.18, 69.97, 72.08, 60.97, 46.94,
    36.71, 31.28, 34.94, 37.25, 64.99, 40.83,
    63.36, 26.13, 25.89, 33.80, 25.12, 60.89,
    43.08, 41.36, 30.66, 32.78, 28.10, 30.42,
    38.29, 37.51, 29.58, 62.67, 56.40, 73.09,
    32.34, 57.92, 89.36, 16.09, 52.05, 54.40,
    75.36, 59.68, 36.97, 22.02, 40.25, 69.87,
    53.63, 21.01, 46.76, 84.00, 88.49, 45.35,
    48.85, 45.06, 30.19, 30.51, 27.35, 48.72,
    21.05, 79.97, 20.08, 62.15, 16.20, 56.69,
    29.13, 27.72, 35.92, 80.50, 98.08, 25.65,
    24.77, 20.32, 62.24, 27.18, 69.97, 72.08,
    47.21, 43.56, 25.74, 23.53, 28.29, 50.16,
    43.11, 40.80, 53.10, 28.08, 54.64, 19.82,
    40.25, 33.02, 32.24, 23.87, 64.41, 18.85,
  ]u"m^3"

  report = twoStageSampling(primary, value, 0.25, 100, 40000)
  @test report isa SamplingReport
  rt = resultTable(report)

  @test ustrip(rt.mean[1]) ≈ 55.8773 atol = 1e-3
  @test rt.m[1] == 6
  @test rt.M[1] == 40000
  @test rt.n[1] == 36
  @test rt.N[1] == 100

  @testset "units" begin
    reportU = twoStageSampling(primary, value, 0.25u"ha", 100, 40000)
    @test resultTable(reportU) == rt
  end

  @test_throws ArgumentError twoStageSampling(primary, value, 0.02, 40, 2)
end
