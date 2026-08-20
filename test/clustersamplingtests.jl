@testset "clustersampling" begin
  # Cluster sampling exercise from Nardini (2023, PPGAAA/UFSM-FW dissertation underlying
  # this package -- see the package docs "About the Authors"): 25 clusters ("conglomerados")
  # of 4 plots each (0.25 ha, 4000 ha stand). Raw data from the original desktop tool's
  # repository (github.com/claitonnardini/Inventario-Florestal.jl, Dados/congl.csv). Our
  # mean/total match the dissertation's reported results (256.99 m^3/ha, 1027968 m^3)
  # exactly.
  cluster = repeat(1:25, inner=4)
  value = [
    84.7, 78.2, 63.0, 88.5, 70.1, 65.4, 59.6, 58.1, 75.6, 77.7, 66.1, 43.5, 81.3, 38.7, 62.8, 31.1,
    99.5, 47.0, 82.1, 77.4, 97.3, 66.1, 57.4, 113.9, 70.2, 43.5, 48.3, 69.3, 81.7, 68.8, 22.6, 78.2,
    60.9, 71.4, 75.2, 81.1, 67.3, 90.5, 72.6, 67.1, 75.4, 88.5, 104.6, 105.5, 71.0, 57.0, 79.3, 71.8,
    61.8, 100.9, 69.4, 61.0, 116.9, 91.0, 36.9, 57.3, 56.0, 22.9, 94.8, 130.2, 68.2, 32.5, 68.2, 58.6,
    55.4, 31.3, 76.3, 40.0, 55.0, 63.0, 56.8, 36.8, 78.8, 59.6, 68.1, 46.8, 33.6, 44.4, 34.1, 79.6,
    48.8, 26.6, 67.5, 51.0, 82.4, 61.8, 23.4, 60.2, 69.8, 49.5, 53.9, 0.0, 58.4, 61.5, 45.5, 50.2,
    27.4, 63.6, 69.4, 60.8,
  ]u"m^3"

  report = clusterSampling(cluster, value, 0.25, 4000)
  @test report isa SamplingReport
  rt = resultTable(report)

  @test ustrip(rt.mean[1]) ≈ 64.248 atol = 1e-3
  @test rt.M[1] == 4
  @test rt.n[1] == 25
  @test rt.N[1] == 4000
  @test rt.intraclassCorrelation[1] ≈ 0.1037 atol = 1e-3
  @test eltype(rt.mean) <: Unitful.Quantity

  @testset "M=1 matches simplecasualsampling" begin
    v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1]u"m^3"
    rc = resultTable(clusterSampling(1:11, v, 0.05, 10))
    rs = simpleCasualSampling(v, 0.05, 10)
    @test ustrip(rc.mean[1]) ≈ ustrip(rs.mean[1]) atol = 1e-6
    @test rc.requiredSampleSize[1] == rs.requiredSampleSize[1]
    @test rc.intraclassCorrelation[1] ≈ 1.0 atol = 1e-9
  end

  @testset "units" begin
    reportU = clusterSampling(cluster, value, 0.25u"ha", 4000u"ha")
    @test resultTable(reportU) == rt
  end

  @test_throws ArgumentError clusterSampling([1, 1, 2, 2, 2], [1.0, 2.0, 3.0, 4.0, 5.0]u"m^3", 0.02, 15)
end
