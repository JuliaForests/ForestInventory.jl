@testset "clustersampling" begin
  # Cluster sampling exercise from Nardini (2023, PPGAAA/UFSM-FW dissertation underlying
  # this package -- see the package docs "About the Authors"): 25 clusters ("conglomerados")
  # of 4 plots each (0.25 ha, 4000 ha stand). Raw data from the original desktop tool's
  # repository (github.com/claitonnardini/Inventario-Florestal.jl, Dados/congl.csv). Our
  # mean/total match the dissertation's reported results (256.99 m^3/ha, 1027968 m^3)
  # exactly.
  # The example of the book (chapter 9, Quadro 9.1): 25 clusters of 4 plots of 0.25 ha, in a stand
  # of 4 000 ha (N = 4 000 possible clusters). The book takes t = 1.99 for the errors, the package
  # the t of the 24 degrees of freedom of the clusters (2.064).
  cluster = clu_cluster
  value = clu_value * u"m^3"

  report = clustersampling(cluster, value, 0.25u"ha", 4000u"ha")
  @test report isa SamplingReport
  rt = resulttable(report)

  @test ustrip(rt.xm[1]) ≈ 64.548 atol = 1e-3
  @test ustrip(rt.s2d[1]) ≈ 432.58 atol = 1e-2
  @test ustrip(rt.s2e[1]) ≈ 46.25 atol = 1e-2
  @test ustrip(rt.s2x[1]) ≈ 478.83 atol = 1e-2
  @test rt.r[1] ≈ 0.0966 atol = 1e-4
  @test ustrip(rt.s2xm[1]) ≈ 6.176 atol = 1e-3
  @test ustrip(rt.sxm[1]) ≈ 2.483 atol = 3e-3
  @test ustrip(rt.Ea[1]) ≈ tq(24) * ustrip(rt.sxm[1])
  @test ustrip(rt.X[1]) ≈ 1_032_768 atol = 1
  @test rt.population[1] == "infinite"    # f = 25/4000 is negligible
  @test rt.nreq[1] == 17                  # the book: 15.79 → 16.83 → 16.65 clusters
  @test rt.M[1] == 4
  @test rt.n[1] == 25
  @test rt.N[1] == 4000
  @test eltype(rt.xm) <: Unitful.Quantity
  @test nrow(clustertable(report)) == 25
  @test ustrip(clustertable(report).xmi[1]) ≈ 78.6 atol = 1e-3

  @testset "the correction multiplies only the variance between clusters (Eq. 9.17)" begin
    # with a population of 100 clusters the sampling fraction 25/100 is not negligible
    r100 = resulttable(clustersampling(cluster, value, 0.25u"ha", 100u"ha"))
    @test r100.population[1] == "finite"
    s2e, s2d = 46.249385, 432.584733
    @test ustrip(r100.s2xm[1]) ≈ (100 - 25) / 100 * s2e / 25 + s2d / (25 * 4) atol = 1e-4
  end

  @testset "M=1 matches simplecasualsampling" begin
    v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1]u"m^3"
    rc = resulttable(clustersampling(1:11, v, 0.05u"ha", 10u"ha"))
    rs = simplecasualsampling(v, 0.05u"ha", 10u"ha")
    @test ustrip(rc.xm[1]) ≈ ustrip(rs.xm[1]) atol = 1e-6
    @test rc.nreq[1] == rs.nreq[1]
    @test rc.r[1] ≈ 1.0 atol = 1e-9
  end

  @testset "units" begin
    reportU = clustersampling(cluster, value, 0.25u"ha", 4000u"ha")
    @test resulttable(reportU) == rt
  end

  @test_throws ArgumentError clustersampling([1, 1, 2, 2, 2], [1.0, 2.0, 3.0, 4.0, 5.0]u"m^3", 0.02u"ha", 15u"ha")
end
