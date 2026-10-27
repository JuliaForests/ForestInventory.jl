@testset "prodansampling" begin
  point = [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
  diameter = [30.56, 30.40, 36.29, 58.25, 29.60, 36.29,
    28.10, 33.20, 31.05, 40.12, 27.88, 34.60]u"cm"
  value = [0.8108, 0.7982, 1.3397, 4.8380, 0.7367, 1.3397,
    0.6534, 1.0102, 0.8811, 1.7765, 0.6011, 1.1893]u"m^3"
  distance = [2.1, 2.8, 3.4, 3.9, 4.3, 5.2, 1.9, 2.6, 3.1, 3.7, 4.1, 4.9]u"m"

  report = prodansampling(point, diameter, value, distance, 10u"ha")
  @test report isa SamplingReport
  rt = resulttable(report)
  pt = pointtable(report)

  @test nrow(pt) == 2
  @test rt.k[1] == 6
  @test rt.n[1] == 2
  @test eltype(rt.N) <: Unitful.Quantity

  @testset "point-1 estimate matches Prodan's (k - 0.5) hand-derivation" begin
    L = 5.2   # farthest (6th) tree's distance, in meters
    circlearea_ha = (π * L^2) / 10_000   # m^2 -> ha
    d = [30.56, 30.40, 36.29, 58.25, 29.60, 36.29] ./ 100   # m
    g = (π / 4) .* d .^ 2                                   # m^2
    v = [0.8108, 0.7982, 1.3397, 4.8380, 0.7367, 1.3397]

    n_expected = 5.5 / circlearea_ha
    g_expected = (sum(g[1:5]) + g[6] / 2) / circlearea_ha
    v_expected = (sum(v[1:5]) + v[6] / 2) / circlearea_ha

    row1 = pt[1, :]
    @test ustrip(row1.N) ≈ n_expected atol = 1e-6
    @test ustrip(row1.G) ≈ g_expected atol = 1e-6
    @test ustrip(row1.xha) ≈ v_expected atol = 1e-6
  end

  @testset "result is the plain arithmetic mean across points" begin
    @test ustrip(rt.N[1]) ≈ mean(ustrip.(pt.N)) atol = 1e-9
    @test ustrip(rt.G[1]) ≈ mean(ustrip.(pt.G)) atol = 1e-9
    @test rt.Gt[1] ≈ rt.G[1] * 10u"ha"
    @test rt.X[1] ≈ rt.xha[1] * 10u"ha"
  end

  @testset "trees need not be pre-sorted by distance" begin
    shuffled = [4, 6, 1, 2, 3, 5]   # point 1's own rows, permuted
    reportShuffled = prodansampling(
      point[shuffled], diameter[shuffled], value[shuffled], distance[shuffled], 10u"ha",
    )
    @test resulttable(reportShuffled) == resulttable(prodansampling(point[1:6], diameter[1:6], value[1:6], distance[1:6], 10u"ha"))
  end

  @testset "units" begin
    reportU = prodansampling(point, diameter, value, distance, 10u"ha")
    @test resulttable(reportU) == rt
  end

  @test_throws ArgumentError prodansampling([1, 1, 2, 2, 2, 2, 2, 2], diameter[1:8], value[1:8], distance[1:8], 10u"ha")
  @test_throws ArgumentError prodansampling([1], diameter[1:1], value[1:1], distance[1:1], 10u"ha")
  @test_throws DimensionMismatch prodansampling(point, diameter, value, distance[1:(end-1)], 10u"ha")
end
