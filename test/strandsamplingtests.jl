@testset "strandsampling" begin
  line = [1, 1, 1, 2, 2, 3]
  diameter = [28.1, 33.2, 31.0, 40.1, 27.9, 34.6]u"cm"
  value = [0.65, 1.01, 0.88, 1.78, 0.60, 1.19]u"m^3"

  report = strandsampling(line, diameter, value, 2.0u"m^2/ha", 20.0u"m", 4, 10u"ha")
  @test report isa SamplingReport
  rt = resulttable(report)
  pt = pointtable(report)

  @test nrow(pt) == 3   # line 4 (zero trees) never appears as a row
  @test rt.n[1] == 4    # but is still counted in the sample size
  @test rt.N[1] == 10
  @test eltype(rt.xha) <: Unitful.Quantity

  @testset "line-1 estimate matches the Kershaw et al. (2017) Table A.6.B derivation" begin
    baf, L = 2.0, 20.0
    k = sqrt(baf / 2500)
    d = [28.1, 33.2, 31.0]
    v = [0.65, 1.01, 0.88]
    tf = 500_000 .* k ./ (d .* L)
    g = (π / 4) .* (d ./ 100) .^ 2

    row1 = pt[1, :]
    @test ustrip(row1.N) ≈ sum(tf) atol = 1e-8
    @test ustrip(row1.G) ≈ sum(tf .* g) atol = 1e-8
    @test ustrip(row1.xha) ≈ sum(tf .* v) atol = 1e-8
  end

  @testset "diameter factor collapses to a per-line constant, independent of Di" begin
    baf, L = 2.0, 20.0
    k = sqrt(baf / 2500)
    d = [28.1, 33.2, 31.0]
    tf = 500_000 .* k ./ (d .* L)
    @test all(x -> isapprox(x, d[1] * tf[1], atol=1e-6), d .* tf)
  end

  @testset "units" begin
    reportU = strandsampling(line, diameter, value, 2.0u"m^2/ha", 20.0u"m", 4, 10u"ha")
    @test resulttable(reportU) == rt
  end

  @testset "order of rows does not matter" begin
    perm = [4, 5, 1, 2, 3, 6]
    reportPerm = strandsampling(line[perm], diameter[perm], value[perm], 2.0u"m^2/ha", 20.0u"m", 4, 10u"ha")
    @test resulttable(reportPerm) == rt
  end

  @testset "zero-count lines are not silently dropped" begin
    reportDropped = strandsampling(line, diameter, value, 2.0u"m^2/ha", 20.0u"m", 5, 10u"ha")
    @test ustrip(resulttable(reportDropped).xha[1]) < ustrip(rt.xha[1])
  end

  @test_throws ArgumentError strandsampling(line, diameter, value, 2.0u"m^2/ha", 20.0u"m", 2, 10u"ha")
  @test_throws ArgumentError strandsampling(line, diameter, value, -2.0u"m^2/ha", 20.0u"m", 4, 10u"ha")
end
