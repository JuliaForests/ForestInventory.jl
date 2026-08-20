@testset "horizontalpointsampling" begin
  point = [1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 5, 5, 5]
  diameter = [25.0, 30.0, 20.0, 28.0, 22.0, 35.0, 30.0, 25.0, 20.0, 22.0, 30.0, 28.0, 26.0]u"cm"
  value = [0.35, 0.55, 0.22, 0.48, 0.28, 0.85, 0.55, 0.35, 0.22, 0.28, 0.55, 0.48, 0.40]u"m^3"

  report = horizontalPointSampling(point, diameter, value, 2.0, 6, 0.1, 10)
  @test report isa SamplingReport
  rt = resultTable(report)
  pt = pointTable(report)

  @test nrow(pt) == 5   # point 6 (zero trees) never appears as a row
  @test rt.n[1] == 6    # but is still counted in the sample size
  @test rt.N[1] == 100
  @test eltype(rt.meanPerHectare) <: Unitful.Quantity

  @testset "basalAreaPerHectare identity: EFi*gi sums to treeCount*baf exactly" begin
    baf = 2.0
    for row in eachrow(pt)
      @test ustrip(row.basalAreaPerHectare) ≈ row.treeCount * baf atol = 1e-9
    end
  end

  @testset "units" begin
    reportU = horizontalPointSampling(point, diameter, value, 2.0, 6, 0.1u"ha", 10u"ha")
    @test resultTable(reportU) == rt
  end

  @testset "diameter still defaults to cm when given bare, unlike value" begin
    reportBareDiam = horizontalPointSampling(point, ustrip.(diameter), value, 2.0, 6, 0.1, 10)
    @test resultTable(reportBareDiam) == rt
  end

  @testset "zero-count points are not silently dropped" begin
    reportDropped = horizontalPointSampling(point, diameter, value, 2.0, 5, 0.1, 10)
    @test ustrip(resultTable(reportDropped).meanPerHectare[1]) > ustrip(rt.meanPerHectare[1])
  end

  @test_throws ArgumentError horizontalPointSampling(point, diameter, value, 2.0, 4, 0.1, 10)
  @test_throws ArgumentError horizontalPointSampling(point, diameter, value, -2.0, 6, 0.1, 10)
end
