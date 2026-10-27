@testset "horizontalpointsampling" begin
  point = [1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 5, 5, 5]
  diameter = [25.0, 30.0, 20.0, 28.0, 22.0, 35.0, 30.0, 25.0, 20.0, 22.0, 30.0, 28.0, 26.0]u"cm"
  value = [0.35, 0.55, 0.22, 0.48, 0.28, 0.85, 0.55, 0.35, 0.22, 0.28, 0.55, 0.48, 0.40]u"m^3"

  report = horizontalpointsampling(point, diameter, value, 2.0u"m^2/ha", 6, 0.1u"ha", 10u"ha")
  @test report isa SamplingReport
  rt = resulttable(report)
  pt = pointtable(report)

  @test nrow(pt) == 5   # point 6 (zero trees) never appears as a row
  @test rt.n[1] == 6    # but is still counted in the sample size
  @test rt.N[1] == 100
  @test eltype(rt.xha) <: Unitful.Quantity

  @testset "G identity: EFi*gi sums to n*baf exactly" begin
    baf = 2.0
    for row in eachrow(pt)
      @test ustrip(row.G) ≈ row.n * baf atol = 1e-9
    end
  end

  @testset "units" begin
    reportU = horizontalpointsampling(point, diameter, value, 2.0u"m^2/ha", 6, 0.1u"ha", 10u"ha")
    @test resulttable(reportU) == rt
  end

  @testset "zero-count points are not silently dropped" begin
    reportDropped = horizontalpointsampling(point, diameter, value, 2.0u"m^2/ha", 5, 0.1u"ha", 10u"ha")
    @test ustrip(resulttable(reportDropped).xha[1]) > ustrip(rt.xha[1])
  end

  @test_throws ArgumentError horizontalpointsampling(point, diameter, value, 2.0u"m^2/ha", 4, 0.1u"ha", 10u"ha")
  @test_throws ArgumentError horizontalpointsampling(point, diameter, value, -2.0u"m^2/ha", 6, 0.1u"ha", 10u"ha")
end
