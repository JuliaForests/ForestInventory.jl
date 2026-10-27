@testset "multistartsystematicsampling" begin
  # Systematic sampling with multiple random starts exercise from Nardini (2023,
  # PPGAAA/UFSM-FW dissertation underlying this package -- see the package docs "About
  # the Authors"): 9 random starts of 8 plots each (0.25 ha, 13000 ha stand). Raw data
  # from the original desktop tool's repository
  # (github.com/claitonnardini/Inventario-Florestal.jl, Dados/mult.csv). Our mean/total
  # match the dissertation's reported results (199.08 m^3/ha, 2588062 m^3) exactly.
  # The example of the book (chapter 10, Quadro 10.1): 9 random starts (clusters) of 8 plots of
  # 0.25 ha, in a stand of 13 000 ha, so N = 6 500 possible starts. The book takes t = 1.996 for the
  # errors, the package the t of the 8 degrees of freedom of the starts (2.306).
  start = mul_start
  value = mul_value * u"m^3"

  report = multistartsystematicsampling(start, value, 0.25u"ha", 13000u"ha")
  @test report isa SamplingReport
  rt = resulttable(report)

  @test ustrip(rt.xm[1]) ≈ 49.77 atol = 1e-2
  @test ustrip(rt.s2d[1]) ≈ 283.169 atol = 1e-2
  @test ustrip(rt.s2e[1]) ≈ 27.487 atol = 1e-2
  @test rt.r[1] ≈ 0.088 atol = 1e-3
  @test ustrip(rt.s2xm[1]) ≈ 6.987 atol = 1e-3
  @test ustrip(rt.sxm[1]) ≈ 2.643 atol = 1e-3
  @test ustrip(rt.Ea[1]) ≈ tq(8) * ustrip(rt.sxm[1])
  @test ustrip(rt.X[1]) ≈ 2_587_988 atol = 100    # the book: 6 500 × 8 × 49.769
  @test rt.nreq[1] == 13                           # the book: 13.49 → 11.84 → 12.29 starts
  @test rt.M[1] == 8
  @test rt.n[1] == 9
  @test rt.N[1] == 6500
  @test nrow(starttable(report)) == 9
  @test ustrip(starttable(report).xmi[1]) ≈ 54.569 atol = 1e-3

  @testset "units" begin
    reportU = multistartsystematicsampling(start, value, 0.25u"ha", 13000u"ha")
    @test resulttable(reportU) == rt
  end
end
