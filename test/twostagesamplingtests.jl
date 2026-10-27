@testset "twostagesampling" begin
  # Two-stage sampling exercise from Nardini (2023, PPGAAA/UFSM-FW dissertation
  # underlying this package -- see the package docs "About the Authors"): 36 primary
  # units sub-sampled with 6 secondary units each (0.25 ha plots, out of N=100 possible
  # primaries and M=40000 possible secondary units per primary). Raw data from the
  # original desktop tool's repository
  # (github.com/claitonnardini/Inventario-Florestal.jl, Dados/de.csv). Our mean/total
  # match the dissertation's reported results (223.51 m^3/ha, 223509259 m^3) exactly.
  # The example of the book (chapter 8): the stand of 1 000 000 ha has N = 100 primary units of
  # 10 000 ha, each with M = 40 000 plots of 0.25 ha. The pilot of Quadro 8.2 has 18 primary units
  # of 6 plots, and the definitive inventory of Quadro 8.3 has 36 (the first 18 are the pilot).
  primary = two_primary
  value = two_value * u"m^3"

  @testset "the pilot inventory of the book" begin
    pilot = resulttable(twostagesampling(primary[1:108], value[1:108], 0.25u"ha", 100, 40000))
    @test ustrip(pilot.xm[1]) ≈ 64.740 atol = 1e-3
    @test ustrip(pilot.s2d[1]) ≈ 697.231 atol = 1e-3    # MQd
    @test ustrip(pilot.s2e[1]) ≈ 389.505 atol = 1e-3
    @test pilot.nreq[1] == 36    # the book: 37.996 → 35.876 → 36.004 primary units
  end

  report = twostagesampling(primary, value, 0.25u"ha", 100, 40000)
  @test report isa SamplingReport
  rt = resulttable(report)

  # In the book MQe is printed as 2 218.968, but its own means give 2 193.009, so its s2e (280.555),
  # s2xm (7.467) and s_xm (2.733) are slightly off; the package gives the ones of the book's data.
  @test ustrip(rt.xm[1]) ≈ 55.877 atol = 1e-3
  @test ustrip(rt.s2d[1]) ≈ 535.637 atol = 1e-3
  @test ustrip(rt.s2e[1]) ≈ 276.228 atol = 1e-3
  @test ustrip(rt.s2xm[1]) ≈ 7.3905 atol = 1e-3
  @test ustrip(rt.sxm[1]) ≈ 2.7186 atol = 1e-3
  @test ustrip(rt.Ea[1]) ≈ tq(35) * ustrip(rt.sxm[1])
  @test ustrip(rt.X[1]) ≈ 223_508_000 rtol = 1e-4
  @test rt.m[1] == 6
  @test rt.M[1] == 40000
  @test rt.n[1] == 36
  @test rt.N[1] == 100
  @test nrow(primarytable(report)) == 36

  @testset "units" begin
    reportU = twostagesampling(primary, value, 0.25u"ha", 100, 40000)
    @test resulttable(reportU) == rt
  end

  @test_throws ArgumentError twostagesampling(primary, value, 0.02u"ha", 40, 2)
end
