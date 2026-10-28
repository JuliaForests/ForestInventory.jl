@testset "simplecasualsampling" begin
  # The example of the book (chapter 5): a population of 450 plots of 0.1 ha (45 ha) of *Pinus*
  # sp.; the pilot inventory of 20 plots (Quadro 5.1) and the definitive one of 70 plots
  # (Quadro 5.2). The book uses the tabulated t = 1.99 for 69 degrees of freedom, the package
  # the exact one (1.995), so the errors and the limits differ in the third digit.
  v = srs_def * u"m^3"

  report = simplecasualsampling(v, 0.1u"ha", 45u"ha"; e=10, α=0.95)
  @test report isa DataFrame
  @test nrow(report) == 1

  @test ustrip(report.xm[1]) ≈ 21.797 atol = 1e-3
  @test report.cv[1] ≈ 39.58 atol = 1e-2
  @test ustrip(report.s2xm[1]) ≈ 0.897 atol = 1e-3
  @test ustrip(report.sxm[1]) ≈ 0.947 atol = 1e-3
  @test ustrip(report.Ea[1]) ≈ tq(69) * ustrip(report.sxm[1])
  @test ustrip(report.Ea[1]) ≈ 1.885 atol = 1e-2    # the book, with t = 1.99
  @test report.Er[1] ≈ 8.65 atol = 3e-2
  @test ustrip(report.xha[1]) ≈ 217.97 atol = 1e-2
  @test ustrip(report.X[1]) ≈ 9808.71 atol = 1e-1
  @test ustrip(report.Xl[1]) ≈ 8960 atol = 3
  @test ustrip(report.Xu[1]) ≈ 10656 atol = 5
  @test report.population[1] == "finite"
  @test report.f[1] ≈ 70 / 450 atol = 1e-3
  @test report.n[1] == 70
  @test report.nreq[1] == 56
  @test report.nadd[1] == 0
  @test report.N[1] == 450

  @testset "the pilot inventory of the book" begin
    # 20 plots: mean 19.605, cv 45.64%, and the required number of plots of the book, 70
    pilot = simplecasualsampling(srs_pilot * u"m^3", 0.1u"ha", 45u"ha"; e=10, α=0.95)
    @test ustrip(pilot.xm[1]) ≈ 19.605 atol = 1e-3
    @test pilot.cv[1] ≈ 45.64 atol = 1e-2
    @test pilot.nreq[1] == 70
    @test pilot.nadd[1] == 50
  end

  # every column that carries a physical unit is a real Unitful quantity column, not a
  # plain number paired with a separate units column -- this is what makes it compatible
  # with `stripunits`/`parseunits` (ForestFoundations.jl)
  @test eltype(report.xm) <: Unitful.Quantity
  @test eltype(report.X) <: Unitful.Quantity
  @test eltype(report.cv) <: Real
  stripped = stripunits(report)
  @test "xm (m^3)" ∈ names(stripped)
  @test stripped[1, "xm (m^3)"] ≈ 21.797 atol = 1e-3
  restored = parseunits(stripped)
  @test restored.xm == report.xm

  @testset "units are propagated, never assumed" begin
    reportU = simplecasualsampling(v, 0.1u"ha", 45u"ha")
    @test reportU == report

    vft = ustrip.(v) * u"ft^3"
    reportFt = simplecasualsampling(vft, 0.1u"ac", 45u"ac")
    @test unit(reportFt.xm[1]) == u"ft^3"
    @test unit(reportFt.xha[1]) == u"ft^3/ac"
  end

  @testset "no unit assumed for a bare-number input" begin
    vbare = ustrip.(v)
    reportBare = simplecasualsampling(vbare, 0.1u"ha", 45u"ha"; e=10, α=0.95)
    @test !(eltype(reportBare.xm) <: Unitful.Quantity)
    @test eltype(reportBare.xm) <: Real
    @test reportBare.xm[1] ≈ 21.797 atol = 1e-3
    # plotarea/totalarea still default to hectares regardless of value's unit
    @test reportBare.N[1] == 450
  end

  @testset "a count of trees, with no unit, is per hectare" begin
    # the number of trees of each plot has no physical unit: the mean per hectare is in 1/ha
    trees = [23, 31, 27, 19, 35, 28, 30, 25, 22, 29]
    rt = simplecasualsampling(trees, 0.05u"ha", 10u"ha")
    @test rt.xm[1] ≈ 26.9
    @test unit(rt.xha[1]) == u"ha^-1"
    @test ustrip(rt.xha[1]) ≈ 538
    @test rt.X[1] ≈ 5380    # trees in the 10 ha
  end

  @testset "the population is infinite only by the sampling fraction" begin
    # 70 plots of 1 ha in 1000 ha: 1 - f = 0.93 < 0.98, finite, whatever the size of the plot
    reportFin = simplecasualsampling(v, 1.0u"ha", 1000u"ha")
    @test reportFin.population[1] == "finite"
    # 70 plots in 10 000 ha: 1 - f = 0.993 >= 0.98, infinite
    reportInf = simplecasualsampling(v, 1.0u"ha", 10000u"ha")
    @test reportInf.population[1] == "infinite"
  end

  @test_throws ArgumentError simplecasualsampling([1.0]u"m^3", 0.05u"ha", 10u"ha")
end
