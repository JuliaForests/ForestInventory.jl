@testset "simplecasualsampling" begin
  # Simple random sampling exercise from Nardini (2023, PPGAAA/UFSM-FW dissertation
  # underlying this package -- see the package docs "About the Authors"): 70 plots of
  # 0.1 ha each, 45 ha stand. The raw plot volumes come from the original desktop tool's
  # repository (github.com/claitonnardini/Inventario-Florestal.jl, Dados/aas.csv), itself
  # based on Péllico Netto & Brena's *Inventário Florestal* (1997). Our mean/total match
  # the dissertation's reported results (217.97 m^3/ha, 9808.71 m^3) exactly; other
  # figures (e.g. requiredSampleSize) can differ slightly since this package solves the
  # required-sample-size equation with its own iterative loop.
  v = [
    18.5, 22.8, 15.8, 9.8, 23.5, 20.7, 28.3, 6.9, 28.7, 26.4, 7.6, 10.3, 21.1, 26.3, 39.9,
    12.7, 11.5, 11.0, 18.9, 31.4, 9.6, 9.8, 7.1, 25.9, 25.7, 24.2, 19.0, 25.9, 37.7, 23.3,
    24.6, 32.3, 23.4, 17.1, 8.9, 32.7, 8.1, 29.2, 24.7, 25.6, 13.2, 25.7, 27.0, 20.4, 13.6,
    32.4, 31.3, 34.6, 9.5, 18.7, 28.9, 21.6, 8.3, 23.2, 26.6, 25.6, 8.0, 22.1, 21.9, 23.3,
    9.9, 30.7, 29.2, 28.0, 13.7, 27.6, 37.9, 27.3, 28.0, 30.7,
  ]u"m^3"

  report = simpleCasualSampling(v, 0.1, 45; e=10, α=0.95)
  @test report isa DataFrame
  @test nrow(report) == 1

  @test ustrip(report.mean[1]) ≈ 21.7971 atol = 1e-3
  @test report.coefficientOfVariation[1] ≈ 39.58 atol = 1e-2
  @test ustrip(report.varianceOfMean[1]) ≈ 0.897789 atol = 1e-3
  @test ustrip(report.standardError[1]) ≈ 0.947518 atol = 1e-4
  @test ustrip(report.absoluteError[1]) ≈ 1.89025 atol = 1e-4
  @test report.relativeError[1] ≈ 8.67 atol = 1e-2
  @test ustrip(report.meanPerHectare[1]) ≈ 217.971 atol = 1e-2
  @test ustrip(report.total[1]) ≈ 9808.71 atol = 1e-1
  @test ustrip(report.lowerLimit[1]) ≈ 8958.10 atol = 1e-1
  @test ustrip(report.upperLimit[1]) ≈ 10659.3 atol = 1e-1
  @test report.population[1] == "finite"
  @test report.samplingFraction[1] ≈ 0.844 atol = 1e-3
  @test report.n[1] == 70
  @test report.requiredSampleSize[1] == 56
  @test report.additionalSampleSize[1] == 0
  @test report.N[1] == 450

  # every column that carries a physical unit is a real Unitful quantity column, not a
  # plain number paired with a separate units column -- this is what makes it compatible
  # with `removeunits`/`restoreunits` (ForestFoundations.jl)
  @test eltype(report.mean) <: Unitful.Quantity
  @test eltype(report.total) <: Unitful.Quantity
  @test eltype(report.coefficientOfVariation) <: Real
  stripped = removeunits(report)
  @test "mean (m^3)" ∈ names(stripped)
  @test stripped[1, "mean (m^3)"] ≈ 21.7971 atol = 1e-3
  restored = restoreunits(stripped)
  @test restored.mean == report.mean

  @testset "units are propagated, never assumed" begin
    reportU = simpleCasualSampling(v, 0.1u"ha", 45u"ha")
    @test reportU == report

    vft = ustrip.(v) * u"ft^3"
    reportFt = simpleCasualSampling(vft, 0.1u"ac", 45u"ac")
    @test unit(reportFt.mean[1]) == u"ft^3"
    @test unit(reportFt.meanPerHectare[1]) == u"ft^3/ac"
  end

  @testset "no unit assumed for a bare-number input" begin
    vbare = ustrip.(v)
    reportBare = simpleCasualSampling(vbare, 0.1, 45; e=10, α=0.95)
    @test !(eltype(reportBare.mean) <: Unitful.Quantity)
    @test eltype(reportBare.mean) <: Real
    @test reportBare.mean[1] ≈ 21.7971 atol = 1e-3
    # plot_area/total_area still default to hectares regardless of value's unit
    @test reportBare.N[1] == 450
  end

  @testset "infinite population via 1-reference-area plot" begin
    reportInf = simpleCasualSampling(v, 1.0, 1000)
    @test reportInf.population[1] == "infinite"
  end

  @test_throws ArgumentError simpleCasualSampling([1.0]u"m^3", 0.05, 10)
end
