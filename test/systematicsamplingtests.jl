@testset "systematicsampling" begin
  # Systematic sampling exercise from Nardini (2023, PPGAAA/UFSM-FW dissertation
  # underlying this package -- see the package docs "About the Authors"): 5 systematic
  # lines of 15 plots each (0.1 ha, 45 ha stand). Raw data from the original desktop
  # tool's repository (github.com/claitonnardini/Inventario-Florestal.jl,
  # Dados/sist.csv). Our mean/total match the dissertation's reported results (221.11
  # m^3/ha, 9949.92 m^3) exactly.
  line1 = [8.0, 8.6, 13.1, 16.6, 21.6, 23.6, 19.7, 30.6, 20.4, 28.0, 30.1, 22.6, 26.7, 31.8, 25.5]
  line2 = [9.4, 12.7, 7.6, 19.0, 15.1, 23.7, 18.4, 29.4, 24.6, 38.6, 27.6, 24.7, 24.8, 28.7, 33.1]
  line3 = [7.3, 9.8, 9.5, 18.6, 15.8, 23.3, 22.8, 27.8, 19.2, 28.3, 23.1, 26.7, 26.9, 33.9, 35.1]
  line4 = [10.1, 11.6, 8.5, 10.8, 18.0, 23.0, 27.1, 28.7, 28.2, 34.9, 28.4, 31.7, 31.8, 28.3, 30.1]
  line5 = [8.7, 8.8, 10.4, 22.9, 16.22, 24.0, 19.5, 30.5, 19.9, 22.9, 21.4, 28.7, 21.8, 30.7, 28.3]
  v = vcat(line1, line2, line3, line4, line5)u"m^3"
  line = repeat(1:5, inner=15)

  report = systematicSampling(v, 0.1, 45; line=line)
  @test nrow(report) == 1

  @test ustrip(report.mean[1]) ≈ 22.1109 atol = 1e-3
  @test ustrip(report.varianceOfMean[1]) ≈ 0.21084 atol = 1e-3
  @test ustrip(report.standardError[1]) ≈ 0.459173 atol = 1e-3
  @test report.samplingInterval[1] == 5
  @test report.n[1] == 75
  @test report.N[1] == 450

  @testset "with lines" begin
    reportNoLine = systematicSampling(v, 0.1, 45)
    @test reportNoLine.samplingInterval[1] == 1
    # the differences straddling each pair of lines are excluded, so the variance differs
    # from the single-line (no grouping) case
    @test ustrip(reportNoLine.varianceOfMean[1]) != ustrip(report.varianceOfMean[1])
  end

  @testset "units" begin
    reportU = systematicSampling(v, 0.1u"ha", 45u"ha"; line=line)
    @test reportU == report
  end

  @testset "no unit assumed for a bare-number input" begin
    reportBare = systematicSampling(ustrip.(v), 0.1, 45; line=line)
    @test reportBare.mean[1] isa Real
    @test !(reportBare.mean[1] isa Unitful.Quantity)
  end

  @test_throws DimensionMismatch systematicSampling(v, 0.1, 45; line=[1, 2, 3])
  @test_throws ArgumentError systematicSampling([1.0]u"m^3", 0.05, 10)
end
