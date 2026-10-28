@testset "systematicsampling" begin
  # Systematic sampling exercise from Nardini (2023, PPGAAA/UFSM-FW dissertation
  # underlying this package -- see the package docs "About the Authors"): 5 systematic
  # lines of 15 plots each (0.1 ha, 45 ha stand). Raw data from the original desktop
  # tool's repository (github.com/claitonnardini/Inventario-Florestal.jl,
  # Dados/sist.csv). Our mean/total match the dissertation's reported results (221.11
  # m^3/ha, 9949.92 m^3) exactly.
  # The example of the book (chapter 7, Quadro 7.1): 75 plots of 0.1 ha in a stand of 45 ha, in 15
  # lines (strips) of 5 plots, which the vector holds in the order of the line. The book uses
  # t = 2.0, the package the exact t of 74 degrees of freedom (1.993).
  v = sys_value * u"m^3"
  line = sys_line

  @test systematicsampling(v, 0.1u"ha", 45u"ha"; line=line) isa SamplingReport

  report = systematicsampling(v, 0.1u"ha", 45u"ha"; line=line).resulttable
  @test nrow(report) == 1

  @test ustrip(report.xm[1]) ≈ 22.11 atol = 1e-2
  @test ustrip(report.s2xm[1]) ≈ 0.1394 atol = 1e-4
  @test ustrip(report.sxm[1]) ≈ 0.3734 atol = 1e-4
  @test ustrip(report.Ea[1]) ≈ tq(74) * ustrip(report.sxm[1])
  @test ustrip(report.Ea[1]) ≈ 0.7468 atol = 1e-2
  @test report.Er[1] ≈ 3.38 atol = 2e-2
  @test ustrip(report.X[1]) ≈ 9950 atol = 1
  @test ustrip(report.Xl[1]) ≈ 9614 atol = 2
  @test ustrip(report.Xu[1]) ≈ 10286 atol = 3
  @test report.cv[1] ≈ std(sys_value) / mean(sys_value) * 100 atol = 1e-2
  @test report.m[1] == 15
  @test report.n[1] == 75
  @test report.N[1] == 450

  @testset "with lines" begin
    reportNoLine = systematicsampling(v, 0.1u"ha", 45u"ha").resulttable
    @test reportNoLine.m[1] == 1
    # the differences straddling each pair of lines are excluded, so the variance differs
    # from the single-line (no grouping) case
    @test ustrip(reportNoLine.s2xm[1]) != ustrip(report.s2xm[1])
    # a single line: the classic Σ(xi - xi+1)² / [2n(n - 1)] of the book (Eq. 7.11), with the correction 1 - f
    d = sum(abs2, diff(sys_value))
    @test ustrip(reportNoLine.s2xm[1]) ≈ d / (2 * 75 * 74) * (1 - 75 / 450)
  end

  @testset "units" begin
    reportU = systematicsampling(v, 0.1u"ha", 45u"ha"; line=line).resulttable
    @test reportU == report
  end

  @testset "no unit assumed for a bare-number input" begin
    reportBare = systematicsampling(ustrip.(v), 0.1u"ha", 45u"ha"; line=line).resulttable
    @test reportBare.xm[1] isa Real
    @test !(reportBare.xm[1] isa Unitful.Quantity)
  end

  @test_throws DimensionMismatch systematicsampling(v, 0.1u"ha", 45u"ha"; line=[1, 2, 3])
  @test_throws ArgumentError systematicsampling([1.0]u"m^3", 0.05u"ha", 10u"ha")
end
