@testset "stratifiedsampling" begin
  # Stratified sampling exercise from Nardini (2023, PPGAAA/UFSM-FW dissertation
  # underlying this package -- see the package docs "About the Authors"): 3 strata (7, 8,
  # and 7 plots, 0.1 ha each) in a 45 ha stand. Raw data from the original desktop tool's
  # repository (github.com/claitonnardini/Inventario-Florestal.jl, Dados/estrat.csv). The
  # dissertation's spreadsheet doesn't record each stratum's individual area, only the
  # 45 ha total, so it's split evenly (15 ha each) here -- consequently our stratified
  # mean/total (213.49 m^3/ha, 9607.23 m^3) differ somewhat from the dissertation's
  # (214.95 m^3/ha, 9672.95 m^3), which used the book's actual per-stratum weights.
  stratum = vcat(fill(1, 7), fill(2, 8), fill(3, 7))
  value = [
    15.8, 7.6, 8.8, 12.5, 11.1, 16.2, 12.2,
    20.4, 30.5, 26.8, 27.2, 28.4, 19.7, 20.4, 23.1,
    21.3, 24.3, 29.2, 21.8, 33.1, 35.8, 26.7,
  ]u"m^3"
  strata_area = Dict(1 => 15.0, 2 => 15.0, 3 => 15.0)

  report = stratifiedSampling(stratum, value, 0.1, strata_area)
  @test report isa SamplingReport
  rt = resultTable(report)
  @test nrow(rt) == 1

  @test ustrip(rt.mean[1]) ≈ 21.3494 atol = 1e-3
  @test ustrip(rt.total[1]) ≈ 9607.23 atol = 1e-1
  @test rt.sampleSizeByStratum[1] == (7, 8, 7)
  @test rt.n[1] == 22
  @test rt.N[1] == 450
  @test eltype(rt.mean) <: Unitful.Quantity

  @test nrow(auxiliaryTable(report)) == 3
  @test auxiliaryTable(report).n == [7, 8, 7]
  @test ustrip(auxiliaryTable(report).mean) ≈ [12.0286, 24.5625, 27.4571] atol = 1e-3

  @test nrow(anova(report)) == 3
  @test anova(report).DOF == [2, 19, 21]

  @testset "strata given out of sorted order still match the right area" begin
    # rows collected in the field are rarely already sorted by stratum -- the dict-based
    # strata_area must resolve by key, not by position, so shuffling stratum/value
    # (together) must not change the result at all.
    shuffled = [15, 3, 20, 8, 1, 22, 11, 5, 17, 9, 2, 19, 6, 14, 21, 10, 4, 18, 7, 13, 16, 12]  # a permutation of 1:22
    rtShuffled = resultTable(stratifiedSampling(stratum[shuffled], value[shuffled], 0.1, strata_area))
    # summing the same floats in a different order can differ in the last bit, so compare
    # numerically rather than requiring exact equality
    @test ustrip(rtShuffled.mean[1]) ≈ ustrip(rt.mean[1])
    @test ustrip(rtShuffled.total[1]) ≈ ustrip(rt.total[1])
    @test rtShuffled.sampleSizeByStratum[1] == rt.sampleSizeByStratum[1]
    @test rtShuffled.N[1] == rt.N[1]
  end

  @testset "one stratum matches simplecasualsampling" begin
    v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1]u"m^3"
    stratum1 = fill(1, 11)
    r1 = resultTable(stratifiedSampling(stratum1, v, 0.05, Dict(1 => 10.0)))
    r2 = simpleCasualSampling(v, 0.05, 10)
    @test ustrip(r1.mean[1]) ≈ ustrip(r2.mean[1]) atol = 1e-6
    @test ustrip(r1.varianceOfMean[1]) ≈ ustrip(r2.varianceOfMean[1]) atol = 1e-6
    @test ustrip(r1.standardError[1]) ≈ ustrip(r2.standardError[1]) atol = 1e-6
    @test r1.requiredSampleSize[1] == r2.requiredSampleSize[1]
    @test r1.population[1] == r2.population[1]
    @test r1.samplingFraction[1] ≈ r2.samplingFraction[1] atol = 1e-6
    # the ANOVA between-strata line is undefined with a single stratum (0 numerator DOF)
    r1full = stratifiedSampling(stratum1, v, 0.05, Dict(1 => 10.0))
    @test ismissing(anova(r1full).F[1])
  end

  @testset "units" begin
    reportU = stratifiedSampling(stratum, value, 0.1u"ha", Dict(1 => 15.0u"ha", 2 => 15.0u"ha", 3 => 15.0u"ha"))
    @test resultTable(reportU) == rt
  end

  @test_throws ArgumentError stratifiedSampling(stratum, value, 0.1, Dict(1 => 12.0, 2 => 8.0))
end
