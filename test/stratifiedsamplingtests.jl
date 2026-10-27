@testset "stratifiedsampling" begin
  # The example of the book (chapter 6): 3 strata of the population of 450 plots of 0.1 ha, with
  # 144, 164 and 142 plots (Wh = 0.320, 0.364, 0.316), and the 22 plots of Quadro 6.2. The book
  # takes the effective degrees of freedom of Satterthwaite, 16.21, up to 17 (t = 2.11).
  stratum = str_stratum
  value = str_value * u"m^3"
  strataarea = Dict(1 => 14.4u"ha", 2 => 16.4u"ha", 3 => 14.2u"ha")

  report = stratifiedsampling(stratum, value, 0.1u"ha", strataarea)
  @test report isa SamplingReport
  rt = resulttable(report)
  @test nrow(rt) == 1

  @test ustrip(rt.xm[1]) ≈ 21.465 atol = 1e-3
  @test ustrip(rt.s2xm[1]) ≈ 0.8391 atol = 1e-3
  @test ustrip(rt.sxm[1]) ≈ 0.9160 atol = 1e-3
  @test ustrip(rt.Ea[1]) ≈ tq(17) * ustrip(rt.sxm[1])
  @test ustrip(rt.Ea[1]) ≈ 1.9328 atol = 2e-3
  @test rt.Er[1] ≈ 9.0 atol = 5e-2
  @test ustrip(rt.X[1]) ≈ 9659 atol = 1
  @test ustrip(rt.Xl[1]) ≈ 8789 atol = 2
  @test ustrip(rt.Xu[1]) ≈ 10529 atol = 2
  @test rt.cv[1] ≈ sqrt(19.452) / 21.465 * 100 atol = 1e-2
  @test rt.population[1] == "finite"
  @test rt.nh[1] == (7, 8, 7)
  @test rt.n[1] == 22
  @test rt.N[1] == 450
  # the book rounds 18.05 plots to 18, the package up to 19
  @test rt.nreq[1] in 18:19
  @test rt.nhreq[1] == (7, 7, 6)
  @test eltype(rt.xm) <: Unitful.Quantity

  aux = auxiliarytable(report)
  @test nrow(aux) == 3
  @test aux.nh == [7, 8, 7]
  @test ustrip(aux.xmh) ≈ [12.029, 24.563, 27.457] atol = 1e-3
  @test ustrip(aux.s2h) ≈ [10.462, 17.483, 30.823] atol = 1e-3
  @test aux.Wh ≈ [0.320, 0.364, 0.316] atol = 1e-3
  @test sum(ustrip(aux.Whs2h)) ≈ 19.452 atol = 1e-2    # the variance of the stratified population

  # the analysis of variance of the book (Quadro 6.1): SQe 951.355, SQd 370.090, F 24.42
  av = anova(report)
  @test nrow(av) == 3
  @test av.gl == [2, 19, 21]
  @test ustrip(av.SQ[2]) ≈ 370.090 atol = 1e-2
  @test ustrip(av.SQ[3]) ≈ 1321.49 atol = 5e-2
  @test av.F[1] ≈ 24.42 atol = 1e-2

  @testset "Neyman's optimum allocation" begin
    ny = resulttable(stratifiedsampling(stratum, value, 0.1u"ha", strataarea; allocation=:neyman))
    Wh, sh = aux.Wh, ustrip(aux.sh)
    E = 0.1 * 21.465
    n1 = tq(17)^2 * sum(Wh .* sh)^2 / (E^2 + tq(17)^2 * sum(Wh .* sh .^ 2) / 450)    # Eq. 6.34
    @test abs(ny.nreq[1] - n1) < 1.5
    @test ny.nhreq[1] == ceil.(Int, ny.nreq[1] * Wh .* sh / sum(Wh .* sh)) |> Tuple
    @test_throws ArgumentError stratifiedsampling(stratum, value, 0.1u"ha", strataarea; allocation=:optimum)
  end

  @testset "strata given out of sorted order still match the right area" begin
    # rows collected in the field are rarely already sorted by stratum -- the dict-based
    # strataarea must resolve by key, not by position, so shuffling stratum/value
    # (together) must not change the result at all.
    shuffled = [15, 3, 20, 8, 1, 22, 11, 5, 17, 9, 2, 19, 6, 14, 21, 10, 4, 18, 7, 13, 16, 12]  # a permutation of 1:22
    rtShuffled = resulttable(stratifiedsampling(stratum[shuffled], value[shuffled], 0.1u"ha", strataarea))
    # summing the same floats in a different order can differ in the last bit, so compare
    # numerically rather than requiring exact equality
    @test ustrip(rtShuffled.xm[1]) ≈ ustrip(rt.xm[1])
    @test ustrip(rtShuffled.X[1]) ≈ ustrip(rt.X[1])
    @test rtShuffled.nh[1] == rt.nh[1]
    @test rtShuffled.N[1] == rt.N[1]
  end

  @testset "one stratum matches simplecasualsampling" begin
    v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1]u"m^3"
    stratum1 = fill(1, 11)
    r1 = resulttable(stratifiedsampling(stratum1, v, 0.05u"ha", Dict(1 => 10.0u"ha")))
    r2 = simplecasualsampling(v, 0.05u"ha", 10u"ha")
    @test ustrip(r1.xm[1]) ≈ ustrip(r2.xm[1]) atol = 1e-6
    @test ustrip(r1.s2xm[1]) ≈ ustrip(r2.s2xm[1]) atol = 1e-6
    @test ustrip(r1.sxm[1]) ≈ ustrip(r2.sxm[1]) atol = 1e-6
    @test r1.nreq[1] == r2.nreq[1]
    @test r1.population[1] == r2.population[1]
    @test r1.f[1] ≈ r2.f[1] atol = 1e-6
    # the ANOVA between-strata line is undefined with a single stratum (0 numerator DOF)
    r1full = stratifiedsampling(stratum1, v, 0.05u"ha", Dict(1 => 10.0u"ha"))
    @test ismissing(anova(r1full).F[1])
  end

  @testset "units" begin
    reportU = stratifiedsampling(stratum, value, 0.1u"ha", strataarea)
    @test resulttable(reportU) == rt
  end

  @test_throws ArgumentError stratifiedsampling(stratum, value, 0.1u"ha", Dict(1 => 14.4u"ha", 2 => 16.4u"ha"))
end
