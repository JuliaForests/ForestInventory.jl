@testset "threepsampling" begin
  predicted = [0.20, 0.45, 0.30, 0.55, 0.15, 0.60, 0.25, 0.40, 0.35, 0.50]u"m^3"
  measured = Vector{Union{Missing,typeof(1.0u"m^3")}}(missing, 10)
  measured[2] = 0.48u"m^3"
  measured[4] = 0.58u"m^3"
  measured[6] = 0.63u"m^3"
  measured[10] = 0.52u"m^3"

  @test threepsampling(predicted, measured, 5u"ha") isa SamplingReport

  rt = threepsampling(predicted, measured, 5u"ha").resulttable
  @test rt isa DataFrame
  @test rt.n[1] == 4
  @test rt.N[1] == 10
  @test eltype(rt.X) <: Unitful.Quantity

  @testset "matches Grosenbaugh's ratio-of-means estimator by hand" begin
    idx = findall(!ismissing, measured)
    K = sum(ustrip.(predicted))
    R = [ustrip(measured[i]) / ustrip(predicted[i]) for i in idx]
    n, N = length(idx), length(predicted)
    R̄ = mean(R)
    meanvar = (var(R) / n) * (1 - n / N)
    total = K * R̄
    se = K * sqrt(meanvar)
    t = -quantile(TDist(n - 1), 0.025)

    @test ustrip(rt.X[1]) ≈ total atol = 1e-9
    @test rt.Rm[1] ≈ R̄ atol = 1e-9
    @test ustrip(rt.K[1]) ≈ K atol = 1e-9
    @test ustrip(rt.sxm[1]) ≈ se atol = 1e-9
    @test ustrip(rt.Ea[1]) ≈ t * se atol = 1e-9
    @test rt.xha[1] ≈ rt.X[1] / 5u"ha"
  end

  @testset "order of rows does not matter" begin
    perm = [10, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    rtPerm = threepsampling(predicted[perm], measured[perm], 5u"ha").resulttable
    @test rtPerm == rt
  end

  @testset "bare (unitless) numbers behave the same as Unitful volumes" begin
    predictedBare = ustrip.(predicted)
    measuredBare = Vector{Union{Missing,Float64}}(missing, 10)
    measuredBare[[2, 4, 6, 10]] = ustrip.(skipmissing(measured))
    rtBare = threepsampling(predictedBare, measuredBare, 5u"ha").resulttable
    @test rtBare.X[1] ≈ ustrip(rt.X[1]) atol = 1e-9
  end

  @testset "near-zero-variance ratios no longer crash _requiredsamplesize" begin
    # every measured tree came in almost exactly at its prediction -- the outcome the
    # method is designed to produce -- which used to solve the required-size iteration
    # below 1 observation and error out of TDist.
    predictedPerfect = fill(1.0, 10) .* u"m^3"
    measuredPerfect = Vector{Union{Missing,typeof(1.0u"m^3")}}(missing, 10)
    measuredPerfect[[1, 3, 5, 7]] .= 1.0u"m^3"
    rtPerfect = threepsampling(predictedPerfect, measuredPerfect, 5u"ha").resulttable
    @test rtPerfect.nreq[1] >= 2
  end

  @test_throws DimensionMismatch threepsampling(predicted, measured[1:9], 5u"ha")

  @testset "fewer than two measured trees" begin
    measuredOne = Vector{Union{Missing,typeof(1.0u"m^3")}}(missing, 10)
    measuredOne[1] = 0.2u"m^3"
    @test_throws ArgumentError threepsampling(predicted, measuredOne, 5u"ha")
  end
end
