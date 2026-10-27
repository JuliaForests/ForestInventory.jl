using ForestInventory
using Test
using DataFrames
using ForestFoundations
using Statistics
using Unitful
using Distributions

include("bookdata.jl")

# Student's t for a confidence level α and df degrees of freedom
tq(df, α=0.95) = -quantile(TDist(df), (1 - α) / 2)

@testset "ForestInventory.jl" begin
  include("simplecasualsamplingtests.jl")
  include("stratifiedsamplingtests.jl")
  include("systematicsamplingtests.jl")
  include("clustersamplingtests.jl")
  include("horizontalpointsamplingtests.jl")
  include("prodansamplingtests.jl")
  include("strandsamplingtests.jl")
  include("threepsamplingtests.jl")
  include("multistartsystematicsamplingtests.jl")
  include("twostagesamplingtests.jl")
  include("independentoccasionssamplingtests.jl")
  include("completereplacementsamplingtests.jl")
  include("partialreplacementsamplingtests.jl")
  include("doublesamplingtests.jl")
end
