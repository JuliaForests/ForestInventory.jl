using ForestInventory
using Test
using DataFrames
using ForestFoundations
using Statistics
using Unitful

@testset "ForestInventory.jl" begin
  include("simplecasualsamplingtests.jl")
  include("stratifiedsamplingtests.jl")
  include("systematicsamplingtests.jl")
  include("clustersamplingtests.jl")
  include("horizontalpointsamplingtests.jl")
  include("multistartsystematicsamplingtests.jl")
  include("twostagesamplingtests.jl")
  include("independentoccasionssamplingtests.jl")
  include("completereplacementsamplingtests.jl")
  include("partialreplacementsamplingtests.jl")
  include("doublesamplingtests.jl")
end
