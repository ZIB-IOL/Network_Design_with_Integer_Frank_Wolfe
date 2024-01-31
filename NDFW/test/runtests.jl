using NDFW
using Boscia
using Test
using FrankWolfe
using Random
using SCIP
import MathOptInterface
const MOI = MathOptInterface

#include("dataset_tests.jl")
include("extreme_point_cut_test.jl")

# for file in readdir(joinpath(@__DIR__, "../examples/"), join=true)
#     if endswith(file, "jl")
#         include(file)
#     end
# end