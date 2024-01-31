module NDFW

using Boscia
using FrankWolfe
using Random
using LinearAlgebra
using SCIP
using Graphs
#using Plots
using GraphRecipes
using GraphIO
using OrderedCollections
using DataFrames
using BinDeps
using DelimitedFiles
using SparseArrays
using StatsBase
using Printf
using Optim
using Gurobi
using JuMP
using SCS
using HiGHS
using Hypatia
using Pajarito
using Bonobo

import MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

# # assorted utility functions
# include("scip/util.jl")
# # load deps, version check
# include("scip/init.jl")
# # wrapper of SCIP library
# include("scip/wrapper.jl")
# # memory management
# include("scip/scip_data.jl")
# # separators
# include("scip/sepa.jl")
# # cut selectors
# include("scip/cut_selector.jl")
# # constraint handlers
# include("scip/conshdlr.jl")
# # constraints from nonlinear expressions
# include("scip/nonlinear.jl")
# # implementation of MOI
# include("scip/MOI_wrapper.jl")
# # convenience functions
# include("scip/convenience.jl")
# # warn about rewrite
# include("scip/compat.jl")

include("network_design_problem/data_structures.jl")
include("utilities.jl")
include("network_design_problem/net_design.jl")
include("shortest_path_methods.jl")
include("lmo_utilities.jl")
include("moi_network_lmo.jl")

#Traffic Assignment Library
include("ta_library/load_network.jl")
include("ta_library/ta_data_processing.jl")
include("ta_library/frank_wolfe.jl")
include("ta_library/modified_methods.jl")


#Methods
include("network_design_problem/methods/direct_solve/direct_scip_solve.jl")
include("network_design_problem/methods/integer_frank_wolfe_method/traffic_flow_polytope.jl")
include("network_design_problem/methods/integer_frank_wolfe_method/network_design_constraints.jl")
include("network_design_problem/methods/integer_frank_wolfe_method/ifw_method.jl")
include("network_design_problem/methods/integer_frank_wolfe_method/ifw_func_and_grad.jl")
include("network_design_problem/methods/direct_solve/conic_method/conic_formulation.jl")
include("network_design_problem/methods/direct_solve/conic_method/conic_method.jl")
include("network_design_problem/methods/bounded_lmo_method/bndlmo.jl")
include("network_design_problem/methods/bounded_lmo_method/bounded_lmo_penalty.jl")
include("network_design_problem/methods/bounded_lmo_method/blmo_func_and_grad.jl")

include("network_design_problem/methods/benders_decomposition_method/benders_decomposition_method.jl")
include("network_design_problem/methods/benders_decomposition_method/decomposition.jl")
include("network_design_problem/methods/benders_decomposition_method/bd_lmo.jl")
include("network_design_problem/methods/benders_decomposition_method/bd_utilities.jl")
include("network_design_problem/methods/benders_decomposition_method/bd_cuts.jl")

include("traffic_assignment_problem/methods/lp_fw_method/lp_ta_method.jl")
include("traffic_assignment_problem/methods/lp_fw_method/lp_func_and_grad.jl")

include("traffic_assignment_problem/methods/sp_fw_method/sp_ta_method.jl")
include("traffic_assignment_problem/methods/sp_fw_method/splmo.jl")
include("traffic_assignment_problem/methods/sp_fw_method/sp_func_and_grad.jl")

include("traffic_assignment_problem/ta_lmo_utilities.jl")

#The following implementations are not working with the current version of the code

#include("network_design_problem/methods/older_implementations/scip_conversions/NetworkOptimizer_v2.jl")
#include("network_design_problem/methods/older_implementations/scip_conversions/network_constants.jl")
#include("network_design_problem/methods/older_implementations/benders_decomposition_method/decompostion.jl")
#include("network_design_problem/methods/older_implementations/benders_decomposition_method/bd_lmo.jl")

#include("network_design_problem/methods/older_implementations/scip_conversions/network_variable.jl")
#include("network_design_problem/methods/older_implementations/scip_conversions/network_objective.jl")
#include("network_design_problem/methods/older_implementations/scip_conversions/network_results.jl")
#include("network_design_problem/methods/older_implementations/scip_conversions/network_constraints.jl")
#include("network_design_problem/methods/older_implementations/penalty_method/penalty_lmo.jl")
#include("network_design_problem/methods/older_implementations/penalty_method/penalty_method.jl")
#include("network_design_problem/methods/older_implementations/penalty_method/network_lmop.jl")
#include("network_design_problem/methods/older_implementations/unused_utilities.jl")


end # module
