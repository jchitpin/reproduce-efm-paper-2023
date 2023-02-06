## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Compute EFM weights in the two example networks by optimization-based
#     methods across different solvers.
# (2) Solutions are manually aggregated and inputted into the paper tables.
# ------------------------------------------------------------------------------

## USER PARAMETERS -------------------------------------------------------------
# Set working directory containing this script. Example:
cd("/home/jchitpin/Documents/PhD/Projects/reproduce-efm-paper-2023/src/")
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using Tables, CSV
using JuMP
using GLPK, Gurobi, SCIP, COSMO, OSQP, CDDLib, ECOS, ProxSDP, Tulip
include.(filter(contains(r".jl$"), readdir("functions"; join=true)))
# ------------------------------------------------------------------------------

## SETTING UP MATHEMATICAL OPTIMIZATION VARIABLES ------------------------------
# Goal is to create system of linear EFM weight equations for example networks
# 1 and 2 with the following format: A * w == v

# Load EFM matrix from MarkovWeightedEFMs (rows are reactions; cols are EFMs)
A1 = CSV.read("../data/efm-matrix-example-1.csv", Tables.matrix, header=false)
A2 = CSV.read("../data/efm-matrix-example-2.csv", Tables.matrix, header=false)

# Steady state fluxes for example 1 and 2
v1 = [2, 2, 2, 2, 2, 2, 2, 4]
v2 = [3, 3, 2, 3, 2, 2, 2, 3, 2];

# Solvers
milp_solvers, qp_solvers, lp_solvers = load_solvers()
# ------------------------------------------------------------------------------

# QUADRATIC PROGRAMMING --------------------------------------------------------
# Reference: Schwartz and Kanehisa (2006) DOI:10.1186/1471-2105-7-186
efms_1_sk = Vector{Vector{Float64}}()
efms_2_sk = Vector{Vector{Float64}}()
for solver in qp_solvers
  push!(efms_1_sk, jump_qp_l2norm(A1, v1, solver))
  push!(efms_2_sk, jump_qp_l2norm(A2, v2, solver))
end
# ------------------------------------------------------------------------------

# QUADRATIC PROGRAMMING MAXIMIZING ACTIVITY OF SHORTEST PATHWAYS ---------------
# Reference: Orman et al. (2011) DOI:10.1016/j.jtbi.2010.11.042
efms_1_or = Vector{Vector{Float64}}()
efms_2_or = Vector{Vector{Float64}}()
for solver in qp_solvers
  push!(efms_1_or, jump_qp_max_shortest_paths(A1, v1, solver))
  push!(efms_2_or, jump_qp_max_shortest_paths(A2, v2, solver))
end
# ------------------------------------------------------------------------------

# LINEAR PROGRAMMING MAXIMIZING ACTIVITY OF SHORTEST PATHWAYS ------------------
# Reference: Rugen et al. (2012) DOI:10.1016/j.ymben.2012.01.009
efms_1_ru = Vector{Vector{Float64}}()
efms_2_ru = Vector{Vector{Float64}}()
for solver in lp_solvers
  push!(efms_1_ru, jump_lp_max_shortest_paths(A1, v1, solver))
  push!(efms_2_ru, jump_lp_max_shortest_paths(A2, v2, solver))
end
# ------------------------------------------------------------------------------

# LINEAR PROGRAMMING MAXIMIZING LONGEST PATHWAYS -------------------------------
# Reference: Ren et al. (2020) DOI:10.1016/j.algal.2019.101767
efms_1_re = Vector{Vector{Float64}}()
efms_2_re = Vector{Vector{Float64}}()
for solver in lp_solvers
  push!(efms_1_re, jump_lp_max_longest_paths(A1, v1, solver))
  push!(efms_2_re, jump_lp_max_longest_paths(A2, v2, solver))
end
# ------------------------------------------------------------------------------

# MIXED INTEGER PROGRAMMING MAXIMIZING SHORTEST ACTIVE PATH---------------------
# Reference: Nookaew et al. (2007) DOI:10.1002/bit.21339
efms_1_no = Vector{Vector{Float64}}()
efms_2_no = Vector{Vector{Float64}}()
for solver in milp_solvers
  push!(efms_1_no, jump_milp_max_active_paths(A1, v1, solver))
  push!(efms_2_no, jump_milp_max_active_paths(A2, v2, solver))
end
# ------------------------------------------------------------------------------

# EXAMPLE 01 WEIGHTS -----------------------------------------------------------
# Ordering EFMs to match manuscript figure 1a (run `res.e` in main-efm-weights-example-markov.jl)
order_1 = [3, 2, 4, 1]

# Figure 1a example weights
A1_ordered = view(A1, :, order_1) # order EFM matrix columns
A1_ordered * [2, 0, 2, 0] .- v1
A1_ordered * [0, 2, 0, 2] .- v1
A1_ordered * [1.1, 0.9, 1.1, 0.9] .- v1
A1_ordered * [1.5, 0.5, 1.5, 0.5] .- v1

# Figure 2a example weights
efms_1_sk_ordered = efms_1_sk .|> x -> x[order_1]
efms_1_or_ordered = efms_1_or .|> x -> x[order_1]
efms_1_ru_ordered = efms_1_ru .|> x -> x[order_1]
efms_1_re_ordered = efms_1_re .|> x -> x[order_1]
efms_1_no_ordered = efms_1_no .|> x -> x[order_1]

# Table 2: Min L2
hcat(efms_1_sk_ordered...) # all column solutions
qp_solvers                 # all solvers yield same solution

# Table 2: Max qSPA
hcat(efms_1_or_ordered...) # all column solutions
qp_solvers                 # all solvers yield same solution

# Table 2: Max lSPA
hcat(efms_1_ru_ordered...)
lp_solvers[[4]]            # [2,0,2,0] solution solvers
lp_solvers[[2,3,6]]        # [0,2,0,2] solution solvers
lp_solvers[[1,5,7,8]]      # [1,1,1,1] solution solvers

# Table 2: Min lSPA
hcat(efms_1_re_ordered...)
lp_solvers[[4]]            # [2,0,2,0] solution solvers
lp_solvers[[2,3,6]]        # [0,2,0,2] solution solvers
lp_solvers[[1,5,7,8]]      # [1,1,1,1] solution solvers

# Table 2: Min milAP
hcat(efms_1_no_ordered...)
milp_solvers[1]            # [2,0,2,0] solution solver

# EXAMPLE 02 WEIGHTS -----------------------------------------------------------
# Ordering EFMs to match manuscript figure 1a
order_2 = [8, 4, 2, 5, 6, 1, 7, 3]

# Figure 1b example weights
A2_ordered = view(A2, :, order_2) # order EFM matrix columns
A2_ordered * [2, 2, 2, 2, 1, 1, 0, 0] .- v2
A2_ordered * [0, 0, 0, 0, 1, 1, 2, 2] .- v2
A2_ordered * [1, 1, 1, 1, 1, 1, 1, 1] .- v2
A2_ordered * [0.6, 0.6, 0.9, 0.9, 1.3, 0.7, 1.1, 1.4] .- v2
A2_ordered * [0.5, 0.5, 1.0, 1.0, 1.5, 0.5, 1.0, 1.5] .- v2

# Figure 2b example weights
efms_2_sk_ordered = efms_2_sk .|> x -> x[order_2]
efms_2_or_ordered = efms_2_or .|> x -> x[order_2]
efms_2_ru_ordered = efms_2_ru .|> x -> x[order_2]
efms_2_re_ordered = efms_2_re .|> x -> x[order_2]
efms_2_no_ordered = efms_2_no .|> x -> x[order_2]

# Table 3: Min L2
hcat(efms_2_sk_ordered...) # all column solutions
qp_solvers                 # all solvers yield same solution

# Table 3: Max qSPA
hcat(efms_2_or_ordered...) # all column solutions
qp_solvers                 # all solvers yield same solution

# Table 3: Max lSPA
hcat(efms_2_ru_ordered...)
lp_solvers[[1,5,7]]          # [2/3,2/3,2/3,2/3,1.0,1.0,4/3,4/3] solution solvers
lp_solvers[[2,4]]          # [1.0,1.0,0.0,0.0,0.0,2.0,2.0,1.0] solution solvers
lp_solvers[[3]]            # [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0] solution solvers
lp_solvers[[6]]            # [0.0,0.0,0.0,0.0,1.0,1.0,2.0,2.0] solution solvers
lp_solvers[[8]]            # [1.12,1.12,1.12,1.12,1.0,1.0,0.88,0.88] solution solvers

# Table 3: Min lSPA
hcat(efms_2_re_ordered...)
lp_solvers[[1,5,7]]        # [2/3,2/3,2/3,2/3,1.0,1.0,4/3,4/3] solution solvers
lp_solvers[[2]]            # [1.0,1.0,0.0,0.0,0.0,2.0,2.0,1.0] solution solvers
lp_solvers[[3]]            # [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0] solution solvers
lp_solvers[[4]]            # [0.0,0.0,1.0,1.0,2.0,0.0,1.0,2.0] solution solvers
lp_solvers[[6]]            # [0.0,0.0,0.0,0.0,1.0,1.0,2.0,2.0] solution solvers
lp_solvers[[8]]            # [1.12,1.12,1.12,1.12,1.0,1.0,0.88,0.88] solution solvers

# Table 3: Min milAP
hcat(efms_2_no_ordered...)
milp_solvers[1]            # [0,0,0,0,1,1,2,2] solution solver

