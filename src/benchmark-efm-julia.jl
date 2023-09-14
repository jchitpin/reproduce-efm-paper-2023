## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Load stoichiometry matrices generated/benchmarked in MATLAB with EFMs
#     computed using FluxModeCalculator.
# (2) Benchmark the same stoichiometry matrices with the Julia code.
# (3) Export the data for plotting.
# ------------------------------------------------------------------------------

## USER PARAMETERS -------------------------------------------------------------
# Set working directory containing this script. Example:
cd("/home/jchitpin/Documents/PhD/Projects/reproduce-efm-paper-2023/src/")
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using MarkovWeightedEFMs
using Tables, CSV, JuMP
using BenchmarkTools
using GLPK, Gurobi, SCIP, COSMO, OSQP, CDDLib, ECOS, ProxSDP, Tulip
include.(filter(contains(r".jl$"), readdir("functions"; join=true)))
# ------------------------------------------------------------------------------

## LOAD SIMULATED STOICHIOMETRY/EFM MATRICES AND FLUXES ------------------------
S01 = CSV.read("../data/stoich-20-30.csv", Tables.matrix, header=false)
S02 = CSV.read("../data/stoich-20-35.csv", Tables.matrix, header=false)
S03 = CSV.read("../data/stoich-20-40.csv", Tables.matrix, header=false)
S04 = CSV.read("../data/stoich-20-45.csv", Tables.matrix, header=false)
S05 = CSV.read("../data/stoich-20-50.csv", Tables.matrix, header=false)
S06 = CSV.read("../data/stoich-20-55.csv", Tables.matrix, header=false)
S07 = CSV.read("../data/stoich-20-60.csv", Tables.matrix, header=false)
S08 = CSV.read("../data/stoich-20-65.csv", Tables.matrix, header=false)
S09 = CSV.read("../data/stoich-20-70.csv", Tables.matrix, header=false)
S10 = CSV.read("../data/stoich-20-75.csv", Tables.matrix, header=false)
S11 = CSV.read("../data/stoich-20-80.csv", Tables.matrix, header=false)
v01 = vec(CSV.read("../data/ss-flux-20-30.csv", Tables.matrix, header=false))
v02 = vec(CSV.read("../data/ss-flux-20-35.csv", Tables.matrix, header=false))
v03 = vec(CSV.read("../data/ss-flux-20-40.csv", Tables.matrix, header=false))
v04 = vec(CSV.read("../data/ss-flux-20-45.csv", Tables.matrix, header=false))
v05 = vec(CSV.read("../data/ss-flux-20-50.csv", Tables.matrix, header=false))
v06 = vec(CSV.read("../data/ss-flux-20-55.csv", Tables.matrix, header=false))
v07 = vec(CSV.read("../data/ss-flux-20-60.csv", Tables.matrix, header=false))
v08 = vec(CSV.read("../data/ss-flux-20-65.csv", Tables.matrix, header=false))
v09 = vec(CSV.read("../data/ss-flux-20-70.csv", Tables.matrix, header=false))
v10 = vec(CSV.read("../data/ss-flux-20-75.csv", Tables.matrix, header=false))
v11 = vec(CSV.read("../data/ss-flux-20-80.csv", Tables.matrix, header=false))
E01 = CSV.read("../data/stoich-20-30-efm-matrix.csv", Tables.matrix, header=false)
E02 = CSV.read("../data/stoich-20-35-efm-matrix.csv", Tables.matrix, header=false)
E03 = CSV.read("../data/stoich-20-40-efm-matrix.csv", Tables.matrix, header=false)
E04 = CSV.read("../data/stoich-20-45-efm-matrix.csv", Tables.matrix, header=false)
E05 = CSV.read("../data/stoich-20-50-efm-matrix.csv", Tables.matrix, header=false)
E06 = CSV.read("../data/stoich-20-55-efm-matrix.csv", Tables.matrix, header=false)
E07 = CSV.read("../data/stoich-20-60-efm-matrix.csv", Tables.matrix, header=false)
E08 = CSV.read("../data/stoich-20-65-efm-matrix.csv", Tables.matrix, header=false)
E09 = CSV.read("../data/stoich-20-70-efm-matrix.csv", Tables.matrix, header=false)
E10 = CSV.read("../data/stoich-20-75-efm-matrix.csv", Tables.matrix, header=false)
E11 = CSV.read("../data/stoich-20-80-efm-matrix.csv", Tables.matrix, header=false)
# ------------------------------------------------------------------------------

## BENCHMARKING JULIA CODE TO ENUMERATE EFMS AND ASSIGN WEIGHTS ----------------
# May need to adjust benchmarking time limits for your CPU hardware
BenchmarkTools.DEFAULT_PARAMETERS.samples = 100 # 100 samples
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5
B1 = @benchmark steady_state_efm_distribution(S01, v01; solver=:Arnoldi)
B2 = @benchmark steady_state_efm_distribution(S02, v02; solver=:Arnoldi)
B3 = @benchmark steady_state_efm_distribution(S03, v03; solver=:Arnoldi)
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 15
B4 = @benchmark steady_state_efm_distribution(S04, v04; solver=:Arnoldi)
B5 = @benchmark steady_state_efm_distribution(S05, v05; solver=:Arnoldi)
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 300 # benchmark stops if time limit exceeded
B6 = @benchmark steady_state_efm_distribution(S06, v06; solver=:Arnoldi)
B7 = @benchmark steady_state_efm_distribution(S07, v07; solver=:Arnoldi, issparse=true) # 17.109 seconds (11.95 GiB)
#B8 = @benchmark steady_state_efm_distribution(S08, v08; solver=:Arnoldi, issparse=true) # took too long
#B8 = @benchmark steady_state_efm_distribution(S08, v08; solver=:Arnoldi, issparse=true) # took too long
#B9 = @benchmark steady_state_efm_distribution(S09, v09; solver=:Arnoldi, issparse=true) # took too long
#B10 = @benchmark steady_state_efm_distribution(S10, v10; solver=:Arnoldi, issparse=true) # took too long
#B11 = @benchmark steady_state_efm_distribution(S11, v11; solver=:Arnoldi, issparse=true) # took too long
# ------------------------------------------------------------------------------

## CONFIRMING JULIA CODE ENUMERATES SAME NUMBER OF EFMS ------------------------
R01 = steady_state_efm_distribution(S01, v01)
R02 = steady_state_efm_distribution(S02, v02)
R03 = steady_state_efm_distribution(S03, v03)
R04 = steady_state_efm_distribution(S04, v04)
R05 = steady_state_efm_distribution(S05, v05)
R06 = steady_state_efm_distribution(S06, v06)
R07 = steady_state_efm_distribution(S07, v07)
#R08 = steady_state_efm_distribution(S08, v08) # infeasible
#R09 = steady_state_efm_distribution(S09, v09) # infeasible
#R10 = steady_state_efm_distribution(S10, v10) # infeasible
#R11 = steady_state_efm_distribution(S11, v11) # infeasible

# E0X matrix was computed from FluxModeCalculator
idx_intersect_01 = compare_efms(E01, R01.e, S01); # 20 EFMs
idx_intersect_02 = compare_efms(E02, R02.e, S02); # 41 EFMs
idx_intersect_03 = compare_efms(E03, R03.e, S03); # 182 EFMs
idx_intersect_04 = compare_efms(E04, R04.e, S04); # 438 EFMs
idx_intersect_05 = compare_efms(E05, R05.e, S05); # 289 EFMs
idx_intersect_06 = compare_efms(E06, R06.e, S06); # 2478 EFMs
idx_intersect_07 = compare_efms(E07, R07.e, S07); # 6860 EFMs
#idx_intersect_08 = compare_efms(E08, R08.e, S08); # 28695 EFMs
#idx_intersect_09 = compare_efms(E09, R09.e, S09); # 34883 EFMs
#idx_intersect_10 = compare_efms(E10, R10.e, S10); # 87207 EFMs
#idx_intersect_11 = compare_efms(E11, R11.e, S11); # 108868 EFMs
# ------------------------------------------------------------------------------

## SOLVING EFM WEIGHTS BY LP OPTIMIZATION-BASED METHOD -------------------------
BenchmarkTools.DEFAULT_PARAMETERS.samples = 100 # 100 samples
_, _, lp_solvers = load_solvers()
C1 = @benchmark jump_lp_max_shortest_paths(E01, v01, lp_solvers[1]) #  μs median.  memory
C2 = @benchmark jump_lp_max_shortest_paths(E02, v02, lp_solvers[1]) #  μs median.  memory
C3 = @benchmark jump_lp_max_shortest_paths(E03, v03, lp_solvers[1]) #  μs median.  memory
C4 = @benchmark jump_lp_max_shortest_paths(E04, v04, lp_solvers[1]) #  μs median.  memory
C5 = @benchmark jump_lp_max_shortest_paths(E05, v05, lp_solvers[1]) #  μs median.  memory
C6 = @benchmark jump_lp_max_shortest_paths(E06, v06, lp_solvers[1]) #  μs median.  memory
C7 = @benchmark jump_lp_max_shortest_paths(E07, v07, lp_solvers[1]) #  μs median.  memory
C8 = @benchmark jump_lp_max_shortest_paths(E08, v08, lp_solvers[1]) #  μs median.  memory
C9 = @benchmark jump_lp_max_shortest_paths(E09, v09, lp_solvers[1]) #  μs median.  memory
C10 = @benchmark jump_lp_max_shortest_paths(E10, v10, lp_solvers[1]) #  μs median.  memory
C11 = @benchmark jump_lp_max_shortest_paths(E11, v11, lp_solvers[1]) #  μs median.  memory
# ------------------------------------------------------------------------------




## EXPORTING COMBINED MATLAB AND JULIA DATA TO TEXT ----------------------------
mat = CSV.read(#
  "../data/benchmark-stoich-fluxmodecalculator.csv",
  Tables.matrix,
  header=true
)
m(x) = sum(x)/length(x) * 1e-9 # mean time in seconds
B = [B1.times, B2.times, B3.times, B4.times, B5.times, B6.times, B7.times, [0], [0], [0], [0]]
C = [C1.times, C2.times, C3.times, C4.times, C5.times, C6.times, C7.times, C8.times, C9.times, C10.times, C11.times]
EFM = [#
  length(R01.e),
  length(R02.e),
  length(R03.e),
  length(R04.e),
  length(R05.e),
  length(R06.e),
  6860,  # hardcoded
  28695, # hardcoded
  34883, # hardcoded
  87207, # hardcoded
  108868 # hardcoded
]
res1 = hcat(EFM, mat[:,2] .+ m.(C))
res2 = hcat(EFM, m.(B))

res1_log = deepcopy(res1)
res1_log[:,1] = log10.(res1_log[:,1])
res1_log[:,2] = log10.(res1_log[:,2])

res2_log = deepcopy(res2[1:7,:])
res2_log[:,1] = log10.(res2_log[:,1])
res2_log[:,2] = log10.(res2_log[:,2])

# Note: Simple linear regression of run time done in PGFPlots
CSV.write(#
  "../data/scatterplot-benchmark-fluxmodecalculator.csv",
  Tables.table(res1),
  header = ["efms", "matlab"],
  delim = '\t'
)
CSV.write(#
  "../data/scatterplot-benchmark-markovweightedefms.jl.csv",
  Tables.table(res2[1:7,:]),
  header = ["efms", "julia"],
  delim = '\t'
)
CSV.write(#
  "../data/scatterplot-benchmark-fluxmodecalculator-log10.csv",
  Tables.table(res1_log),
  header = ["efms", "matlab"],
  delim = '\t'
)
CSV.write(#
  "../data/scatterplot-benchmark-markovweightedefms.jl-log10.csv",
  Tables.table(res2_log),
  header = ["efms", "julia"],
  delim = '\t'
)
# ------------------------------------------------------------------------------




