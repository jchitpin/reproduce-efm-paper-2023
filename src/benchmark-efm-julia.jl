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
B1 = @benchmark steady_state_efm_distribution(S01, v01)
B2 = @benchmark steady_state_efm_distribution(S02, v02)
B3 = @benchmark steady_state_efm_distribution(S03, v03)
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10
B4 = @benchmark steady_state_efm_distribution(S04, v04)
B5 = @benchmark steady_state_efm_distribution(S05, v05)
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 300 # benchmark stops if time limit exceeded
B6 = @benchmark steady_state_efm_distribution(S06, v06)
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 600 # benchmark stops if time limit exceeded
B7 = @benchmark steady_state_efm_distribution(S07, v07) # 34 GiB memory required
#B8 = @benchmark steady_state_efm_distribution(S08, v08) # infeasible
#B9 = @benchmark steady_state_efm_distribution(S09, v09) # infeasible
#B10 = @benchmark steady_state_efm_distribution(S10, v10) # infeasible
#B11 = @benchmark steady_state_efm_distribution(S11, v11) # infeasible
# ------------------------------------------------------------------------------

## CONFIRMING JULIA CODE ENUMERATES SAME NUMBER OF EFMS ------------------------
R01 = steady_state_efm_distribution(S01, v01)
R02 = steady_state_efm_distribution(S02, v02)
R03 = steady_state_efm_distribution(S03, v03)
R04 = steady_state_efm_distribution(S04, v04)
R05 = steady_state_efm_distribution(S05, v05)
R06 = steady_state_efm_distribution(S06, v06)
R07 = steady_state_efm_distribution(S07, v07) # 34 GiB memory required
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

## EXPORTING COMBINED MATLAB AND JULIA DATA TO TEXT ----------------------------
mat = CSV.read(#
  "../data/benchmark-stoich-fluxmodecalculator.csv",
  Tables.matrix,
  header=true
)
m(x) = sum(x)/length(x) * 1e-9 # mean time in seconds
B = [B1.times, B2.times, B3.times, B4.times, B5.times, B6.times, B7.times, [0], [0], [0], [0]]
EFM = [#
  length(R01.e),
  length(R02.e),
  length(R03.e),
  length(R04.e),
  length(R05.e),
  length(R06.e),
  length(R07.e),
   28695, # hardcoded
   34883, # hardcoded
   87207, # hardcoded
   108868 # hardcoded
]
res1 = hcat(EFM, mat[:,2])
res2 = hcat(EFM, m.(B))

# Simple linear regression of time on number of EFMs
linreg(x, y) = hcat(fill!(similar(x), 1), x) \ y
linreg(res1[:,1], res1[:,2]) # slope and intercept of FluxModeCalculator
linreg(res2[:,1], res2[:,2]) # slope and intercept of MarkovWeightedEFMs.jl

CSV.write(#
  "../data/scatterplot-benchmark-fluxmodecalculator.csv",
  Tables.table(res1),
  header = ["efms", "matlab"],
  delim = '\t'
)
CSV.write(#
  "../data/scatterplot-benchmark-markovweightedefms.jl.csv",
  Tables.table(res2[1:6,:]),
  header = ["efms", "julia"],
  delim = '\t'
)
# ------------------------------------------------------------------------------


