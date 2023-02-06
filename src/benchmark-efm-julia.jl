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
# Set benchmarking parameters for 100 samples
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 12_000 # benchmark stops if time limit exceeded
BenchmarkTools.DEFAULT_PARAMETERS.samples = 100
B1 = @benchmark steady_state_efm_distribution(S01, v01)
B2 = @benchmark steady_state_efm_distribution(S02, v02)
B3 = @benchmark steady_state_efm_distribution(S03, v03)
B4 = @benchmark steady_state_efm_distribution(S04, v04)
B5 = @benchmark steady_state_efm_distribution(S05, v05)
B6 = @benchmark steady_state_efm_distribution(S06, v06)
#B7 = @benchmark steady_state_efm_distribution(S07, v07) # infeasible
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
#R07 = steady_state_efm_distribution(S07, v07) # infeasible
#R08 = steady_state_efm_distribution(S08, v08) # infeasible
#R09 = steady_state_efm_distribution(S09, v09) # infeasible
#R10 = steady_state_efm_distribution(S10, v10) # infeasible
#R11 = steady_state_efm_distribution(S11, v11) # infeasible

# It so happens that FluxModeCalculator returns duplicate EFMs while
# the Julia code returns the exact number of EFMs. However, the Julia code
# fails to scale beyond approximately 2000 EFMs. The memory requirement to
# generate the cycle-history Markov chain increases exponentially with EFM
# count and the code seems to fail when computing the steady state
# probabilities using the QuantEcon package
idx_intersect_01 = compare_efms(E01, R01.e, S01); # 20 EFMs
idx_intersect_02 = compare_efms(E02, R02.e, S02); # 41 EFMs
idx_intersect_03 = compare_efms(E03, R03.e, S03); # 175 not 182 EFMs
idx_intersect_04 = compare_efms(E04, R04.e, S04); # 423 not 438 EFMs
idx_intersect_05 = compare_efms(E05, R05.e, S05); # 279 not 289 EFMs
idx_intersect_06 = compare_efms(E06, R06.e, S06); # 1900 not 2478 EFMs
idx_intersect_07 = compare_efms(E07, R06.e, S07); # 3650 not 6860 EFMs
idx_intersect_08 = compare_efms(E08, R06.e, S08); # 13122 not 28695 EFMs
idx_intersect_09 = compare_efms(E09, R06.e, S09); # 11538 not 34883 EFMs
idx_intersect_10 = compare_efms(E10, R06.e, S10); # 28178 not 87207 EFMs
idx_intersect_11 = compare_efms(E11, R06.e, S11); # 27897 not 108868 EFMs
# ------------------------------------------------------------------------------

## EXPORTING COMBINED MATLAB AND JULIA DATA TO TEXT ----------------------------
mat = CSV.read(#
  "../data/benchmark-stoich-fluxmodecalculator.csv",
  Tables.matrix,
  header=true
)
m(x) = sum(x)/length(x) * 10e-10 # mean time in seconds
n(x) = length(x)
B = [B2.times, B1.times, B3.times, B4.times, B5.times, B6.times, [0], [0], [0], [0], [0]]
EFM = [#
   n(idx_intersect_01),
   n(idx_intersect_02),
   n(idx_intersect_03),
   n(idx_intersect_04),
   n(idx_intersect_05),
   n(idx_intersect_06),
   3650, # hardcoded
   13122, # hardcoded
   11538, # hardcoded
   28178, # hardcoded
   27897 # hardcoded
]
res1 = hcat(EFM, mat[:,2])
res2 = hcat(EFM, m.(B))

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


