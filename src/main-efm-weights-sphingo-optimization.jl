## DESCRIPTION -----------------------------------------------------------------
# This script does the following (for wildtype/disease sphingolipid dataset):
# (1)  Compute EFM weights using LP, MILP, and QP formulations described
#      in literature using various solvers.
# (2a) Compute the log10(squared reconstruction error) for the Markov solution
# (2b) Compute the log10(squared reconstruction error) for optimization methods
# (3)  Export log10(squared reconstruction error) for optimization methods
# (4)  Export EFM weights from each method/solver.
# ------------------------------------------------------------------------------
#
## USER PARAMETERS -------------------------------------------------------------
# Set working directory of this script. Example:
cd("/home/jchitpin/Documents/PhD/Projects/reproduce-efm-paper-2023/src/")
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using Tables, CSV, Plots
using JuMP, BenchmarkTools
using GLPK, Gurobi, SCIP, COSMO, OSQP, CDDLib, ECOS, ProxSDP, Tulip
include.(filter(contains(r".jl$"), readdir("functions"; join=true)))
# ------------------------------------------------------------------------------

## LOAD RELEVANT DATA ----------------------------------------------------------
fluxes_wt = vec(CSV.read("../data/fluxes-wt.csv", Tables.matrix, header=false))
fluxes_ad = vec(CSV.read("../data/fluxes-ad.csv", Tables.matrix, header=false))
efms = CSV.read("../data/efm-matrix-sphingo.csv", Tables.matrix, header=false)
w_raw_wt = vec(CSV.read("../data/efm-weights-wt-raw-markov.csv", Tables.matrix, header=false))
w_raw_ad = vec(CSV.read("../data/efm-weights-ad-raw-markov.csv", Tables.matrix, header=false))
w_log_wt = vec(CSV.read("../data/efm-weights-wt-log-markov.csv", Tables.matrix, header=false))
w_log_ad = vec(CSV.read("../data/efm-weights-ad-log-markov.csv", Tables.matrix, header=false))
p_raw_wt = vec(CSV.read("../data/efm-proport-wt-raw-markov.csv", Tables.matrix, header=false))
p_raw_ad = vec(CSV.read("../data/efm-proport-ad-raw-markov.csv", Tables.matrix, header=false))
p_log_wt = vec(CSV.read("../data/efm-proport-wt-log-markov.csv", Tables.matrix, header=false))
p_log_ad = vec(CSV.read("../data/efm-proport-ad-log-markov.csv", Tables.matrix, header=false))
# ------------------------------------------------------------------------------

## SETTING UP MATHEMATICAL OPTIMIZATION VARIABLES ------------------------------
# Renaming variables to follow convention: Ax=v
A = copy(efms)
v1 = fluxes_wt
v2 = fluxes_ad

# Solvers
milp_solvers, qp_solvers, lp_solvers = load_solvers()
# ------------------------------------------------------------------------------

# QUADRATIC PROGRAMMING --------------------------------------------------------
# Reference: Schwartz and Kanehisa (2006) DOI:10.1186/1471-2105-7-186
sk_wt = optimization_sk(A, v1, qp_solvers)
sk_ad = optimization_sk(A, v2, qp_solvers)

b_sk_wt = @benchmark jump_qp_l2norm(A, v1, qp_solvers[1]) # 1.823 ms median. 1.41 MiB memory
b_sk_ad = @benchmark jump_qp_l2norm(A, v2, qp_solvers[1]) # 2.086 ms median. 1.58 MiB
# ------------------------------------------------------------------------------

# QUADRATIC PROGRAMMING MAXIMIZING ACTIVITY OF SHORTEST PATHWAYS ---------------
# Reference: Orman et al. (2011) DOI:10.1016/j.jtbi.2010.11.042
or_wt = optimization_or(A, v1, qp_solvers)
or_ad = optimization_or(A, v2, qp_solvers)
b_or_wt = @benchmark jump_qp_max_shortest_paths(A, v1, qp_solvers[1]) # 1.948 ms median. 1.62 MiB memory
b_or_ad = @benchmark jump_qp_max_shortest_paths(A, v2, qp_solvers[1]) # 1.998 ms median. 1.60 MiB memory
# ------------------------------------------------------------------------------

# LINEAR PROGRAMMING MAXIMIZING ACTIVITY OF SHORTEST PATHWAYS ------------------
# Reference: Rugen et al. (2012) DOI:10.1016/j.ymben.2012.01.009
ru_wt = optimization_ru(A, v1, lp_solvers)
ru_ad = optimization_ru(A, v2, lp_solvers)
b_ru_wt = @benchmark jump_lp_max_shortest_paths(A, v1, lp_solvers[1]) # 654.919 μs median. 639.01 KiB memory
b_ru_ad = @benchmark jump_lp_max_shortest_paths(A, v2, lp_solvers[1]) # 652.590 μs median. 639.01 KiB memory

# ------------------------------------------------------------------------------

# LINEAR PROGRAMMING MAXIMIZING LONGEST PATHWAYS -------------------------------
# Reference: Ren et al. (2020) DOI:10.1016/j.algal.2019.101767
re_wt = optimization_re(A, v1, lp_solvers)
re_ad = optimization_re(A, v2, lp_solvers)
b_re_wt = @benchmark jump_lp_max_longest_paths(A, v1, lp_solvers[1]) # 646.388 μs median. 639.01 KiB memory
b_re_ad = @benchmark jump_lp_max_longest_paths(A, v2, lp_solvers[1]) # 644.420 μs median. 639.01 KiB memory
# ------------------------------------------------------------------------------

# MIXED INTEGER PROGRAMMING MAXIMIZING SHORTEST ACTIVE PATH---------------------
# Reference: Nookaew et al. (2007) DOI:10.1002/bit.21339
no_wt = optimization_no(A, v1, milp_solvers)
no_ad = optimization_no(A, v2, milp_solvers)
b_no_wt = @benchmark optimization_no(A, v1, milp_solvers) # did not run because license expired.
b_no_ad = @benchmark optimization_no(A, v2, milp_solvers) # did not run because license expired.
# ------------------------------------------------------------------------------

# EXPORT LOG10 RECONSTRUCTION ERRORS -------------------------------------------
# Markov precision over individual fluxes
error_wt_ad_mc_ind = vec(# error is -14.13 and -13.62 for wildtype/disease
  CSV.read(#
    "../data/efm-error-ind-wt-ad-markov.csv", Tables.matrix, header=false
  )
)

# Markov precision over total fluxes
error_wt_ad_mc_tot = vec(# error is -15.35 and -Inf for wildtype/disease
  CSV.read(#
    "../data/efm-error-total-wt-ad-markov.csv", Tables.matrix, header=false
  )
)

# Markov and optimization-based error (wildtype)
combined_error_wt_ind = hcat(#
  [#
    [-error_wt_ad_mc_ind[1]; repeat([""], 9)],
    [""; -sk_wt.efms_raw_error_ind; repeat([""], 7)],
    [""; -or_wt.efms_raw_error_ind; repeat([""], 7)],
    [""; ""; -ru_wt.efms_raw_error_ind],
    [""; ""; -re_wt.efms_raw_error_ind],
    [""; ""; ""; -no_wt.efms_raw_error_ind; repeat([""], 6)]
  ]...
)
combined_error_wt_ind = permutedims(string.(combined_error_wt_ind))
combined_error_wt_ind = hcat(string.(collect(1:6)), combined_error_wt_ind)
CSV.write(#
  "../data/error-wt-ind.csv", Tables.table(combined_error_wt_ind), header=false
)
combined_error_wt_tot = hcat(#
  [#
    [-error_wt_ad_mc_tot[1]; repeat([""], 9)],
    [""; -sk_wt.efms_raw_error_tot; repeat([""], 7)],
    [""; -or_wt.efms_raw_error_tot; repeat([""], 7)],
    [""; ""; -ru_wt.efms_raw_error_tot],
    [""; ""; -re_wt.efms_raw_error_tot],
    [""; ""; ""; -no_wt.efms_raw_error_tot; repeat([""], 6)]
  ]...
)
combined_error_wt_tot = permutedims(string.(combined_error_wt_tot))
combined_error_wt_tot = hcat(string.(collect(1:6)), combined_error_wt_tot)
CSV.write(#
  "../data/error-wt-tot.csv", Tables.table(combined_error_wt_tot), header=false
)

# Markov and optimization-based error (disease)
combined_error_ad_ind = hcat(#
  [#
    [-error_wt_ad_mc_ind[2]; repeat([""], 9)],
    [""; -sk_ad.efms_raw_error_ind; repeat([""], 7)],
    [""; -or_ad.efms_raw_error_ind; repeat([""], 7)],
    [""; ""; -ru_ad.efms_raw_error_ind],
    [""; ""; -re_ad.efms_raw_error_ind],
    [""; ""; ""; -no_ad.efms_raw_error_ind; repeat([""], 6)]
  ]...
)
combined_error_ad_ind = permutedims(string.(combined_error_ad_ind))
combined_error_ad_ind = hcat(string.(collect(1:6)), combined_error_ad_ind)
CSV.write(#
  "../data/error-ad-ind.csv", Tables.table(combined_error_ad_ind), header=false
)
combined_error_ad_tot = hcat(#
  [#
    [-error_wt_ad_mc_tot[2]; repeat([""], 9)],
    [""; -sk_ad.efms_raw_error_tot; repeat([""], 7)],
    [""; -or_ad.efms_raw_error_tot; repeat([""], 7)],
    [""; ""; -ru_ad.efms_raw_error_tot],
    [""; ""; -re_ad.efms_raw_error_tot],
    [""; ""; ""; -no_ad.efms_raw_error_tot; repeat([""], 6)]
  ]...
)
combined_error_ad_tot = permutedims(string.(combined_error_ad_tot))
combined_error_ad_tot = hcat(string.(collect(1:6)), combined_error_ad_tot)
CSV.write(#
  "../data/error-ad-tot.csv", Tables.table(combined_error_ad_tot), header=false
)
# ------------------------------------------------------------------------------

## EXPORT LOG10 EFM WEIGHTS ACROSS ALL METHODS ---------------------------------
combined_wt = [sk_wt, or_wt, ru_wt, re_wt, no_wt]
combined_ad = [sk_ad, or_ad, ru_ad, re_ad, no_ad]
aggregate_efm_values(#
  w_log_wt,
  combined_wt,
  :efms_raw_log,
  "../data/scatterplot-total-efm-weights-wt-log.csv"
)
aggregate_efm_values(#
  w_log_ad,
  combined_ad,
  :efms_raw_log,
  "../data/scatterplot-total-efm-weights-ad-log.csv",
  sortperm(w_log_wt)
)
# ------------------------------------------------------------------------------

## EXPORT LOG10 EFM PROPORTIONS ACROSS ALL METHODS -----------------------------
combined_wt = [sk_wt, or_wt, ru_wt, re_wt, no_wt]
combined_ad = [sk_ad, or_ad, ru_ad, re_ad, no_ad]
aggregate_efm_values(#
  p_log_wt,
  combined_wt,
  :efms_prop_log,
  "../data/scatterplot-total-efm-proport-wt-log.csv"
)
aggregate_efm_values(#
  p_log_ad,
  combined_ad,
  :efms_prop_log,
  "../data/scatterplot-total-efm-proport-ad-log.csv",
  sortperm(p_log_wt)
)
# ------------------------------------------------------------------------------

## EXPORT FREQUENCY OF ZERO EFM WEIGHTS ----------------------------------------
bars_wt = aggregate_efm_zero_freq(w_log_wt, combined_wt)
CSV.write(#
  "../data/barplot-zero-weight-freq-wt.csv",
  Tables.table(hcat(1:size(bars_wt)[1], bars_wt)),
  header = ["Label", "S1", "S2", "S3", "S4", "S5"],
  delim = '\t'
)

bars_ad = aggregate_efm_zero_freq(w_log_ad, combined_ad, sortperm(w_log_wt))
CSV.write(#
  "../data/barplot-zero-weight-freq-ad.csv",
  Tables.table(hcat(1:size(bars_ad)[1], bars_ad)),
  header = ["Label", "S1", "S2", "S3", "S4", "S5"],
  delim = '\t'
)
# ------------------------------------------------------------------------------

## ZERO WEIGHT PROBABILITY -----------------------------------------------------
# Mean probability of all optimization methods/solvers assigning a zero weight
mean_zeros_wt = sum(sum(bars_wt, dims=2)) / size(bars_wt,1) # 0.32
mean_zeros_ad = sum(sum(bars_ad, dims=2)) / size(bars_ad,1) # 0.35

vec(sum(bars_wt .* 21, dims=1)) ./ [2, 2, 8, 8, 1] # average zero weights per method
vec(sum(bars_ad .* 21, dims=1)) ./ [2, 2, 8, 8, 1] # average zero weights per method
# ------------------------------------------------------------------------------

## PERCENTAGE OF FLUX EXPLAINED BY EACH EFM WEIGHT -----------------------------
# First element is Markov solution, rest are the other methods/solvers
num_efms_explain_flux_wt = efm_flux_percentage(w_raw_wt, combined_wt, A, 0.95)
num_efms_explain_flux_ad = efm_flux_percentage(w_raw_ad, combined_ad, A, 0.95)
# Mean of 13.333 EFMs explain 95% of fluxes
sum(num_efms_explain_flux_wt[2:end]) / length(num_efms_explain_flux_wt[2:end])
# Mean of 13 EFMs explain 95% of fluxes
sum(num_efms_explain_flux_ad[2:end]) / length(num_efms_explain_flux_ad[2:end])
# ------------------------------------------------------------------------------


## DIFFERENCE BETWEEN MARKOV WILDTYPE AND DISEASE EFM WEIGHTS ------------------
# Run main-efm-weights-sphingo-markov.jl first to get this text file
dat = CSV.read("../data/table-markov-sphingo.csv", Tables.matrix, header=true)

aa = (w_raw_wt .* vec(sum(A, dims=1)))
aa = aa ./ sum(aa)
sum(sort(aa, rev=true)[1:9]) # Top 9 EFMs in wildtype explain 85.8% of fluxes
id = sortperm(aa, rev=true)
bb = (w_raw_ad .* vec(sum(A, dims=1)))
bb = bb ./ sum(bb)
sum(bb[id][1:9]) # Top 9 EFMs in wildtype explain 86.7% of fluxes of disease

mc_top_wt = efm_flux_percentage_top(w_raw_wt, A, 0.95) # 13 EFMs
mc_top_ad = efm_flux_percentage_top(w_raw_ad, A, 0.95) # 13 EFMs
mc_top_shared = intersect(Set(mc_top_wt), Set(mc_top_ad)) # top 11/13 shared

# Top EFMs unique to wildtype
mc_top_wt[findall(!in(mc_top_shared), mc_top_wt)] # EFMs 18 and 35

# Top EFMs unique to disease
mc_top_wt[findall(!in(mc_top_shared), mc_top_ad)] # EFMs 11 and 13

# Number of extreme fold changes (greater or equal than 2 fold)
f(x,y,i) = log2.((x ./ y))[i]
sum(abs.(f(w_raw_ad, w_raw_wt, id)) .> 2) # 17
id = .!isinf.(f(sk_ad.efms_raw[1], sk_wt.efms_raw[1], sortperm(sk_wt.efms_raw[1])))
sum(abs.(f(sk_ad.efms_raw[1], sk_wt.efms_raw[1], sortperm(sk_wt.efms_raw[1]))[id]) .> 2) # 10
id = .!isinf.(f(or_ad.efms_raw[1], or_wt.efms_raw[1], sortperm(or_wt.efms_raw[1])))
sum(abs.(f(or_ad.efms_raw[1], or_wt.efms_raw[1], sortperm(or_wt.efms_raw[1]))[id]) .> 2) # 8
id = .!isinf.(f(ru_ad.efms_raw[1], ru_wt.efms_raw[1], sortperm(ru_wt.efms_raw[1])))
sum(abs.(f(ru_ad.efms_raw[1], ru_wt.efms_raw[1], sortperm(ru_wt.efms_raw[1]))[id]) .> 2) # 10
id = .!isinf.(f(re_ad.efms_raw[1], re_wt.efms_raw[1], sortperm(re_wt.efms_raw[1])))
sum(abs.(f(re_ad.efms_raw[1], re_wt.efms_raw[1], sortperm(re_wt.efms_raw[1]))[id]) .> 2) # 10
id = .!isinf.(f(no_ad.efms_raw[1], no_wt.efms_raw[1], sortperm(no_wt.efms_raw[1])))
sum(abs.(f(no_ad.efms_raw[1], no_wt.efms_raw[1], sortperm(no_wt.efms_raw[1]))[id]) .> 2) # 6

# Export data (X-axis is EFMs sorted by flux contribution; Y-axis is fold change)
f(x,y,i) = log2.((x ./ y))[i]
g(x,y) = log10.(sort(x .* y))/sum(x .* y)
id = sortperm(w_raw_wt)

res = [#
  g(w_raw_wt, vec(sum(A, dims=1))) f(w_raw_ad, w_raw_wt, id) repeat(["Markov"], length(w_raw_wt));
  g(sk_wt.efms_raw[1], vec(sum(A, dims=1))) f(sk_ad.efms_raw[1], sk_wt.efms_raw[1], sortperm(sk_wt.efms_raw[1])) repeat(["L2_norm_COSMO"], length(w_raw_wt));
  g(or_wt.efms_raw[1], vec(sum(A, dims=1))) f(or_ad.efms_raw[1], or_wt.efms_raw[1], sortperm(or_wt.efms_raw[1])) repeat(["qp_max_spa_COSMO"], length(w_raw_wt));
  g(ru_wt.efms_raw[3], vec(sum(A, dims=1))) f(ru_ad.efms_raw[3], ru_wt.efms_raw[3], sortperm(ru_wt.efms_raw[3])) repeat(["lp_max_spa_SCIP"], length(w_raw_wt));
  g(re_wt.efms_raw[3], vec(sum(A, dims=1))) f(re_ad.efms_raw[3], re_wt.efms_raw[3], sortperm(re_wt.efms_raw[3])) repeat(["lp_min_spa_SCIP"], length(w_raw_wt));
  g(no_wt.efms_raw[1], vec(sum(A, dims=1))) f(no_ad.efms_raw[1], no_wt.efms_raw[1], sortperm(no_wt.efms_raw[1])) repeat(["milp_Gurobi"], length(w_raw_wt));
]


## Number of zero folds because one or both EFM weights were undefined:
## (Manually identified from ..data/scatterplot-fc-flux-contribution.csv
# Markov           |  0
# L2_norm_COSMO:   | 25
# qp_max_spa_COSMO | 27
# lp_max_spa_SCIP  | 32
# lp_min_spa_SCIP  | 34
# milp_Gurobi      | 34

CSV.write(#
  "../data/scatterplot-fc-flux-contribution.csv",
  Tables.table(res),
  header = ["x", "y", "method"],
  delim = '\t'
)

# Figure 4 panel c plot
p1 = plot(g(w_raw_wt, vec(sum(A, dims=1))), f(w_raw_ad, w_raw_wt, id), seriestype = :scatter, legend = :none, markercolor = :red, xlim = (-7, 0), ylim = (-20, 10), xlabel = "EFMs sorted by how much flux they explain log₁₀(%)", ylabel = "Fold change log₂(AD/WT)")
plot!(g.(sk_wt.efms_raw[[1]], Ref(vec(sum(A, dims=1)))), f.(sk_ad.efms_raw[[1]], sk_wt.efms_raw[[1]], sortperm.(sk_wt.efms_raw[[1]])), seriestype = :scatter, markercolor = :blue)
plot!(g.(or_wt.efms_raw[[1]], Ref(vec(sum(A, dims=1)))), f.(or_ad.efms_raw[[1]], or_wt.efms_raw[[1]], sortperm.(or_wt.efms_raw[[1]])), seriestype = :scatter, markercolor = :green)
plot!(g.(ru_wt.efms_raw[[3]], Ref(vec(sum(A, dims=1)))), f.(ru_ad.efms_raw[[3]], ru_wt.efms_raw[[3]], sortperm.(ru_wt.efms_raw[[3]])), seriestype = :scatter, markercolor = :orange)
plot!(g.(re_wt.efms_raw[[3]], Ref(vec(sum(A, dims=1)))), f.(re_ad.efms_raw[[3]], re_wt.efms_raw[[3]], sortperm.(re_wt.efms_raw[[3]])), seriestype = :scatter, markercolor = :purple)
plot!(g.(no_wt.efms_raw, Ref(vec(sum(A, dims=1)))), f.(no_ad.efms_raw, no_wt.efms_raw, sortperm.(no_wt.efms_raw)), seriestype = :scatter, markercolor = :yellow)
plot!([-7; 0], [2.0; 2.0], lw=1, lc=:black, legend=false, line=(:dash))
plot!([-7; 0], [-2.0; -2.0], lw=1, lc=:black, legend=false, line=(:dash))
# ------------------------------------------------------------------------------

