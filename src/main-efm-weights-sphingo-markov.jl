## DESCRIPTION -----------------------------------------------------------------
# This script does the following (for wildtype/disease sphingolipid dataset):
# (1) Enumerate and compute EFM weights by cycle-history Markov chain method.
# (2) Export EFM weights.
# (3) Export EFM matrix.
# (4) Compute and export log10(∑(squared reconstruction error)/|v|)
#     for the Markov solution.
# ------------------------------------------------------------------------------

## USER PARAMETERS -------------------------------------------------------------
# Set working directory containing this script. Example:
cd("/home/jchitpin/Documents/PhD/Projects/reproduce-efm-paper-2023/src/")
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using MarkovWeightedEFMs
using Tables, CSV, NumericIO, JuMP, DelimitedFiles
include.(filter(contains(r".jl$"), readdir("functions"; join=true)))
# ------------------------------------------------------------------------------

## LOAD STOICHIOMETRY AND FLUX DATA --------------------------------------------
stoich_wt = CSV.read("../data/stoich-corrected.csv", Tables.matrix, header=false)
stoich_ad = CSV.read("../data/stoich-corrected.csv", Tables.matrix, header=false)
fluxes_wt = vec(CSV.read("../data/fluxes-wt.csv", Tables.matrix, header=false))
fluxes_ad = vec(CSV.read("../data/fluxes-ad.csv", Tables.matrix, header=false))
mets = vec(CSV.read("../data/metabolites.csv", Tables.matrix, header=false))
# ------------------------------------------------------------------------------

## COMPUTING EFM PROBABILITIES BY MARKOV CHAIN ---------------------------------
res_wt = steady_state_efm_distribution(stoich_wt, fluxes_wt)
res_ad = steady_state_efm_distribution(stoich_ad, fluxes_ad)
# ------------------------------------------------------------------------------

## EXPORTING EFM PROPORTIONS ---------------------------------------------------
# Raw weights
CSV.write(#
  "../data/efm-proport-wt-raw-markov.csv",
  Tables.table(res_wt.p),
  header=false
)
CSV.write(#
  "../data/efm-proport-ad-raw-markov.csv",
  Tables.table(res_ad.p),
  header=false
)
# Log10(weights)
CSV.write(#
  "../data/efm-proport-wt-log-markov.csv",
  Tables.table(log10.(res_wt.p)),
  header=false
)
CSV.write(#
  "../data/efm-proport-ad-log-markov.csv",
  Tables.table(log10.(res_ad.p)),
  header=false
)
# ------------------------------------------------------------------------------

## EXPORTING EFM WEIGHTS -------------------------------------------------------
# Raw weights
CSV.write(#
  "../data/efm-weights-wt-raw-markov.csv",
  Tables.table(res_wt.w),
  header=false
)
CSV.write(#
  "../data/efm-weights-ad-raw-markov.csv",
  Tables.table(res_ad.w),
  header=false
)
# Log10(weights)
CSV.write(#
  "../data/efm-weights-wt-log-markov.csv",
  Tables.table(log10.(res_wt.w)),
  header=false
)
CSV.write(#
  "../data/efm-weights-ad-log-markov.csv",
  Tables.table(log10.(res_ad.w)),
  header=false
)
# ------------------------------------------------------------------------------

## EXPORTING EFM MATRIX --------------------------------------------------------
# EFM matrix (rows are reactions; columns are EFMs)
A_wt = reshape_efm_vector(res_wt.e, stoich_wt)
A_ad = reshape_efm_vector(res_ad.e, stoich_ad)
@assert(A_wt == A_ad)
CSV.write(#
  "../data/efm-matrix-sphingo.csv",
  Tables.table(A_wt),
  header=false
)
# ------------------------------------------------------------------------------

## LOG10 SQUARED RECONSTRUCTION ERROR FOR INDIVIDUAL FLUXES --------------------
# log₁₀(∑((Ax-v)²)/|v|) = -14.08 (wildtype) (with older version of DifferentialEquations -14.13)
error_wt = log10((sum(((A_wt * res_wt.w .- fluxes_wt).^2)) / length(fluxes_wt)))

# log₁₀(∑((Ax-v)²)/|v|) = -13.58 (disease) (with older version of DifferentialEquations -13.62)
error_ad = log10((sum(((A_ad * res_ad.w .- fluxes_ad).^2)) / length(fluxes_ad)))
CSV.write(#
  "../data/efm-error-ind-wt-ad-markov.csv",
  Tables.table([error_wt, error_ad]),
  header=false
)
# ------------------------------------------------------------------------------

## LOG10 SQUARED RECONSTRUCTION ERROR FOR TOTAL FLUXES -------------------------
# log₁₀(|∑(Ax) - ∑(v)|) = -Inf (wildtype) (with older version of DifferentialEquations -15.35)
error_wt_total = log10(abs(sum(A_wt * res_wt.w) - sum(fluxes_wt)))
# log₁₀(|∑(Ax) - ∑(v)|) = -Inf (disease)
error_ad_total = log10(abs(sum(A_ad * res_ad.w) - sum(fluxes_ad)))
CSV.write(#
  "../data/efm-error-total-wt-ad-markov.csv",
  Tables.table([error_wt_total, error_ad_total]),
  header=false
)
# ------------------------------------------------------------------------------

## SUMMARY TABLE FOR MARKOV WEIGHTS FOR WILDTYPE AND DISEASE -------------------
markov_table(res_wt.w, res_ad.w, res_wt.e, mets, "../data/table-markov-sphingo.csv")
matrix2table("../data/table-markov-sphingo.csv", "../data/table-markov-sphingo.tex")
# ------------------------------------------------------------------------------

## EFM METABOLITES/REACTIONS ---------------------------------------------------
efms_mets = [mets[res_wt.e[i]] for i in 1:length(res_wt.e)]
efms_mets = [join(efms_mets[i], " \$\\rightarrow\$ ") for i in 1:length(efms_mets)]
efms_mets = [string(i) * " & " * efms_mets[i] * "\\\\ \\hline" for i in 1:length(efms_mets)]
efms_mets[end] = efms_mets[end][1:53]
writedlm("../data/efm-metabolites.tex", efms_mets)
# ------------------------------------------------------------------------------


