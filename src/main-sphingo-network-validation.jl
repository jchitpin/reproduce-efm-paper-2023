## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Validate the computational sphingolipid model proposed by Wronowska et al.
#     (Computational modeling of sphingolipid metabolism, BMC Systems Biology,
#     DOI 10.1186/s12918-015-0176-9).
#     Note that the inhibitory kinetic parameters were not provided in the paper
#     and were reverse-engineered to satisfy flux conservation and reported
#     steady state lipid concentrations using the provided kinetic parameters.
#     Note too that the supplementary information incorrectly claimed certain
#     reactions were modelled with mass-action versus Michaelis-Menten kinetics.
#     The correct reaction was chosen based on the specified kinetic parameters
#     provided in the supplementary table.
# (2) Export the stoichiometry matrix and list of metabolites from the model.
#     Note that enzyme-catalyzed reversible reactions are explicitly written as
#     pairs of irreversible reactions with appropriate signs. Reactions modelled
#     by mass-action kinetics are treated as a net one-way reaction because the
#     flux in one direction tends to be orders of magnitude larger than the
#     other. The stoichiometry matrices are identical in the wildtype and
#     disease conditions.
# (3) Compute and export the steady state fluxes and corrected stoichiometry
#     matrix from the computational model of wildtype and disease conditions.
# ------------------------------------------------------------------------------

## USER PARAMETERS -------------------------------------------------------------
# Set working directory containing this script. Example:
cd("/home/jchitpin/Documents/PhD/Projects/reproduce-efm-paper-2023/src/")
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using DifferentialEquations, CSV, Tables, Plots, SBML, JuMP
include.(filter(contains(r".jl$"), readdir("functions"; join=true)))
# ------------------------------------------------------------------------------

## EXPORTING STOICHIOMETRY MATRIX AND METABOLITES ------------------------------
# (1) Must ensure the stoichiometry matrix represents a closed-loop network.
# Note: Reactions 1 and 6 in the stoichiometry matrix break the closed-loop
# assumption. These reactions are set to zero in variable stoich and stoich.csv.
# The wildtype and disease fluxes for these reactions are already zero or were
# set to zero (Vm1 = 0 and not 0.004 for disease state).
export_model_to_csv("../data") # metabolites.csv and stoich-uncorrected.csv
mets = CSV.read("../data/metabolites.csv", Tables.matrix, header=false)
stoich = CSV.read("../data/stoich-uncorrected.csv", Tables.matrix, header=false)
stoich_wt = copy(stoich) # stoichiometry matrix identical between conditions
stoich_ad = copy(stoich)
# ------------------------------------------------------------------------------

## VALIDATING STEADY STATE CONCENTRATIONS OF ODE MODEL -------------------------
# Initial concentrations (wildtype steady state from supplementary information)
u0 = ss_wt_conc()
tspan = (0.0,1000.0)

# Construct system of ODEs for the wildtype model
prob_wt = DifferentialEquations.ODEProblem(model_wt!, u0, tspan)
sol_wt = DifferentialEquations.solve(prob_wt);
Plots.plot(#
  sol_wt,
  tspan = tspan,
  legendfontsize = 4,
  label = mets,
  title = "Time series of lipid family concentrations",
  xlabel = "Time (min)",
  ylabel = "Concentration (nmol/mg)",
  framestyle = :box
)

# Compare wildtype conc. in supp. info (left) with numerical values (right)
hcat(1:length(u0), u0, sol_wt[end])

# Construct system of ODEs for the disease model
prob_ad = DifferentialEquations.ODEProblem(model_ad!, u0, tspan)
sol_ad = DifferentialEquations.solve(prob_ad);

Plots.plot(#
  sol_ad,
  tspan = tspan,
  legendfontsize = 4,
  label = mets,
  title = "Time series of lipid family concentrations",
  xlabel = "Time (min)",
  ylabel = "Concentration (nmol/mg)",
  framestyle = :box
)
# ------------------------------------------------------------------------------

## EXPORTING WILDTYPE AND DISEASE FLUXES ---------------------------------------
# Remove all reactions in the stoichiometry matrices containing all zeroes
# These reactions contain zero fluxes in the kinetic models for both conditions
idx_remove = findall([all(stoich[:,j] .== 0) for j in 1:size(stoich,2)])
stoich_clean = stoich[:,setdiff(1:size(stoich,2), idx_remove)]
stoich_wt = copy(stoich_clean)
stoich_ad = copy(stoich_clean)

# Wildtype fluxes
fluxes_wt = model_wt_fluxes(sol_wt[end])
fluxes_wt_heatmap = deepcopy(abs.(fluxes_wt))
fluxes_wt = fluxes_wt[setdiff(1:length(fluxes_wt), idx_remove)]

# Disease fluxes
fluxes_ad = model_ad_fluxes(sol_ad[end])
fluxes_ad_heatmap = deepcopy(abs.(fluxes_ad))
fluxes_ad = fluxes_ad[setdiff(1:length(fluxes_ad), idx_remove)]

# (1) Export fluxes for heatmap
x = repeat(1:length(fluxes_wt_heatmap), 2)
y1 = repeat([1], length(fluxes_wt_heatmap))
y2 = repeat([2], length(fluxes_wt_heatmap))
y = [y1; y2]
# Re-order fluxes by compartment
c1 = sort([2, 3, 4, 5, 6]) # ER
c2 = sort([8, 9, 10, 11, 12]) # nucleus
c3 = sort([17, 18, 19, 20]) # mito
c4 = sort([25, 26, 27, 28, 29, 30, 31, 32, 33, 46, 48, 62, 63, 64, 65, 66, 67, 68, 69]) # IM/OM
c5 = sort([41, 42, 43, 44, 52, 53, 54, 55, 57, 58, 59, 60]) # GA and GACF
c6 = sort([47, 49, 50]) # Lysosome
c7 = sort([56, 24, 23, 37, 22, 16, 36, 34, 15, 35, 14, 7, 13, 21, 61, 1, 38, 39, 40, 45, 51]) # cell
compart = [c1; c2; c3; c4; c5; c6; c7]
# Construct heatmap file
heat = hcat(x, y, [fluxes_wt_heatmap[compart]; fluxes_ad_heatmap[compart]])
heat_log10 = deepcopy(heat)
rm_idx = heat_log10[:, 3] .== 0
heat_log10[.!rm_idx,3] = log10.(heat_log10[.!rm_idx,3])
CSV.write("../data/fluxes-heatmap-log10-scale.csv", Tables.table(heat_log10), header = [:x, :y, :z])

# (2) Must ensure fluxes are positive for computing transition probabilities
# by flipping the flux signs and reaction directionalities in the stoichiometry
# matrix.
# Wildtype fluxes
idx_neg_wt = findall(<(0), fluxes_wt) # reaction indices with negative fluxes
stoich_wt[:,idx_neg_wt] = -1 * stoich_wt[:,idx_neg_wt] # flip reaction directions
fluxes_wt[idx_neg_wt] .= -1 * fluxes_wt[idx_neg_wt] # negative fluxes now positive
CSV.write("../data/fluxes-wt.csv", Tables.table(fluxes_wt), header = false)

# Disease fluxes
idx_neg_ad = findall(<(0), fluxes_ad) # reaction indices with negative fluxes
stoich_ad[:,idx_neg_ad] = -1 * stoich_ad[:,idx_neg_ad] # flip reaction directions
fluxes_ad[idx_neg_ad] .= -1 * fluxes_ad[idx_neg_ad] # flip negatives fluxes
CSV.write("../data/fluxes-ad.csv", Tables.table(fluxes_ad), header = false)

# Wildtype/disease stoichiometry matrix with corrected reaction directions
@assert(stoich_wt == stoich_ad)
CSV.write("../data/stoich-corrected.csv", Tables.table(stoich_wt), header = false)
# ------------------------------------------------------------------------------


