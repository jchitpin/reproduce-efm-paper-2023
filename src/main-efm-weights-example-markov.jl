## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Enumerate and compute EFM weights in the two example networks by
#     cycle-history Markov chain method.
# (2) Export EFM matrix.
# (3) Solutions are manually aggregated and inputted into the paper tables.
# ------------------------------------------------------------------------------

## USER PARAMETERS -------------------------------------------------------------
# Set working directory containing this script. Example:
cd("/home/jchitpin/Documents/PhD/Projects/reproduce-efm-paper-2023/src/")
# ------------------------------------------------------------------------------

## JULIA PACKAGES AND FUNCTIONS ------------------------------------------------
using MarkovWeightedEFMs
using Tables, CSV
using GLMakie
# ------------------------------------------------------------------------------

## EXAMPLE 01 ------------------------------------------------------------------
# Stoichiometry matrix (rows are metabolites; cols are reactions)
S1 = [#
     -1 -1  0  0  0  0  0  1
      1  0 -1  0  0  0  0  0
      0  1  1 -1 -1  0  0  0
      0  0  0  1  0 -1  0  0
      0  0  0  0  1  0 -1  0
      0  0  0  0  0  1  1 -1
]
v1 = [2, 2, 2, 2, 2, 2, 2, 4]

all(S1 * v1 .== 0) # steady state requirement met


# Enumerate EFMs, compute their probabilities and weights
res1 = steady_state_efm_distribution(S1, v1)
res1.e # EFM state sequences
res1.p # EFM probabilities (in same order as EFMs)
res1.w # EFM weights (in same order as EFMs)

# Plot cycle-history Markov chain
T1 = stoich_to_transition(S1, v1)
GLMakie.activate!() # OpenGL backend for plotting
tree_plot(#
  T1,
  1, # root state. Arbitrary/optional and does not change EFM weights
  show_all = true # optional. Shows cycles from green leaves to root
)

# Reshape EFM sequences into matrix and export (for optimization methods)
A1 = reshape_efm_vector(res1.e, S1)
CSV.write(#
  "../data/efm-matrix-example-1.csv",
  Tables.table(A1),
  header=false
)
# ------------------------------------------------------------------------------

## EXAMPLE 02 ------------------------------------------------------------------
# Stoichiometry matrix (rows are metabolites; cols are reactions)
S2 = [#
     -1 -1  1  0  1  0  0  0  1
      1  0 -1 -1  0  1  0  0  0
      0  0  0  1 -1 -1 -1  1  0
      0  1  0  0  0  0  1 -1 -1
]

# Steady state fluxes
v2 = [3, 3, 2, 3, 2, 2, 2, 3, 2];

# Enumerate EFMs, compute their probabilities and weights
res2 = steady_state_efm_distribution(S2, v2);
res2.e # EFM state sequences
res2.p # EFM probabilities (in same order as EFMs)
res2.w # EFM weights (in same order as EFMs)

A = reshape_efm_vector(res2.e, S2);
A * res2.w - v2

# Plot cycle-history Markov chain
T2 = stoich_to_transition(S2, v2)
GLMakie.activate!() # OpenGL backend for plotting
tree_plot(#
  T2,
  1, # root state. Arbitrary/optional and does not change EFM weights
  show_all = true # optional. Shows cycles from green leaves to root
)

# Reshape EFM sequences into matrix and export (for optimization methods)
A2 = reshape_efm_vector(res2.e, S2)
CSV.write(#
  "../data/efm-matrix-example-2.csv",
  Tables.table(A2),
  header=false
)
# ------------------------------------------------------------------------------





