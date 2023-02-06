## DESCRIPTION -----------------------------------------------------------------
# This script does the following:
# (1) Installs all Julia packages necessary to run the scripts.
#     Note: GLMakie may return an error if OpenGL is not installed.
#     See https://github.com/JuliaPlots/GLMakie.jl for troubleshooting.
#     It is only required for one optional plotting function in
#     MarkovWeightedEFMs.jl.
#     Note: Gurobi is the only licensed solver that must be manually installed
# ------------------------------------------------------------------------------

## INSTALLING MarkovWeightedEFMs.jl -------------------------------------------
using Pkg
Pkg.add(url="https://github.com/jchitpin/MarkovWeightedEFMs.jl.git")
# ------------------------------------------------------------------------------

## INSTALLING PACKAGE DEPENDENCIES ---------------------------------------------
using Pkg
Pkg.add(#
  [
     "JuMP", # mathematical optimization interface
     "GLPK", "Gurobi", "SCIP", "COSMO", # solvers
     "OSQP", "CDDLib", "ECOS", "ProxSDP", "Tulip", # solvers
     "CSV", "Tables", # importing/exporting data as text files/matrices
     "DelimitedFiles", # exporting data
     "RowEchelon", # for some matrix property calculations
     "PrettyTables", # printing formatted tables to REPL
     "DifferentialEquations", "SBML", "Plots", # validating sphingolipid model
     "BenchmarkTools", # benchmarking
     "GLMakie", # plotting backend for MarkovWeightedEFMs
     "NumericIO" # formatting numbers for printing tables
  ]
)
# ------------------------------------------------------------------------------

