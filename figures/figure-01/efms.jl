using MarkovWeightedEFMs
using MATLAB
using QuantEcon
path = "/home/jchitpin/Documents/PhD/"
include(path * "Code/Julia/MarkovWeightedEFMs.jl/src/higher-order-catalyst.jl")
include(path * "Code/Julia/MarkovWeightedEFMs.jl/src/higher-order-generalization.jl")
dir_matlab = "/home/jchitpin/Documents/PhD/Code/MATLAB/efmtool/"


S = [#
 -1 -1 -1 -1  1  0  0  0  0  0  0  1  1  0  0  1  0  0  1
  1  0  0  0 -1 -1  1  0  0  0  0  0  0  1  0  0  0  0  0
  0  0  0  0  0  1 -1 -1 -1 -1  1  0  0  0  1  0  1  0  0
  0  1  0  0  0  0  0  1  0  0 -1 -1  0  0  0  0  0  0  0
  0  0  1  0  0  0  0  0  1  0  0  0 -1 -1 -1  0  0  0  0
  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0 -1 -1 -1  0
  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1
]
S = Float64.(S)

efms = call_efmtool(S, dir_matlab) # 26... too many

S = [#
 -1 -1 -1  1  0  0  0  0  0  1  1  0  0
  1  0  0 -1 -1  1  0  0  0  0  0  1  0
  0  0  0  0  1 -1 -1 -1  1  0  0  0  1
  0  1  0  0  0  0  1  0 -1 -1  0  0  0
  0  0  1  0  0  0  0  1  0  0 -1 -1 -1
]
S = Float64.(S)
efms = call_efmtool(S, dir_matlab) # 16... too many

S = [#
 -1 -1  1  0  1  0  0  0  1
  1  0 -1 -1  0  1  0  0  0
  0  0  0  1 -1 -1 -1  1  0
  0  1  0  0  0  0  1 -1 -1
]
S = Float64.(S)
efms = call_efmtool(S, dir_matlab) # 8 efms and underdetermined




