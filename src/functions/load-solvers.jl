function load_solvers()
  milp_solvers = [# MILP solvers with quadratic constraints
    Gurobi.Optimizer,
  ]
  qp_solvers = [#
    COSMO.Optimizer,
    OSQP.Optimizer
  ]
  lp_solvers = [#
    OSQP.Optimizer,
    Gurobi.Optimizer,
    SCIP.Optimizer,
    CDDLib.Optimizer{Float64},
    ECOS.Optimizer,
    GLPK.Optimizer,
    ProxSDP.Optimizer,
    Tulip.Optimizer
  ]

  return milp_solvers, qp_solvers, lp_solvers
end

