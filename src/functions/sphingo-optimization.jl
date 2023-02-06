function optimization_inner(#
  ϕ::Vector{Vector{Float64}},
  A::Matrix{Int64},
  b::Vector{Float64}
)
  # Set any extremely small negative EFM weights to 0
  [ϕ[i][ϕ[i] .<0] .= 0 for i in 1:length(ϕ)]

  # Compute EFM proportions
  efms_prop = [ϕ[i] ./ sum(ϕ[i]) for i in 1:length(ϕ)]

  # Log10 the EFM weights and proportions for plotting purposes
  efms_log10 = ϕ .|> x -> log10.(x)
  efms_prop_log10 = efms_prop .|> x -> log10.(x)

  # Mean squared log10 reconstruction error over individual fluxes
  error_ind = [#
    #log10(sum(((A * ϕ[i] .- b).^2)) / length(ϕ[i])) for i in 1:length(ϕ)
    log10((sum(((A * ϕ[i] .- b).^2)) / length(b))) for i in 1:length(ϕ)
  ]

  # Log10 reconstruction error over total fluxes
  error_tot = [#
    log10(abs(sum(A * ϕ[i]) - sum(b))) for i in 1:length(ϕ)
  ]

  res = (#
    efms_raw=ϕ,
    efms_prop=efms_prop,
    efms_raw_log=efms_log10,
    efms_prop_log=efms_prop_log10,
    efms_raw_error_ind=error_ind,
    efms_raw_error_tot=error_tot
  )
  return res
end

function optimization_sk(A::Matrix{Int64}, b::Vector{Float64}, qp_solvers)
  efms = Vector{Vector{Float64}}()
  for solver in qp_solvers
    push!(efms, jump_qp_l2norm(A, b, solver))
  end
  return optimization_inner(efms, A, b)
end

function optimization_or(A::Matrix{Int64}, b::Vector{Float64}, qp_solvers)
  efms = Vector{Vector{Float64}}()
  for solver in qp_solvers
    push!(efms, jump_qp_max_shortest_paths(A, b, solver))
  end
  return optimization_inner(efms, A, b)
end

function optimization_ru(A::Matrix{Int64}, b::Vector{Float64}, lp_solvers)
  efms = Vector{Vector{Float64}}()
  for solver in lp_solvers
    push!(efms, jump_lp_max_shortest_paths(A, b, solver))
  end
  return optimization_inner(efms, A, b)
end

function optimization_re(A::Matrix{Int64}, b::Vector{Float64}, lp_solvers)
  efms = Vector{Vector{Float64}}()
  for solver in lp_solvers
    push!(efms, jump_lp_max_longest_paths(A, b, solver))
  end
  return optimization_inner(efms, A, b)
end

function optimization_no(A::Matrix{Int64}, b::Vector{Float64}, milp_solvers)
  efms = Vector{Vector{Float64}}()
  for solver in milp_solvers
    push!(efms, jump_milp_max_active_paths(A, b, solver))
  end
  return optimization_inner(efms, A, b)
end

