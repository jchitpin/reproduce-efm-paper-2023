function jump_qp_max_shortest_paths(#
  A::Matrix{<:Real},
  v::Vector{<:Real},
  solver
)
  efms_length = sum(A, dims=1)
  model = Model(solver)
  JuMP.@variable(model, x[1:size(A,2)])
  @constraint(model, x .>= 0.0)
  @constraint(model, A * x .- v .== 0)
  @objective(#
    model,
    Min,
    sum((efms_length[i] * x[i])^2 for i = 1:length(x))
  )
  JuMP.optimize!(model)
  return value.(x)
end

function jump_qp_l2norm(#
  A::Matrix{<:Real},
  v::Vector{<:Real},
  solver
)
  model = Model(solver)
  JuMP.@variable(model, x[1:size(A,2)])
  @constraint(model, x .>= 0.0)
  @constraint(model, A * x .- v .== 0)
  @objective(#
    model,
    Min,
    sum(x[i]^2 for i = 1:length(x))
  )
  JuMP.optimize!(model)

  return value.(x)
end

function jump_milp_max_active_paths(#
  A::Matrix{<:Real},
  v::Vector{<:Real},
  solver
)
  model = Model(solver)
  JuMP.@variable(model, x[1:size(A,2)])
  JuMP.@variable(model, beta[1:size(A,2)], Bin)
  @constraint(model, x .>= 0.0)
  @constraint(model, A * (x .* beta) .- v .== 0)
  @objective(#
    model,
    Min,
    sum(beta)
  )
  JuMP.optimize!(model)
  return value.(x)
end

function jump_lp_max_shortest_paths(#
  A::Matrix{<:Real},
  v::Vector{<:Real},
  solver
)
  efms_length = sum(A, dims=1)
  model = Model(solver)
  JuMP.@variable(model, x[1:size(A,2)])
  @constraint(model, x .>= 0.0)
  @constraint(model, A * x .- v .== 0)
  @objective(#
    model,
    Min,
    sum((efms_length[i] * x[i]) for i = 1:length(x))
  )
  JuMP.optimize!(model)
  return value.(x)
end

function jump_lp_max_longest_paths(#
  A::Matrix{<:Real},
  v::Vector{<:Real},
  solver
)
  efms_length = sum(A, dims=1)
  model = Model(solver)
  JuMP.@variable(model, x[1:size(A,2)])
  @constraint(model, x .>= 0.0)
  @constraint(model, A * x .- v .== 0)
  @objective(#
    model,
    Max,
    sum((efms_length[i] * x[i]) for i = 1:length(x))
  )
  JuMP.optimize!(model)
  return value.(x)
end

