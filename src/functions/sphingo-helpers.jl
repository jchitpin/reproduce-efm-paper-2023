# Any zero EFMs that are logged become -Inf. An EFM with a raw weight of 10^-10
# and zero are converted to a very small value (-900) and used in the following
# section to compute the frequency of zero EFM weights.
function aggregate_efm_values(#
  mc::Vector{Float64},
  combined,
  type::Symbol,
  fname::String,
  idx_mc::Vector{Int64}=sortperm(mc)
)
  @assert(length(combined) == 5, "Expected results from 5 objective function.")

  # Order data by Markov and transform and remove infinities
  #idx_mc = sortperm(mc)
  if type ∈ [:efms_raw_log, :efms_prop_log]
    g1(x) = replace(y -> isinf(y) ? -900 : y, x)
    g2(x) = replace(y -> y <= -10 ? -900 : y, x) # very small EFMs are effectively zero and are omitted from this plot
    t(x) = g2(g1(x))[idx_mc]
  else
    @assert(false, "type symbol value incorrect.")
  end

  # Indices to access combined
  if type == :efms_raw
    idx = 1
  elseif type == :efms_prop
    idx = 2
  elseif type == :efms_raw_log
    idx = 3
  elseif type == :efms_prop_log
    idx = 4
  end

  # EFM values
  combined_vals = vcat(#
    [#
      mc[idx_mc];
      t.(combined[1][idx]);
      t.(combined[2][idx]);
      t.(combined[3][idx]);
      t.(combined[4][idx]);
      t.(combined[5][idx]);
    ]...
  )

  # The following is hardcoded
  idx_vals = repeat(1:55, 22)
  efms_idx = repeat(idx_mc, 22)
  methods = [#
    repeat(["Markov"], 55);
    repeat(["L2_norm"], 55*2);
    repeat(["qp_max_spa"], 55*2);
    repeat(["lp_max_spa"], 55*8);
    repeat(["lp_min_spa"], 55*8);
    repeat(["milp"], 55);
  ]
  solvers = [#
    repeat(["Markov"], 55);
    repeat(["COSMO"], 55);
    repeat(["OSQP"], 55);
    repeat(["COSMO"], 55);
    repeat(["OSQP"], 55);
    repeat(["GLPK"], 55);
    repeat(["Gurobi"], 55);
    repeat(["SCIP"], 55);
    repeat(["CDDLib"], 55);
    repeat(["ECOS"], 55);
    repeat(["OSQP"], 55);
    repeat(["ProxSDP"], 55);
    repeat(["Tulip"], 55);
    repeat(["GLPK"], 55);
    repeat(["Gurobi"], 55);
    repeat(["SCIP"], 55);
    repeat(["CDDLib"], 55);
    repeat(["ECOS"], 55);
    repeat(["OSQP"], 55);
    repeat(["ProxSDP"], 55);
    repeat(["Tulip"], 55);
    repeat(["Gurobi"], 55);
  ]

  # Export
  CSV.write(#
    fname,
    Tables.table(#
      hcat(#
        idx_vals, combined_vals, methods, solvers, efms_idx
      )
    ),
    header = ["x", "y", "class", "solver", "efm"],
    delim = '\t'
  )
end

function aggregate_efm_zero_freq(#
  mc::Vector{Float64},
  combined,
  idx_mc::Vector{Int64}=sortperm(mc)
)

  #idx_mc = sortperm(mc)
  g1(x) = replace!(y -> isinf(y) ? -900 : y, x)
  #g2(x) = replace!(y -> y <= -10 ? -900 : y, x) # do not count very small EFMs as zero for this plot!
  s(x) = findall(==(-900), g1(x)[idx_mc])
  function u(x::Vector{Float64})
    y = zeros(length(x))
    y[s(x)] .= 1
    return y
  end

  bars = hcat(#
    sum(u.(combined[1][3])),
    sum(u.(combined[2][3])),
    sum(u.(combined[3][3])),
    sum(u.(combined[4][3])),
    sum(u.(combined[5][3]))
  )
  return bars ./ 21 # 21 solvers across all methods
end

function efm_flux_percentage(#
  mc::Vector{Float64},
  combined,
  efm::Matrix{<:Real},
  cutoff::Float64
)

  @assert(cutoff >= 0.0 && cutoff <= 1.0, "cutoff must be a percentage.")

  # EFM flux percentages (ordered)
  t(x) = x .* vec(sum(efm, dims=1))
  u(x) = x / sum(x)
  v(x) = sort(x, rev=true)
  combined_vals = [#
    [v(u(t(mc)))];
    v.(u.(t.(combined[1].efms_raw)));
    v.(u.(t.(combined[2].efms_raw)));
    v.(u.(t.(combined[3].efms_raw)));
    v.(u.(t.(combined[4].efms_raw)));
    v.(u.(t.(combined[5].efms_raw)));
  ]

  # Top x EFMs explaining |cutoff|
  function f(x)
    i = 0
    j = 0.0
    while j <= cutoff
      i += 1
      j = j + x[i]
    end
    return i
  end
  num_efms_explaining_cutoff = f.(combined_vals)
  return num_efms_explaining_cutoff
end

function efm_flux_percentage_top(#
  mc::Vector{Float64},
  efm::Matrix{<:Real},
  cutoff::Float64
)
  @assert(cutoff >= 0.0 && cutoff <= 1.0, "cutoff must be a percentage.")

  t(x) = x .* vec(sum(efm, dims=1))
  u(x) = x / sum(x)
  v(x) = sort(x, rev=true)
  w(x) = sortperm(x, rev=true)

  function f(x)
    i = 0
    j = 0.0
    while j <= cutoff
      i += 1
      j = j + x[i]
    end
    return i
  end
  top_i_efms = f(v(u(t(mc))))

  return w(mc)[1:top_i_efms] # indices of top EFMs explaining cutoff
end

function compare_efms(#
  matlab::Matrix{Int64},
  julia::Vector{Vector{Int16}},
  S::Matrix{<:Real}
)
  efms_matlab = reshape_efm_matrix(matlab, S)

  print("FluxModeCalculator reports $(length(efms_matlab)) EFMs.\n")
  print("ElementaryFluxModes.jl reports $(length(julia)) EFMs. ")
end

function markov_table(#
  wt_w,
  ad_w,
  efms::Vector{Vector{Int64}},
  mets::Vector{String15},
  fname::String
)
  col_1_idx = string.(collect(1:length(wt_w)))
  col_2_idx = formatted.(wt_w, :SCI, ndigits=4, charset=:ASCII)
  col_3_idx = formatted.(ad_w, :SCI, ndigits=4, charset=:ASCII)
  col_4_idx = formatted.(log2.(ad_w ./ wt_w), :SCI, ndigits=4, charset=:ASCII)
  col_5_idx = ad_w .- wt_w
  col_5_idx = formatted.(ad_w .- wt_w, :SCI, ndigits=4, charset=:ASCII)
  col_6_idx = string.(compartment_spanned(efms, mets))
  col_7_idx = string.((length.(efms) .- 1))

  results = hcat(#
    col_1_idx, col_2_idx, col_3_idx, col_4_idx, col_5_idx, col_6_idx, col_7_idx
  )

  CSV.write(#
    fname,
    Tables.table(#
      results
    ),
    header = ["EFM", "{\$w_{wt}\$}", "{\$w_{ad}\$}", "{\$\\log_{2}(w_{ad}/w_{wt})\$}", "{\$w_{ad} - w_{wt}\$}", "\\# compartments", "EFM length"],
    delim = '\t'
  )
end

function compartment_spanned(efms::Vector{Vector{Int64}}, mets::Vector{String15})
  efms_length = Vector{Int64}(undef, length(efms))
  for i in 1:length(efms)
    efms_length[i] = length(#
      unique([split.(mets[efms[i]], ".")[j][1] for j in 1:length(efms[i])])
    )
  end
  return efms_length
end

