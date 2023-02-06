function export_model_to_csv(exportdir)

  # Metabolite names indexed by supplementary information
  mets = hcat(#
    "GA.C1P", "OM.C1P",
    "ER.CER", "IM.CER", "L.CER", "GA.CER", "GACF.CER", "M.CER", "N.CER",
    "OM.CER",
    "GA.GluCER", "GACF.GluCER",
    "L.GSL", "GA.GSL", "OM.GSL",
    "GA.LacCER",
    "C.S1P", "ER.S1P", "IM.S1P", "GA.S1P", "M.S1P", "N.S1P", "OM.S1P",
    "ER.SM", "IM.SM", "L.SM", "GA.SM", "N.SM", "OM.SM",
    "C.SPH", "ER.SPH", "IM.SPH", "L.SPH", "GA.SPH", "M.SPH", "N.SPH", "OM.SPH"
  )

  reaction_pairs = [
    "-" "ER.CER"
    "ER.CER" "ER.SPH"
    "ER.SPH" "ER.CER"
    "ER.SPH" "ER.S1P"
    "ER.S1P" "ER.SPH"
    "ER.S1P" "-"
    "ER.CER" "N.CER"
    "N.CER" "N.SM"
    "N.SM" "N.CER"
    "N.CER" "N.SPH"
    "N.SPH" "N.S1P"
    "N.S1P" "N.SPH"
    "N.SM" "ER.SM"
    "N.SPH" "C.SPH"
    "N.S1P" "C.S1P"
    "C.S1P" "M.S1P"
    "M.S1P" "M.SPH"
    "M.SPH" "M.S1P"
    "M.SPH" "M.CER"
    "M.CER" "M.SPH"
    "M.SPH" "C.SPH"
    "C.SPH" "ER.SPH"
    "ER.SM" "IM.SM"
    "ER.SM" "OM.SM"
    "OM.SM" "OM.CER"
    "OM.CER" "OM.SPH"
    "OM.SPH" "OM.S1P"
    "IM.SM" "IM.CER"
    "IM.CER" "IM.SPH"
    "IM.SPH" "IM.S1P"
    "IM.CER" "OM.CER"
    "IM.SPH" "OM.SPH"
    "IM.S1P" "OM.S1P"
    "IM.S1P" "C.S1P"
    "IM.SPH" "C.SPH"
    "C.SPH" "C.S1P"
    "C.S1P" "ER.S1P"
    "ER.CER" "M.CER"
    "ER.CER" "GA.CER"
    "ER.CER" "GACF.CER"
    "GACF.CER" "GACF.GluCER"
    "GACF.GluCER" "GA.GluCER"
    "GA.GluCER" "GA.LacCER"
    "GA.LacCER" "GA.GSL"
    "GA.GSL" "OM.GSL"
    "OM.GSL" "L.GSL"
    "L.GSL" "L.CER"
    "OM.SM" "L.SM"
    "L.SM" "L.CER"
    "L.CER" "L.SPH"
    "L.SPH" "C.SPH"
    "C.SPH" "GA.SPH"
    "GA.SPH" "GA.S1P"
    "GA.S1P" "GA.SPH"
    "GA.CER" "GA.SM"
    "GA.SM" "OM.SM"
    "GA.SM" "GA.CER"
    "GA.CER" "GA.SPH"
    "GA.CER" "GA.C1P"
    "GA.C1P" "GA.CER"
    "GA.C1P" "OM.C1P"
    "OM.C1P" "OM.CER"
    "OM.CER" "OM.C1P"
    "-" "OM.CER"
    "-" "OM.C1P"
    "-" "OM.CER"
    "-" "OM.SPH"
    "-" "OM.S1P"
    "OM.S1P" "OM.SPH"
  ]

  S = zeros(Int, length(mets), 69)
  for j in 1:size(reaction_pairs, 1)
    if j ∉ [6, 68, 67, 66, 65, 64] # these reactions always have zero flux
      if idx_i(reaction_pairs[j,1]) != 0
        S[idx_i(reaction_pairs[j,1]), j] = -1
      end
      if idx_i(reaction_pairs[j,2]) != 0
        S[idx_i(reaction_pairs[j,2]), j] = 1
      end
    end
  end

  # Remove reaction 1 and 6 from the stoichiometry matrix
  # Their fluxes are zero anyways in wildtype/disease condition
  S[:,1] .= 0 # this reaction breaks the closed-loop network assumption
  S[:,6] .= 0 # this reaction breaks the closed-loop network assumption

  # Export
  CSV.write(#
    join([exportdir, "/stoich-uncorrected.csv"]), Tables.table(S), header=false
  )
  CSV.write(#
    join([exportdir, "/metabolites.csv"]), Tables.table(mets), header=false
  )
end

# Get index of metabolite name
function idx_i(name::String)
  mets = vec(hcat(#
    "GA.C1P", "OM.C1P",
    "ER.CER", "IM.CER", "L.CER", "GA.CER", "GACF.CER", "M.CER", "N.CER",
    "OM.CER",
    "GA.GluCER", "GACF.GluCER",
    "L.GSL", "GA.GSL", "OM.GSL",
    "GA.LacCER",
    "C.S1P", "ER.S1P", "IM.S1P", "GA.S1P", "M.S1P", "N.S1P", "OM.S1P",
    "ER.SM", "IM.SM", "L.SM", "GA.SM", "N.SM", "OM.SM",
    "C.SPH", "ER.SPH", "IM.SPH", "L.SPH", "GA.SPH", "M.SPH", "N.SPH", "OM.SPH"
  ))

  if name == "-"
    return 0
  end
  name = Regex(join(["^", name, "\$"]))
  midx = findall(x -> occursin(name, x), mets)
  if isempty(midx)
    @warn("Check that the metabolite names are spelled correctly.")
    return nothing
  else
    if length(midx) > 1
      @assert(#
        length(midx) <= 1,
        "Check that the metabolite names are spelled correctly."
      )
    else
      return midx[1]
    end
  end
end

