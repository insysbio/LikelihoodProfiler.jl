module ComponentArraysExt

using LikelihoodProfiler, ComponentArrays

function LikelihoodProfiler._infer_parameter_labels(optpars::ComponentArrays.ComponentArray)
  ks = collect(keys(optpars))
  n = length(optpars)
  if length(ks) == n && all(k -> (k isa Symbol || k isa AbstractString), ks)
    inferred_labels = [k isa Symbol ? k : Symbol(k) for k in ks]
    allunique(inferred_labels) || return nothing
    return inferred_labels
  end
  return nothing
end

end #module
