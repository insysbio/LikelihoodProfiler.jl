module ComponentArraysExt

using LikelihoodProfiler, ComponentArrays

function LikelihoodProfiler._infer_container_labels(x::ComponentArrays.ComponentArray)
  ks = collect(keys(x))
  n = length(x)
  if length(ks) == n && all(k -> (k isa Symbol || k isa AbstractString), ks)
    inferred_labels = [k isa Symbol ? k : Symbol(k) for k in ks]
    allunique(inferred_labels) || return nothing
    return inferred_labels
  end
  return nothing
end

end #module
