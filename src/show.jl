function Base.show(io::IO, mime::MIME"text/plain", ep::EndPoint; extra_spaces::Bool = false)
  spaces = extra_spaces ? "    " : ""

  println(io, ":$(ep.direction) Endpoint with status :$(ep.status)")
  println(io, spaces*"  value: $(ep.value)")
  println(io, spaces*"  profilePoints count: $(length(ep.profilePoints))")
  println(io, spaces*"  status: :$(ep.status)")
  println(io, spaces*"  direction: :$(ep.direction)")
  println(io, spaces*"  counter: $(ep.counter)")
  println(io, spaces*"  supreme: $(ep.supreme)")
end

function Base.show(io::IO, mime::MIME"text/plain", pe::ParamInterval)
  left_interval_point = if pe.result[1].status === :BORDER_FOUND_BY_SCAN_TOL || pe.result[1].status === :BORDER_FOUND_BY_LOSS_TOL
    round(pe.result[1].value, sigdigits=4)
  elseif pe.result[1].status === :SCAN_BOUND_REACHED
    bound_1 = pe.input.scan_bounds[1]
    "<$(round(bound_1, sigdigits=4))"
  else
    pe.result[1].value
  end
  right_interval_point = if pe.result[2].status === :BORDER_FOUND_BY_SCAN_TOL || pe.result[2].status === :BORDER_FOUND_BY_LOSS_TOL
    round(pe.result[2].value, sigdigits=4)
  elseif pe.result[2].status === :SCAN_BOUND_REACHED
    bound_2 = pe.input.scan_bounds[2]
    ">$(round(bound_2, sigdigits=4))"
  else
    pe.result[2].value
  end

  println(io, "ParamInterval [$left_interval_point, $right_interval_point] by :$(pe.method)")
  #println(io, "  input: $(pe.input)")
  println(io, "  loss_init: $(pe.loss_init)")
  println(io, "  method: :$(pe.method)")
  println(io, "  result:")
  print(io, "    1: ")
  show(io, mime, pe.result[1], extra_spaces = true)
  print(io, "    2: ")
  show(io, mime, pe.result[2], extra_spaces = true)
end
