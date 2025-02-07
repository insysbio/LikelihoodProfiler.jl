# interpolate endpoints when profiling
function interpolate_endpoint(profile_values::ProfileValues)
  # future check if dense
  obj_level = get_obj_level(profile_values)
  interp = LinearInterpolation(profile_values.x[end-1:end], profile_values.obj[end-1:end])
  return interp(obj_level)
end