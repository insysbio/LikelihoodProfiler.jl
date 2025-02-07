
@recipe function f(sol::PLSolution)
  ls = length(sol)

  layout --> ls
  for i in 1:ls
    @series begin
      subplot := i
      xguide := "x[$(i)]"
      yguide := "likelihood L(x)"
      sol[i]
    end
  end
  return nothing
end

@recipe function f(pv::ProfileValues; steps=true, threshold=true, endpoints=false)

  @series begin
    color --> :blue
    linewidth --> 3
    endpoints ? (label --> "CI interval") : (label --> "profile")
    (pv.x, pv.obj)
  end
  if steps 
    @series begin
      seriestype --> :scatter
      endpoints ? (label --> "CI endpoints") : (label --> "profiler steps")
      (pv.x, pv.obj)
    end
  end 
  if threshold
    obj_level = get_obj_level(pv)
    @series begin
      seriestype --> :hline
      label --> "threshold"
      linewidth --> 2
      linestyle --> :dash
      color --> :green
      [obj_level]
    end
  end
end