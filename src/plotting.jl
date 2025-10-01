
@recipe function f(sol::ProfileLikelihoodSolution)
  ls = length(sol)

  layout --> ls
  for i in 1:ls
    @series begin
      subplot := i
      xguide --> "x[$(i)]"
      yguide --> "likelihood L(x)"
      sol[i]
    end
  end
  return nothing
end

@recipe function f(c::ProfileCurve; steps=true, threshold=hasthreshold(c.plprob), endpoints=false)

  @series begin
    color --> :blue
    linewidth --> 3
    endpoints ? (label --> "CI interval") : (label --> "profile")
    (c.x, c.obj)
  end
  if steps 
    @series begin
      seriestype --> :scatter
      endpoints ? (label --> "CI endpoints") : (label --> "profiler steps")
      (c.x, c.obj)
    end
  end 
  if threshold
    obj_level = c.obj_level
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