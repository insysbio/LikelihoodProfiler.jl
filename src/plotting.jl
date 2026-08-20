
@recipe function f(sol::ProfileLikelihoodSolution; xtransform=identity)
  nprofiles = length(sol)
  lbls = profile_labels(sol)

  ncols = min(nprofiles, 3)
  nrows = cld(nprofiles, ncols)

  layout --> (nrows, ncols)
  size --> (360 * ncols, 300 * nrows)
  legend --> (nprofiles == 1 ? :topright : false)
  margin --> (3, :mm)

  for i in 1:nprofiles
    xlbl = isnothing(lbls) ? "x[$(i)]" : string(lbls[i])
    @series begin
      subplot := i
      xguide --> xlbl
      yguide --> "objective function"
      xtransform := xtransform
      sol[i]
    end
  end
  return nothing
end

@recipe function f(c::ProfileCurve; steps=true, threshold=hasthreshold(c.plprob), endpoints=false, xtransform=identity)

  xvals = map(xtransform, c.x)

  @series begin
    color --> :blue
    linewidth --> 3
    endpoints ? (label --> "CI interval") : (label --> "profile")
    (xvals, c.obj)
  end
  if steps 
    @series begin
      seriestype --> :scatter
      endpoints ? (label --> "CI endpoints") : (label --> "profiler steps")
      (xvals, c.obj)
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
