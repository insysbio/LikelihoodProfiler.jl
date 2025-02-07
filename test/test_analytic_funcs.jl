using LikelihoodProfiler
using Test, Optimization, OptimizationNLopt, ForwardDiff, OrdinaryDiffEq

const step = 0.3
const atol = step/2

include("problems/analytic_funcs.jl")

function test_plmethod(method, funcs_dict)
  for (_fname,_f) in funcs_dict 
    @testset "$(_fname)" begin
      if !haskey(_f, :grad!) && !haskey(_f, :hess!)
        optf = OptimizationFunction(_f[:func], Optimization.AutoForwardDiff())
      else
        optf = OptimizationFunction(_f[:func]; grad=_f[:grad!], hess=_f[:hess!])
      end
      optprob = OptimizationProblem(optf, _f[:optim])
      plprob = PLProblem(optprob, _f[:optim], _f[:profile_range]; threshold=_f[:threshold])

      sol = profile(plprob, method)
      for i in eachindex(_f[:optim])
        ret = get_retcodes(sol[i])
        ci = get_endpoints(sol[i])
        @test _f[:retcode][i][1] == ret[1] 
        @test _f[:retcode][i][2] == ret[2] 
        _f[:retcode][i][1] == :Identifiable && (@test isapprox(ci[1], _f[:ci][i][1]; atol))
        _f[:retcode][i][2] == :Identifiable && (@test isapprox(ci[2], _f[:ci][i][2]; atol))
      end
    end
  end
end

@testset "Analytic funcs. Fixed-step OptimizationProfiler with derivative-free optimizer" begin

  method = OptimizationProfiler(optimizer = NLopt.LN_NELDERMEAD(), stepper = FixedStep(; initial_step=step))
  test_plmethod(method, funcs_dict)

end

@testset "Analytic funcs. Fixed-step OptimizationProfiler with gradient-based optimizer" begin

  method = OptimizationProfiler(optimizer = Optimization.LBFGS(), stepper = FixedStep(; initial_step=step))
  test_plmethod(method, funcs_dict)

end

@testset "Analytic funcs. IntegrationProfiler with full hessian" begin

  method = IntegrationProfiler(integrator = FBDF(autodiff=false), integrator_opts = (dtmax=step,), matrix_type = :hessian)
  test_plmethod(method, funcs_dict)
  
end

#=
@testset "Analytic funcs. CICOProfiler" begin

  method = CICOProfiler(optimizer = :LN_NELDERMEAD, scan_tol = 1e-3)
  test_plmethod(method, funcs_dict)
  
end
=#