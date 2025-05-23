using LikelihoodProfiler
using Test, Optimization, OptimizationNLopt, ForwardDiff, OrdinaryDiffEq, CICOBase

const step = 0.3
const atol = step/2

include(joinpath(@__DIR__, "../models/AnalyticFuncs/analytic_funcs.jl"))

function test_plmethod(method, funcs_dict)
  for (_fname,_f) in funcs_dict 
    @testset "$(_fname)" begin
      if !haskey(_f, :grad!) && !haskey(_f, :hess!)
        optf = OptimizationFunction(_f[:func], Optimization.AutoForwardDiff())
      else
        optf = OptimizationFunction(_f[:func]; grad=_f[:grad!], hess=_f[:hess!])
      end
      if !haskey(_f, :p)
        optprob = OptimizationProblem(optf, _f[:optim])
      else
        optprob = OptimizationProblem(optf, _f[:optim], _f[:p])
      end

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

  method = IntegrationProfiler(
    integrator = FBDF(autodiff = AutoFiniteDiff()), 
    integrator_opts = (dtmax=step,), 
    matrix_type = :hessian
  )
  test_plmethod(method, funcs_dict)
  
end



@testset "Analytic funcs. IntegrationProfiler with identity matrix" begin

  method = IntegrationProfiler(
    integrator = FBDF(autodiff = AutoFiniteDiff()), 
    integrator_opts = (dtmax=step,), 
    matrix_type = :identity,
    gamma=1.0
  )
  test_plmethod(method, funcs_dict)
  
end

@testset "Analytic funcs. IntegrationProfiler with identity matrix + reoptimize" begin

  method = IntegrationProfiler(
    integrator = Rosenbrock23(autodiff = AutoFiniteDiff()),
    integrator_opts = (dtmax=0.3,),
    matrix_type = :identity,
    gamma=0.2,            # (!!!) select "bad" gamma
    reoptimize=true,
    optimizer = Optimization.LBFGS()
  )
  test_plmethod(method, funcs_dict)

end

# @testset "Analytic funcs. IntegrationProfiler with Fisher matrix" begin

#   method = IntegrationProfiler(
#     integrator = AutoVern7(Rodas5()), 
#     integrator_opts = (dtmax=step,), 
#     matrix_type = :fisher,
#     gamma=1e-1
#   )
#   test_plmethod(method, funcs_dict)
  
# end


# @testset "Analytic funcs. CICOProfiler" begin

#   method = CICOProfiler(optimizer = :LN_SBPLX, scan_tol = 1e-3)
#   test_plmethod(method, funcs_dict)
  
# end

