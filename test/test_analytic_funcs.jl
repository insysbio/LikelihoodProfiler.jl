using LikelihoodProfiler
using Test, Optimization, OptimizationNLopt, ForwardDiff, OrdinaryDiffEq, CICOBase

const step = 0.3
const atol = step/2

include(joinpath(@__DIR__, "../models/AnalyticFuncs/analytic_funcs.jl"))

function test_parameter_target(method, funcs_dict)
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
      profile_lower = first.(_f[:profile_range])
      profile_upper = last.(_f[:profile_range])
      plprob = ProfileLikelihoodProblem(optprob, _f[:optim]; profile_lower, profile_upper, threshold=_f[:threshold])

      sol = solve(plprob, method)
      for i in eachindex(_f[:optim])
        ret = retcodes(sol[i])
        ci = endpoints(sol[i])
        @test _f[:retcode][i][1] == ret[1] 
        @test _f[:retcode][i][2] == ret[2] 
        _f[:retcode][i][1] == :Identifiable && (@test isapprox(ci[1], _f[:ci][i][1]; atol))
        _f[:retcode][i][2] == :Identifiable && (@test isapprox(ci[2], _f[:ci][i][2]; atol))
      end
    end
  end
end

function test_function_target(method, funcs_dict)
  f = funcs_dict[:f_2p]
  optf = OptimizationFunction(f[:func], Optimization.AutoForwardDiff())
  optprob = OptimizationProblem(optf, f[:optim])

  g_sum = OptimizationFunction((x, p) -> x[1] + x[2], Optimization.AutoForwardDiff())
  g_diff = OptimizationFunction((x, p) -> x[1] - x[2], Optimization.AutoForwardDiff())
  plprob = ProfileLikelihoodProblem(
    optprob,
    f[:optim],
    [g_sum, g_diff];
    profile_lower = [-20.0, -20.0],
    profile_upper = [20.0, 20.0],
    threshold = f[:threshold]
  )

  sol = solve(plprob, method)

  expected_cis = (
    (7 - sqrt(8.0), 7 + sqrt(8.0)),
    (-1 - sqrt(8.0), -1 + sqrt(8.0))
  )

  for i in eachindex(expected_cis)
    ret = retcodes(sol[i])
    ci = endpoints(sol[i])
    @test ret[1] == :Identifiable
    @test ret[2] == :Identifiable
    @test isapprox(ci[1], expected_cis[i][1]; atol)
    @test isapprox(ci[2], expected_cis[i][2]; atol)
  end

  f_im = funcs_dict[:f_2p_1im]
  optf_im = OptimizationFunction(f_im[:func], Optimization.AutoForwardDiff())
  optprob_im = OptimizationProblem(optf_im, f_im[:optim])

  g_first = OptimizationFunction((x, p) -> x[1], Optimization.AutoForwardDiff())
  g_second = OptimizationFunction((x, p) -> x[2], Optimization.AutoForwardDiff())
  plprob_im = ProfileLikelihoodProblem(
    optprob_im,
    f_im[:optim],
    [g_first, g_second];
    profile_lower = [-20.0, -20.0],
    profile_upper = [20.0, 20.0],
    threshold = f_im[:threshold]
  )

  sol_im = solve(plprob_im, method)

  ret_first = retcodes(sol_im[1])
  ci_first = endpoints(sol_im[1])
  @test ret_first[1] == :Identifiable
  @test ret_first[2] == :Identifiable
  @test isapprox(ci_first[1], 1.0; atol)
  @test isapprox(ci_first[2], 5.0; atol)

  ret_second = retcodes(sol_im[2])
  ci_second = endpoints(sol_im[2])
  @test ret_second[1] == :NonIdentifiable
  @test ret_second[2] == :NonIdentifiable
  @test ci_second[1] === nothing
  @test ci_second[2] === nothing
end

@testset "Analytic funcs. Fixed-step OptimizationProfiler with derivative-free optimizer" begin

  method = OptimizationProfiler(optimizer = NLopt.LN_NELDERMEAD(), stepper = FixedStep(; initial_step=step))
  test_parameter_target(method, funcs_dict)

end

@testset "Analytic funcs. Fixed-step OptimizationProfiler with gradient-based optimizer" begin

  method = OptimizationProfiler(optimizer = Optimization.LBFGS(), stepper = FixedStep(; initial_step=step))
  test_parameter_target(method, funcs_dict)

end

@testset "Analytic funcs. FunctionTarget with CICOProfiler" begin

  method = CICOProfiler(optimizer=:LN_NELDERMEAD)
  test_function_target(method, funcs_dict)

end

#=
@testset "Analytic funcs. Adaptive LineSearchStep OptimizationProfiler with gradient-based optimizer" begin

  method = OptimizationProfiler(optimizer = Optimization.LBFGS(), stepper = LineSearchStep(; initial_step=step, direction=:Secant, linesearch=InterpolationLineSearch()))
  test_parameter_target(method, funcs_dict)

end
=#

@testset "Analytic funcs. IntegrationProfiler with full hessian" begin

  method = IntegrationProfiler(
    integrator = FBDF(autodiff = AutoFiniteDiff()), 
    integrator_opts = (dtmax=step,), 
    matrix_type = :hessian
  )
  test_parameter_target(method, funcs_dict)
  
end

@testset "Analytic funcs. FunctionTarget with IntegrationProfiler" begin

  method = IntegrationProfiler(
    integrator = FBDF(autodiff = AutoFiniteDiff()),
    integrator_opts = (dtmax=step,),
    matrix_type = :hessian
  )
  test_function_target(method, funcs_dict)

end

@testset "Analytic funcs. IntegrationProfiler with identity matrix" begin

  method = IntegrationProfiler(
    integrator = FBDF(autodiff = AutoFiniteDiff()), 
    integrator_opts = (dtmax=step,), 
    matrix_type = :identity,
    gamma=1.0
  )
  test_parameter_target(method, funcs_dict)
  
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
  test_parameter_target(method, funcs_dict)

end

# @testset "Analytic funcs. IntegrationProfiler with Fisher matrix" begin

#   method = IntegrationProfiler(
#     integrator = AutoVern7(Rodas5()), 
#     integrator_opts = (dtmax=step,), 
#     matrix_type = :fisher,
#     gamma=1e-1
#   )
#   test_parameter_target(method, funcs_dict)
  
# end


# @testset "Analytic funcs. CICOProfiler" begin

#   method = CICOProfiler(optimizer = :LN_SBPLX, scan_tol = 1e-3)
#   test_parameter_target(method, funcs_dict)
  
# end

