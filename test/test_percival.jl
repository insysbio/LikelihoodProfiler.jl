using Test, LikelihoodProfiler
include("cases_func.jl")

# test each alg on all functions
function test_percival(
  func_dict::AbstractDict = test_funcs;
  ctol::Float64 = 1e-5,
  kwargs...
)
  @testset "get_endpoint() for CICO_PERCIVAL" begin
    for (f_name, f) in func_dict
      #println("Testing $f_name")
      @testset "Case $f_name" begin
        for i in eachindex(f.x0)
          for (j,dir) in enumerate([:lower, :upper])
            ep = get_endpoint(
              f.x0,
              i,
              x->f.func(x) - f.loss_crit,
              Val(:CICO_PERCIVAL),
              dir;
              ctol,
              kwargs...
            )

            if ep[3] == :BORDER_FOUND_BY_SCAN_TOL
              @test ep[3] == f.status[i][j]
              @test isapprox(ep[1], f.endpoints[i][j], atol=ctol)
            else
              @test ep[3] == :SCAN_BOUND_REACHED
              #@test ep[3] in [:stalled, :SCAN_BOUND_REACHED]
            end
          end
        end
      end
    end
  end
end

ret = get_endpoint(
  [3., 0.5, 8., 2., 2.], # initial parameters' values
  5, # number of parameter to scan
  x->f_5p_3im(x) - 9.0, # lambda(theta) - labmbda_min - delta_lambda

  Val(:CICO_PERCIVAL);
  scan_bound = -9.,

  atol = 1e-6,
  rtol = 1e-6,
  ctol = 1e-5, 

  max_eval = 100000,
  max_time = 30.0,
  max_iter = 2000
)