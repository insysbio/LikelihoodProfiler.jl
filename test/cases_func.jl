
f_1p(x) = 5.0 + (x[1]-3.0)^2 # [3.]

f_2p_1im(x) = 5.0 + (x[1]-3.0)^2 + 0.0*x[2]  # [3., missing]

f_2p(x) = 5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2  # [3., 4.]

f_3p_1im(x) = 5.0 + (x[1]-3.0)^2 + (x[2]/x[3]-4.0)^2 # [3., missing, missing]

f_3p_1im_dep(x) = 5.0 + (x[1]-3.0)^2 + (x[1]-x[2]-1.0)^2 + 0*x[3]^2 # [3., 2., missing]

f_4p_2im(x) = 5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2 + 0.0*x[3] + 0.0*x[4] # [3., 4., missing, missing]

f_4p_3im(x) = 5.0 + (x[1]-3.0)^2 + (x[2]/x[3]-4.0)^2 + 0.0*x[4] # [3., missing, missing, missing]

f_1p_ex(x) = 5.0 + (x[1]-1e-8)^2 # [1e-8, missing]

f_5p_3im(x) = 5.0 + (x[1]-3.0)^2 + (exp(x[2])-1.0)^2 + (x[3]/x[4]-4.0)^2 + 0.0*x[5] # [3., 0., missing, missing, missing]

f_3p_im(x) = 5.0 + (x[1]-3.0)^2 + (exp(x[2])-1.0)^2 + 0.0*x[3] # [3.0, 0., missing]

# test each algorithm on get_interval()
function test_alg_interval(
  alg::NamedTuple,
  func_dict::AbstractDict = test_funcs;
  bounds::Tuple{Float64,Float64} = (-Inf,Inf),
  tol::Float64 = 1e-2,
  kwargs...
)
  @testset "get_interval() for $(alg.algorithm)" begin
    for (f_name, f) in func_dict
      #println("Testing $f_name")
      @testset "Case $f_name" begin
        for i in eachindex(f.x0)
          ep = get_interval(
            f.x0,
            i,
            f.func,
            :CICO_ONE_PASS;
            theta_bounds=fill(bounds,length(f.x0)),
            scan_tol=1e-8,
            local_alg = alg.algorithm,
            loss_crit = f.loss_crit,
            silent = true,
            kwargs...
          )
          should_skip = f_name in alg.skip
          @test ep.result[1].status == f.status[i][1] skip = should_skip
          @test ep.result[2].status == f.status[i][2] skip = should_skip
          if isa(f.endpoints[i][1], Nothing)
            @test isa(ep.result[1].value, Nothing) skip = should_skip
          else
            @test isapprox(ep.result[1].value, f.endpoints[i][1], atol=tol) skip = should_skip
          end
          if isa(f.endpoints[i][2], Nothing)
            @test isa(ep.result[2].value, Nothing) skip = should_skip
          else
            @test isapprox(ep.result[2].value, f.endpoints[i][2], atol=tol) skip = should_skip
          end
        end
      end
    end
  end
end

# test algorithms from the list on get_optimal()
function test_alg_optimal(
  alg::NamedTuple;
  bounds::Tuple{Float64,Float64} = (-Inf,Inf),
  scale::Symbol = :direct,
  scan_tol::Union{Float64,Nothing} = nothing,
  loss_tol::Float64 = 0.,
)
  @testset "get_optimal() for $(alg.algorithm)" begin
    for (f_name, f) in test_funcs
      should_skip = f_name in alg.skip
      @testset "Case $f_name" begin
        result = get_optimal(
          f.x0,
          f.func;
          theta_bounds = fill(bounds, length(f.x0)),
          scale = fill(scale, length(f.x0)),
          scan_tol = scan_tol,
          loss_tol = loss_tol,
          local_alg = alg.algorithm,
          silent = true
        )
        #println("Result: ", result)
        # check loss_tol
        if loss_tol !== 0.
          @test result.ret == :FTOL_REACHED skip = should_skip
          @test isapprox(result.loss, f.loss_optim, atol = loss_tol * 10.) skip = should_skip
        end

        if scan_tol !== nothing
          for i in eachindex(f.x0)
            if f.x_optim !== nothing
              @test result.ret == :XTOL_REACHED skip = should_skip
              #@test isapprox(result.params[i], f.x_optim[i], atol = scan_tol * 10) skip = should_skip
            end
          end
        end

      end
    end
  end
end

# functions dict
test_funcs = Dict(
  :f_1p => (
    func = f_1p,
    x0 = [2.],
    endpoints = [(1.,5.)], 
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL)],
    loss_crit = 9.,
    loss_optim = 5.,
    x_optim = [3.]
  ),

  :f_2p_1im => (
    func = f_2p_1im,
    x0 = [4.,1.],
    endpoints = [(1.,5.),(nothing,nothing)], 
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED)],
    loss_crit = 9.,
    loss_optim = 5.,
    x_optim = [3., nothing]
  ),

  :f_2p => (
    func = f_2p,
    x0 = [4.,5.],
    endpoints = [(1.,5.),
                 (2.,6.)], 
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL)],
    loss_crit = 9.,
    loss_optim = 5.,
    x_optim = [3., 4.]
  ),

  :f_3p_1im => (
    func = f_3p_1im,
    x0 = [4.,4.,1.1],
    endpoints = [(1.,5.),
                 (nothing,nothing),
                 (nothing,nothing)], 
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED)],
    loss_crit = 9.,
    loss_optim = 5.,
    x_optim = [3., nothing, nothing]
  ),

  :f_3p_1im_dep => (
    func = f_3p_1im_dep, 
    x0 = [4., 3., 2.1],
    endpoints = [(1.,5.),
                 (2.0-2.0*sqrt(2.),2.0+2.0*sqrt(2.)),
                 (nothing,nothing)], 
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED)],
    loss_crit = 9.,
    loss_optim = 5.,
    x_optim = [3., 2., nothing]
  ),

  :f_4p_2im => (
    func = f_4p_2im,
    x0 = [4.,5.,1.,1.],
    endpoints = [(1.,5.),
                 (2.,6.),
                 (nothing,nothing),
                 (nothing,nothing)], 
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED)],
    loss_crit = 9.,
    loss_optim = 5.,
    x_optim = [3., 4., nothing, nothing]
  ),

  :f_4p_3im => (
    func = f_4p_3im,
    x0 = [4.,4.,1.1,1.1],
    endpoints = [(1.,5.),
                 (nothing,nothing),
                 (nothing,nothing),
                 (nothing,nothing)], 
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED)],
    loss_crit = 9.,
    loss_optim = 5.,
    x_optim = [3., nothing, nothing, nothing]
  ),

  :f_1p_ex => (
    func = f_1p_ex,
    x0 = [1.5, 2.],
    endpoints = [(-2+1e-8,2+1e-8), (nothing, nothing)], 
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),(:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED)],
    loss_crit = 9.,
    loss_optim = 5.,
    x_optim = [1e-8, nothing]
  ),

  :f_5p_3im => (
    func = f_5p_3im,
    x0 = [4., 0.5, 8., 2., 2.],
    endpoints = [(1.,5.),
                 (nothing,log(3)),
                 (nothing,nothing),
                 (nothing,nothing),
                 (nothing,nothing)], 
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:SCAN_BOUND_REACHED,:BORDER_FOUND_BY_SCAN_TOL),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED)],
    loss_crit = 9.,
    loss_optim = 5.,
    x_optim = [3., nothing, nothing, nothing, nothing]
  ),

  :f_3p_im => (
    func = f_3p_im,
    x0 = [4.,0.5,1,],
    endpoints = [(1.,5.),
                 (nothing,log(3)),
                 (nothing,nothing)],
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:SCAN_BOUND_REACHED,:BORDER_FOUND_BY_SCAN_TOL),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED)],
    loss_crit = 9.,
    loss_optim = 5.,
    x_optim = [3., 0., nothing]
  )
)
