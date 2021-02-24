
f_1p(x) = 5.0 + (x[1]-3.0)^2 # [100., missing, missing, missing]

f_2p_1im(x) = 5.0 + (x[1]-3.0)^2 + 0.0*x[2]  # [3., missing, missing, missing]

f_2p(x) = 5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2  # [3., 4., missing, missing]

f_3p_1im(x) = 5.0 + (x[1]-3.0)^2 + (x[2]/x[3]-4.0)^2 # [3., missing, missing, missing]

f_3p_1im_dep(x) = 5.0 + (x[1]-3.0)^2 + (x[1]-x[2]-1.0)^2 + 0*x[3]^2 # [3., 2., missing]

f_4p_2im(x) = 5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2 + 0.0*x[3] + 0.0*x[4] # [3., 4., missing, missing]

f_4p_3im(x) = 5.0 + (x[1]-3.0)^2 + (x[2]/x[3]-4.0)^2 + 0.0*x[4] # [3., missing, missing, missing]

f_1p_ex(x) = 5.0 + (x[1]-1e-8)^2 # [1e-8, missing, missing, missing]

f_5p_3im(x) = 5.0 + (x[1]-3.0)^2 + (exp(x[2])-1.0)^2 + (x[3]/x[4]-4.0)^2 + 0.0*x[5] # [3., 0., missing, missing, missing]

f_3p_im(x) = 5.0 + (x[1]-3.0)^2 + (exp(x[2])-1.0)^2 + 0.0*x[3] # [3.0, 0., missing]

# test each alg on all functions
function test_alg(alg::Symbol, func_dict::AbstractDict = test_funcs; bounds::Tuple{Float64,Float64} = (-Inf,Inf), tol::Float64 = 1e-2)
  @testset "Testing $alg" begin
    for (f_name, f) in func_dict
      println("Testing $f_name")
      for i in eachindex(f.x0)
        ep = get_interval(
          f.x0,
          i,
          f.func,
          :CICO_ONE_PASS,
          theta_bounds=fill(bounds,length(f.x0)),
          scan_tol=1e-8,
          local_alg = alg,
          loss_crit = f.loss_crit
        )
        @test ep.result[1].status == f.status[i][1]
        @test ep.result[2].status == f.status[i][2]  
        !isa(f.endpoints[i][1], Nothing) && @test isapprox(ep.result[1].value, f.endpoints[i][1], atol=tol)
        !isa(f.endpoints[i][1], Nothing) && @test isapprox(ep.result[2].value, f.endpoints[i][2], atol=tol)
      end
    end
  end
end

# functions dict
test_funcs = Dict(
  :f_1p => (
    func = f_1p,         
    x0 = [3.],
    endpoints = [(1.,5.)], 
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL)],
    loss_crit = 9.
  ),

  :f_2p_1im => (
    func = f_2p_1im,     
    x0 = [3.,1.],
    endpoints = [(1.,5.),
                 (nothing,nothing)], 
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED)],
    loss_crit = 9.
  ),

  :f_2p => (
    func = f_2p,         
    x0 = [3.,4.],
    endpoints = [(1.,5.),
                 (2.,6.)], 
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL)],
    loss_crit = 9.
  ),

  :f_3p_1im => (
    func = f_3p_1im,     
    x0 = [3.,4.,1.1],
    endpoints = [(1.,5.),
                 (nothing,nothing),
                 (nothing,nothing)], 
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED)],
    loss_crit = 9.
  ),

  :f_3p_1im_dep => (
    func = f_3p_1im_dep, 
    x0 = [3., 2., 2.1],
    endpoints = [(1.,5.),
                 (2.0-2.0*sqrt(2.),2.0+2.0*sqrt(2.)),
                 (nothing,nothing)], 
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED)],
    loss_crit = 9.
  ),

  :f_4p_2im => (
    func = f_4p_2im,     
    x0 = [3.,4.,1.,1.],
    endpoints = [(1.,5.),
                 (2.,6.),
                 (nothing,nothing),
                 (nothing,nothing)], 
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED)],
    loss_crit = 9.
  ),

  :f_4p_3im => (
    func = f_4p_3im,     
    x0 = [3.,4.,1.1,1.1],
    endpoints = [(1.,5.),
                 (nothing,nothing),
                 (nothing,nothing),
                 (nothing,nothing)], 
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED)],
    loss_crit = 9.    
  ),

  :f_1p_ex => (
    func = f_1p_ex,      
    x0 = [1.5],
    endpoints = [(-2+1e-8,2+1e-8)], 
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL)],
    loss_crit = 9. 
  ),

  :f_5p_3im => (
    func = f_5p_3im,     
    x0 = [3., 0.5, 8., 2., 2.],
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
    loss_crit = 9.      
  ),

  :f_3p_im => (
    func = f_3p_im,      
    x0 = [3.,1.,1,],
    endpoints = [(1.,5.),
                 (nothing,log(3)),
                 (nothing,nothing)],    
    status = [(:BORDER_FOUND_BY_SCAN_TOL,:BORDER_FOUND_BY_SCAN_TOL),
              (:SCAN_BOUND_REACHED,:BORDER_FOUND_BY_SCAN_TOL),
              (:SCAN_BOUND_REACHED,:SCAN_BOUND_REACHED)],
    loss_crit = 9.
  )
)