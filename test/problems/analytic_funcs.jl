
funcs_dict = Dict(

#=
  :f_1p => Dict(
    :func => (x,p) -> 5.0 + (x[1]-3.0)^2,
    :optim => [3.],
    :threshold => 4.0,
    :ci => [(1.,5.)],
    :profile_range => [(-20.,20.)],
    :retcode => [(:Identifiable,:Identifiable)]
  ),
=#
  :f_1p_ex => Dict(
    :func => (x,p) -> 5.0 + (x[1]-1e-8)^2 + 0.0*x[2], 
    :optim => [1e-8, 2.],
    :grad! => (G,x,p) -> begin G[1] = 2.0*(x[1]-1e-8); G[2] = 0.0 end,
    :hess! => (H,x,p) -> begin H[1,1] = 2.0; H[1,2] = 0.0; H[2,1] = 0.0; H[2,2] = 0.0 end,
    :threshold => 4.0,
    :ci => [(-2+1e-8,2+1e-8), 
                   (nothing, nothing)],
    :profile_range => [(-20.,20.), (-20.,20.)],
    :retcode => [(:Identifiable,:Identifiable),
                (:NonIdentifiable,:NonIdentifiable)]
  ),

  :f_2p_1im => Dict(
    :func => (x,p) -> 5.0 + (x[1]-3.0)^2 + 0.0*x[2], 
    :grad! => (G,x,p) -> begin G[1] = 2.0*(x[1]-3.0); G[2] = 0.0 end,
    :hess! => (H,x,p) -> begin H[1,1] = 2.0; H[1,2] = 0.0; H[2,1] = 0.0; H[2,2] = 0.0 end,
    :optim => [3.,1.],
    :threshold => 4.0,
    :ci => [(1.,5.),
                   (nothing,nothing)],
    :profile_range => [(-20.,20.),(-20.,20.)],
    :retcode => [(:Identifiable,:Identifiable),
                (:NonIdentifiable,:NonIdentifiable)]
  ),

  :f_2p => Dict(
    :func => (x,p) -> 5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2,
    :grad! => (G,x,p) -> begin G[1] = 2.0*(x[1]-3.0); G[2] = 2.0*(x[2]-4.0) end,
    :hess! => (H,x,p) -> begin H[1,1] = 2.0; H[1,2] = 0.0; H[2,1] = 0.0; H[2,2] = 2.0 end,
    :optim => [3.,4.],
    :threshold => 4.0,
    :ci => [(1.,5.),
                   (2.,6.)],
    :profile_range => [(-20.,20.),(-20.,20.)],
    :retcode => [(:Identifiable,:Identifiable),
                (:Identifiable,:Identifiable)]
  ),

  :f_3p_1im => Dict(
    :func =>  (x,p) -> 5.0 + (x[1]-3.0)^2 + (x[1]-x[2]-1.0)^2 + 0*x[3]^2,
    :optim => [3., 2., 2.1],
    :grad! => (G,x,p) -> begin G[1] = 2.0*(x[1]-3.0) + 2.0*(x[1]-x[2]-1.0); G[2] = -2.0*(x[1]-x[2]-1.0); G[3] = 0.0 end,
    :hess! => (H,x,p) -> begin H[1,1] = 4.0; H[1,2] = -2.0; H[1,3] = 0.0; H[2,1] = -2.0; H[2,2] = 2.0; H[2,3] = 0.0; H[3,1] = 0.0; H[3,2] = 0.0; H[3,3] = 0.0 end,
    :threshold => 4.0,
    :ci => [(1.,5.),
                  (2.0-2.0*sqrt(2.),2.0+2.0*sqrt(2.)),
                  (nothing,nothing)], 
    :profile_range => [(-20.,20.),(-20.,20.),(-20.,20.)],
    :retcode => [(:Identifiable,:Identifiable),
                (:Identifiable,:Identifiable),
                (:NonIdentifiable,:NonIdentifiable)]
  ),

  :f_3p_1im2 => Dict(
    :func => (x,p) -> 5.0 + (x[1]-3.0)^2 + (exp(x[2])-1.0)^2 + 0.0*x[3], 
    :optim => [3.,0.,1,],
    :grad! => (G,x,p) -> begin f(x) = 5.0 + (x[1]-3.0)^2 + (exp(x[2])-1.0)^2 + 0.0*x[3]; ForwardDiff.gradient!(G, f, x) end,
    :hess! => (H,x,p) -> begin f(x) = 5.0 + (x[1]-3.0)^2 + (exp(x[2])-1.0)^2 + 0.0*x[3]; ForwardDiff.hessian!(H, f, x) end,
    :threshold => 4.0,
    :ci => [(1.,5.),
                   (nothing,log(3)),
                   (nothing,nothing)],
    :profile_range => [(-20.,20.),(-20.,20.),(-20.,20.)],
    :retcode => [(:Identifiable,:Identifiable),
                (:NonIdentifiable,:Identifiable),
                (:NonIdentifiable,:NonIdentifiable)]
  ),
#=
  :f_3p_2im => Dict(
    :func => (x,p) -> 5.0 + (x[1]-3.0)^2 + (x[2]/x[3]-4.0)^2, 
    :optim => [3.,4.,1.1],
    :threshold => 4.0,
    :ci => [(1.,5.),
                   (nothing,nothing),
                   (nothing,nothing)], 
    :retcode => [(:Identifiable,:Identifiable),
                (:NonIdentifiable,:NonIdentifiable),
                (:NonIdentifiable,:NonIdentifiable)]
  ),
=#
  :f_4p_2im => Dict(
    :func => (x,p) -> 5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2 + 0.0*x[3] + 0.0*x[4], 
    :optim => [3.,4.,1.,1.],
    :grad! => (G,x,p) -> begin f(x) = 5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2 + 0.0*x[3] + 0.0*x[4]; ForwardDiff.gradient!(G, f, x) end,
    :hess! => (H,x,p) -> begin f(x) = 5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2 + 0.0*x[3] + 0.0*x[4]; ForwardDiff.hessian!(H, f, x) end,
    :threshold => 4.0,
    :ci => [(1.,5.),
                   (2.,6.),
                  (nothing,nothing),
                  (nothing,nothing)],
    :profile_range => [(-20.,20.),(-20.,20.),(-20.,20.),(-20.,20.)],
    :retcode => [(:Identifiable,:Identifiable),
                (:Identifiable,:Identifiable),
                (:NonIdentifiable,:NonIdentifiable),
                (:NonIdentifiable,:NonIdentifiable)]
  ),
  :rosenbrock => Dict(
    :func => (x,p) -> 5.0 + (1.0 - x[1])^2 + 100.0*(x[2] - x[1]^2)^2, 
    :optim => [1.,1.],
    :threshold => 4.0,
    :ci => [(-1.,3.),
            (-0.174,9.)],
    :profile_range => [(-10.,10.),(-10.,10.)],
    :retcode => [(:Identifiable,:Identifiable),
                (:Identifiable,:Identifiable)]
  ),
#=
  :f_4p_3im => Dict(
    :func => (x,p) -> 5.0 + (x[1]-3.0)^2 + (x[2]/x[3]-4.0)^2 + 0.0*x[4], 
    :optim => [3.,4.,1.1,1.1],
    :threshold => 4.0,
    :ci => [(1.,5.),
                   (nothing,nothing),
                   (nothing,nothing),
                   (nothing,nothing)],
    :retcode => [(:Identifiable,:Identifiable),
                (:NonIdentifiable,:NonIdentifiable),
                (:NonIdentifiable,:NonIdentifiable),
                (:NonIdentifiable,:NonIdentifiable)]
  ),
=#
#=
  :f_5p_3im => Dict(
    :optim => [3., 0.5, 8., 2., 2.],  
    :func => (x,p) -> 5.0 + (x[1]-3.0)^2 + (exp(x[2])-1.0)^2 + (x[3]/x[4]-4.0)^2 + 0.0*x[5], 
    :threshold => 4.0,
    :ci => [(1.,5.),
                   (nothing,log(3)),
                   (nothing,nothing),
                   (nothing,nothing),
                   (nothing,nothing)],
    :retcode => [(:Identifiable,:Identifiable),
                (:NonIdentifiable,:Identifiable),
                (:NonIdentifiable,:NonIdentifiable),
                (:NonIdentifiable,:NonIdentifiable),
                (:NonIdentifiable,:NonIdentifiable)]
  )
=#
)

