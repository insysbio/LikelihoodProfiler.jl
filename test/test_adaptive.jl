using LineSearch, SciMLBase

f(x,p) = 5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2
optim = [3.,4.]
prob = NonlinearProblem(f, optim)

#=
grad(G,x,p) = begin G[1] = 2.0*(x[1]-3.0); G[2] = 2.0*(x[2]-4.0) end
hess!(H,x,p) = begin H[1,1] = 2.0; H[1,2] = 0.0; H[2,1] = 0.0; H[2,2] = 2.0 end
threshold = 4.0,
ci = [(1.,5.), (2.,6.)]
=#

x = optim
fx = f(x,[])
ls_alg = LiFukushimaLineSearch(; nan_maxiters = nothing)
ls_cache = init(prob, ls_alg, fx, x)

ls_sol = solve!(ls_cache, x, [1,1])
α = ls_sol.step_size