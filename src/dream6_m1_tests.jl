# equality test
out_inq = params_intervals(p0, 8, obj0, obj; logscale = direct_scale, local_alg= :LN_SBPLX,
bounds_params = bounds_params_half_open, bounds_id=bounds_id[8],
constraints_type=:equality, tol_glob=1e-1, tol_loc=1e-3, max_iter=1e5)

# inequality test
out_inq = params_intervals(p0, 8, obj0, obj; logscale = direct_scale, local_alg= :LN_SBPLX,
bounds_params = bounds_params_half_open, bounds_id=bounds_id[8],
constraints_type=:inequality, tol_glob=1e-1, max_iter=1e5)
