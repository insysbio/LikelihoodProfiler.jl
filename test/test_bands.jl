res1 = get_interval(
    [3.0],
    (x)->x[1]^2,
    f_1p,
    :CICO_ONE_PASS;
    scan_bounds=(1e-9,1e9),
    #theta_bounds = [(1e-2,1e2)],
    loss_crit = 8.,
    local_alg = :LN_NELDERMEAD
)
@test isapprox(res1.result[1].value, (3-sqrt(3))^2,atol=1e-2)
@test isapprox(res1.result[2].value, (3+sqrt(3))^2,atol=1e-2)

res2 = get_interval(
    [3.0],
    (x)->x[1]^2,
    f_1p,
    :CICO_ONE_PASS;
    scan_bounds=(1e-9,1e9),
    #theta_bounds = [(1e-2,1e2)],
    loss_crit = 8.,
    local_alg = :LD_MMA
)
@test isapprox(res2.result[1].value, (3-sqrt(3))^2,atol=1e-2)
@test isapprox(res2.result[2].value, (3+sqrt(3))^2,atol=1e-2)

res3 = get_interval(
    [3.0],
    (x) -> log10(x[1]^2), 
    f_1p,
    :CICO_ONE_PASS,
    scan_bounds=(1e-9,1e9),
    #theta_bounds = [(1e-2,1e2)],
    loss_crit = 8.,
    local_alg = :LD_MMA
)
@test isapprox(res3.result[1].value, 2*log10(3-sqrt(3)),atol=1e-2)
@test isapprox(res3.result[2].value, 2*log10(3+sqrt(3)),atol=1e-2)