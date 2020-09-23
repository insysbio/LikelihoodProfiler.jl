
# not working correction needed
res1 = get_endpoint(
    [3.0],
    f_1p,
    (x)->x[1]^2,
    :CICO_ONE_PASS,
    :left;
    loss_crit = 8.
)

res1 = get_endpoint(
    [3.0],
    (x::Vector{Float64}) -> (log10(x[1]^2), f_1p(x)),
    :CICO_ONE_PASS,
    :left;
    loss_crit = 8.
)
