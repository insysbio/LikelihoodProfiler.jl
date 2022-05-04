using LikelihoodProfiler, Plots

# functions with zero derivative at one of the endpoints
f1(x) = 3x[1]^4 - 6x[1]^3 + 5
f2(x) = -3(x[1]-5)^4 + 5(x[1]-5)^2 + 5
f3(x) = 3x[1]^3 - 6x[1]^2 + 5
f4(x) = 3x[1]^4 - 6x[1]^3 -  x[1]^2 +  x[1] + 5
f5(x) = 3x[1]^4 - 6x[1]^3 - 3x[1]^2 + 4x[1] + 5

# functions non-differentiable at one of the endpoints
function f6(x)
  if (1.0 <= x[1] < 5.0) 
    return ((x[1]-3)^2 + 5.0)
  else
    return abs(4.5x[1] - 13.5)
  end
end

res1 = get_interval(
  [1.5],
  1,
  f1,
  :CICO_ONE_PASS;
  loss_crit = 5.,
  silent = true
)

res2 = get_interval(
  [5.0],
  1,
  f2,
  :CICO_ONE_PASS;
  loss_crit = 7.,
  silent = true
)

res3 = get_interval(
  [1.5],
  1,
  f3,
  :CICO_ONE_PASS;
  loss_crit = 5.,
  silent = true
)

res4 = get_interval(
  [1.5],
  1,
  f4,
  :CICO_ONE_PASS;
  loss_crit = 5.,
  silent = true
)

res6 = get_interval(
  [3.0],
  1,
  f6,
  :CICO_ONE_PASS;
  loss_crit = 9.,
  silent = true
)