using ADNLPModels, Percival

nlp = ADNLPModel(
    x -> (x[1] - 1)^2 + 100 * (x[2] - x[1]^2)^2,
    [-1.2; 1.0],
    x -> [x[1]^2 + x[2]^2],
    [-Inf],
    [1.0]
)

output = percival(nlp, subsolver=:LN_NELDERMEAD)
println(output)

#=
Generic Execution stats
  status: first-order stationary
  objective value: 0.04567480871343478
  primal feasibility: 4.992273261450464e-11
  dual feasibility: 5.75232680813261e-14
  solution: [0.786415154182609  0.6176983125457493]
  multipliers: [-0.12149655698662709]
  iterations: 6
  elapsed time: 7.631999969482422
=#