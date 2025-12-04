using Test
using LikelihoodProfiler

# ----------------------------
# ParameterTarget constructors
# ----------------------------
@testset "ParameterTarget constructors" begin
    @test_throws DimensionMismatch ParameterTarget(; idxs=1:2,
        profile_lower=[-10], profile_upper=[10,10])

    @test_throws ArgumentError ParameterTarget(; idxs=[1],
        profile_lower=[7], profile_upper=[4])

    @test_throws ArgumentError ParameterTarget(; idxs=[1,1],
        profile_lower=[-1, -1], profile_upper=[1, 1])

    @test_throws ArgumentError ParameterTarget(; idxs=Int[],
        profile_lower=Float64[], profile_upper=Float64[])

    @test_throws ArgumentError ParameterTarget(; idxs=[0],
        profile_lower=[-1], profile_upper=[1])

    # positive construction with vector bounds 
    pt = ParameterTarget(; idxs=[1,3], profile_lower=[-2.0, -2.0], profile_upper=[2.0, 2.0])
    @test pt.profile_lower == [-2.0, -2.0]
    @test pt.profile_upper == [ 2.0,  2.0]
end

# ----------------------------
# FunctionTarget constructors
# ----------------------------
@testset "FunctionTarget constructors" begin
    f1 = OptimizationFunction((x,p)->x[1])
    f2 = OptimizationFunction((x,p)->x[2])

    # positive construction with vector bounds
    ft = FunctionTarget(; fs=[f1,f2], profile_lower=[-1.0, -1.0], profile_upper=[1.0, 1.0])
    @test ft.profile_lower == [-1.0, -1.0]
    @test ft.profile_upper == [ 1.0,  1.0]

    # length mismatch
    @test_throws DimensionMismatch FunctionTarget(; fs=[f1,f2],
        profile_lower=[-1.0], profile_upper=[1.0, 1.0])

    # non-empty
    @test_throws ArgumentError FunctionTarget(; fs=OptimizationFunction[],
        profile_lower=Float64[], profile_upper=Float64[])
end

# ----------------------------
# ProfileLikelihoodProblem (explicit target)
# ----------------------------
@testset "PL: explicit target" begin
    f = OptimizationFunction((x,p)->x[1]^2 + x[2]^2)
    optprob = OptimizationProblem(f, [1.0, 2.0])

    t_ok = ParameterTarget(; idxs=1:2, profile_lower=[-5.0, -2.0], profile_upper=[4.0, 1.0])

    # u0/optpars length mismatch
    @test_throws DimensionMismatch ProfileLikelihoodProblem(optprob, [0.0, 0.0, 0.0], t_ok)

    # negative threshold
    @test_throws ArgumentError ProfileLikelihoodProblem(optprob, [0.0, 0.0], t_ok; threshold=-1.0)

    # non-finite bounds are rejected at target construction time
    @test_throws ArgumentError ParameterTarget(; idxs=1:2,
        profile_lower=[-Inf, -2.0], profile_upper=[4.0, 1.0])

    # happy path + threshold precedence (threshold wins over conf_level/df)
    prob = ProfileLikelihoodProblem(optprob, [0.0, 0.0], t_ok; conf_level=0.5, df=99, threshold=2.5)
    @test prob.threshold == 2.5
end

# ----------------------------
# PL sugar: parameter profiling
# ----------------------------
@testset "PL sugar: parameters" begin
    f = OptimizationFunction((x,p)->x[1]^2 + x[2]^2 + x[3]^2)
    # No lb/ub in the problem; provide via sugar
    optprob = OptimizationProblem(f, [1.0, 2.0, 3.0])

    # profile a single index with scalar bounds (scalar allowed via sugar)
    prob1 = ProfileLikelihoodProblem(optprob, [0.0, 0.0, 0.0];
        idxs=2, profile_lower=-2.0, profile_upper=3.0)
    @test prob1.target isa ParameterTarget
    @test prob1.target.idxs == [2]
    @test prob1.target.profile_lower == [-2.0]
    @test prob1.target.profile_upper == [ 3.0]

    # profile a subset with per-index bounds (of subset length)
    prob2 = ProfileLikelihoodProblem(optprob, [0.0, 0.0, 0.0];
        idxs=[1,3], profile_lower=[-5.0, -1.0], profile_upper=[2.0, 4.0])
    @test prob2.target.idxs == [1,3]
    @test prob2.target.profile_lower == [-5.0, -1.0]
    @test prob2.target.profile_upper == [ 2.0,  4.0]

    # duplicate idxs rejected (caught in ParameterTarget)
    @test_throws ArgumentError ProfileLikelihoodProblem(optprob, [0.0, 0.0, 0.0];
        idxs=[1,1], profile_lower=[-1.0, -1.0], profile_upper=[1.0, 1.0])

    # initial value must lie within bounds
    @test_throws ArgumentError ProfileLikelihoodProblem(optprob, [0.0, 5.0, 0.0];
        idxs=[2], profile_lower=0.0, profile_upper=4.0)

    # use full-length bounds and slice by idxs
    prob3 = ProfileLikelihoodProblem(optprob, [0.0, 0.0, 0.0];
        idxs=[1,3], profile_lower=[-5.0, -10.0, -1.0], profile_upper=[5.0, 10.0, 1.0])
    @test prob3.target.profile_lower == [-5.0, -1.0]
    @test prob3.target.profile_upper == [ 5.0,  1.0]

    # accept non-Real entries for unselected indices (ignored after slicing)
    full_lb = Any[-5.0, nothing, -1.0]
    full_ub = Any[ 5.0, nothing,  1.0]
    prob4 = ProfileLikelihoodProblem(optprob, [0.0, 0.0, 0.0];
        idxs=[1,3], profile_lower=full_lb, profile_upper=full_ub)
    @test prob4.target.profile_lower == [-5.0, -1.0]
    @test prob4.target.profile_upper == [ 5.0,  1.0]

    # idxs out of range
    @test_throws ArgumentError ProfileLikelihoodProblem(optprob, [0.0, 0.0, 0.0];
        idxs=0, profile_lower=-1.0, profile_upper=1.0)
    @test_throws ArgumentError ProfileLikelihoodProblem(optprob, [0.0, 0.0, 0.0];
        idxs=[1,4], profile_lower=-1.0, profile_upper=1.0)
end

# ----------------------------
# PL sugar: function profiling
# ----------------------------
@testset "PL sugar: functions" begin
    f = OptimizationFunction((x,p)->sum(abs2, x))
    optprob = OptimizationProblem(f, [1.0, 2.0, 3.0])

    g1 = OptimizationFunction((x,p)->x[1] + x[2])
    g2 = OptimizationFunction((x,p)->x[2] - x[3])

    # require bounds
    @test_throws ArgumentError ProfileLikelihoodProblem(optprob, [0.0, 0.0, 0.0], g1)

    # scalar bounds replicate across functions (allowed via sugar)
    probF1 = ProfileLikelihoodProblem(optprob, [0.0, 0.0, 0.0], [g1,g2];
        profile_lower=-2.0, profile_upper=2.0)
    @test probF1.target isa FunctionTarget
    @test length(probF1.target.fs) == 2
    @test probF1.target.profile_lower == [-2.0, -2.0]
    @test probF1.target.profile_upper == [ 2.0,  2.0]

    # vector bounds must match number of functions
    @test_throws DimensionMismatch ProfileLikelihoodProblem(optprob, [0.0, 0.0, 0.0], [g1,g2];
        profile_lower=[-1.0], profile_upper=[1.0, 1.0])
end
