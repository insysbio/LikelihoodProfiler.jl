"""
    struct ProfilePoint
        value::Float64               # x value of profile point
        loss::Float64                # y value of profile point (loss function value)
        params::Array{Float64, 1}    # vector of optimal values of `loss_func` arguments
        ret::Symbol                  # `NLOpt.optimize()` return value 
        counter::Union{Int, Nothing} # number of `loss_func` evaluations
    end

Structure storing one point from profile function.
ret values: `:FORCED_STOP`, `:MAXEVAL_REACHED`, `:FTOL_REACHED`

"""
struct ProfilePoint
    value::Float64
    loss::Float64
    params::Array{Float64, 1}
    ret::Symbol
    counter::Union{Int, Nothing}
end

"""
    struct EndPoint
        value::Union{Real, Nothing}           # value of endpoint or `nothing`
        profilePoints::Array{ProfilePoint, 1} # vector of profile points
        status::Symbol                        # status
        direction::Symbol                     # `:right` or `:left`
        counter::Int                          # number of `loss_func` evaluations
        supreme::Union{Real, Nothing}         # maximal value inside profile interval
    end

Structure storing one endpoint of confidence interval.

status values: `:BORDER_FOUND_BY_SCAN_TOL`, `:BORDER_FOUND_BY_LOSS_TOL`,
 `:SCAN_BOUND_REACHED`, `:MAX_ITER_STOP`, `:LOSS_ERROR_STOP`
"""
struct EndPoint
    value::Union{Real, Nothing}
    profilePoints::Array{ProfilePoint, 1}
    status::Symbol
    direction::Symbol
    counter::Int
    supreme::Union{Real, Nothing}
end
