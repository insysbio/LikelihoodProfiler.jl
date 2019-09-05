"""
    struct ProfilePoint
        value::Float64
        loss::Float64
        params::Array{Float64, 1}
        ret::Symbol
        counter::Union{Int, Nothing}
    end

Structure storing one point from profile function.
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
        value::Float64
        profilePoints::Array{ProfilePoint, 1}
        status::Symbol
        direction::Symbol
        counter::Int
        supreme::Union{Float64, Nothing}
    end
Structure storing end point for confidence interval.
"""
struct EndPoint
    value::Union{Float64, Nothing}
    profilePoints::Array{ProfilePoint, 1}
    status::Symbol
    direction::Symbol
    counter::Int
    supreme::Union{Float64, Nothing}
end
