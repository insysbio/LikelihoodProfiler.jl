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

#=
@recipe function f(::Type{Array{ProfilePoint,1}}, pps::Array{ProfilePoint,1})
    x = [pps[i].value for i in 1:length(pps)]
    y = [pps[i].loss for i in 1:length(pps)]
    plot!(x, y)
end
=#

function plot2(pps::Array{ProfilePoint,1})
    x = [pps[i].value for i in 1:length(pps)]
    y = [pps[i].loss for i in 1:length(pps)]
    plot(x, y)
end

"""
    struct EndPoint
        value::Float64
        profilePoints::Array{ProfilePoint, 1}
        status::Symbol
        direction::Symbol
        counter::Int
    end
Structure storing end point for confidence interval.
"""
struct EndPoint
    value::Union{Float64, Nothing}
    profilePoints::Array{ProfilePoint, 1}
    status::Symbol
    direction::Symbol
    counter::Int
end

function plot2(ep::EndPoint)
    # plot profile points
    p = plot2(ep.profilePoints)
    # plot endpoint
    #plot(p, [ep.value], [ep.profilePoints[1].loss], seriestype=:scatter)
end
