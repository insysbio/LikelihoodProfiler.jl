
function garm(x::Float64, scale::Symbol = :direct)
    if scale == :direct
        return x
    elseif scale == :log
        return log10(x)
    else
        throw(DomainError(scale, "scale type is not supported"))
    end
end
garm(x_vec::Vector{Float64}, scale::Symbol = :direct) = garm.(x_vec, scale)
function ungarm(x::Float64, scale::Symbol = :direct)
    if scale == :direct
        return x
    elseif scale == :log
        return exp10(x)
    else
        throw(DomainError(scale, "scale type is not supported"))
    end
end
ungarm(x_vec::Vector{Float64}, scale::Symbol = :direct) = ungarm.(x_vec, scale)

# not used
function ungarm(pp::ProfilePoint, scale::Array{Symbol,1})
    ProfilePoint(
        pp.loss,
        ungarm.(pp.params, scale),
        pp.ret
    )
end
