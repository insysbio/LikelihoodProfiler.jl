f_1p(x::Vector{Float64}) = 5.0 + (x[1]-100.0)^2

f_2p_1im(x::Vector{Float64}) = 5.0 + (x[1]-3.0)^2 + 0.0*x[2]

f_2p(x::Vector{Float64}) = 5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2

f_3p_1im(x::Vector{Float64}) = 5.0 + (x[1]-3.0)^2 + (x[2]/x[3]-4.0)^2

f_4p_2im(x::Vector{Float64}) = 5.0 + (x[1]-3.0)^2 + (x[2]-4.0)^2 + 0.0*x[3] + 0.0*x[4]

f_1p_ex(x::Vector{Float64}) = 5.0 + (x[1]-1e-8)^2
