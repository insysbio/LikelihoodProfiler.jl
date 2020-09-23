using NLPModels, Percival

function get_interval_percival(
    theta_init::Vector{Float64},
    theta_num::Int,
    loss_func,

    loss_crit::Float64 = 0.0,
)
    ln = length(theta_init)
    m_left = ADNLPModel(x->x[theta_num], theta_init, x->[loss_func(x)], [-Inf], [loss_crit])
    m_right = ADNLPModel(x->-x[theta_num], theta_init, x->[loss_func(x)], [-Inf], [loss_crit])

    return (percival(m_left),percival(m_right))
end
