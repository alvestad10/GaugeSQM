"""
    Setup of the simulation
"""
struct Setup{uType,tType,modelT <: model}
    u0::uType
    dt::Float64
    tspan::tType
    model::modelT
    NTr::Integer
end

"""
    Integrator containing all informaiton to update the scheme
"""
struct Integrator{S<:Solver}
    opts::Setup
    alg::S
end
