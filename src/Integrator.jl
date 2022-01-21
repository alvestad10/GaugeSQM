"""
    Setup of the simulation
"""
struct Problem{uType,tType,modelT <: model,pType}
    U0::uType
    dt::Float64
    tspan::tType
    model::modelT
    NTr::Integer
    p::pType
end

"""
    Integrator containing all informaiton to update the scheme
"""
mutable struct Integrator{uType, pType, ηType ,S<:Solver}
    U::uType
    dt::Float64
    p::pType
    η::ηType
    opts::Problem
    alg::S

    function Integrator(U0::uType,dt,p::pType,opts::Problem,alg::algType) where {uType <: GaugeFields, pType, algType}
        U = copy(U0)

        num_of_basis = U.NC > 1 ? integrator.U.NC^2-1 : 1
        η = sqrt(2*dt) .* randn(num_of_basis,U.NV)
        new{uType,pType,typeof(η),algType}(U,dt,p,η,opts,alg)
    end
end
