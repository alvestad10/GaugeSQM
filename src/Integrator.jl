"""
    Setup of the simulation
"""
struct Problem{uType,tType,modelT <: model, fType, oType}
    U0::uType
    dt::Float64
    tspan::tType
    model::modelT
    NTr::Integer
    f::fType
    observable::oType
end

"""
    Integrator containing all informaiton to update the scheme
"""
mutable struct Integrator{uType, fType, ηType ,S<:Solver,R<:Regulators}
    U::uType
    dt::Float64
    f::fType
    η::ηType
    opts::Problem
    alg::S
    Regulators::R

    function Integrator(U0::uType,dt,opts::Problem,alg::algType, regs::R) where {uType <: GaugeFields, algType, R <: Regulators}
        @unpack f = opts
        
        U = copy(U0)

        num_of_basis = U.NC > 1 ? U.NC^2-1 : 1
        η = sqrt(2*dt) .* randn(num_of_basis,U.NV)
        new{uType,typeof(f),typeof(η),algType,R}(U,dt,f,η,opts,alg,regs)
    end
end
