"""
    Setup of the simulation
"""
struct Problem{uType,tType,modelT <: model, fType, jType, oType}
    U0::uType
    dt::Float64
    tspan::tType
    model::modelT
    NTr::Integer
    f::fType
    j::jType
    observable::oType
end

"""
    Integrator containing all informaiton to update the scheme
"""
mutable struct Integrator{uType, fType, jType, ηType ,S<:Solver,R<:Regulators}
    U::uType
    dt::Float64
    f::fType
    j::jType
    η::ηType
    opts::Problem
    alg::S
    Regulators::R

    function Integrator(U0::uType,dt,opts::Problem,alg::algType, regs::R) where {uType <: AbstractGaugeFields, algType, R <: Regulators}
        @unpack f,j = opts
        
        U = copy(U0)

        num_of_basis = U.NC > 1 ? U.NC^2-1 : 1
        η = sqrt(2*dt) .* randn(num_of_basis,U.NV)
        new{uType,typeof(f),typeof(j),typeof(η),algType,R}(U,dt,f,j,η,opts,alg,regs)
    end

    function Integrator(integrator::Integrator{uType,fType,jType,ηType,sType,R}, ff) where {uType,fType,jType,ηType,sType,R}
        new{uType,typeof(ff),jType,ηType,sType,R}(integrator.U,integrator.dt,ff,integrator.j,integrator.η,integrator.opts,integrator.alg,integrator.Regulators)
    end
end
