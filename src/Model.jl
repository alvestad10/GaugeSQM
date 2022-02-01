
abstract type model end

"""
Polyakov chain model; 1-dimensional
"""
struct PolyakovChainModel{T <: SUn,βType,β2Type} <: model
    β::βType
    N::Integer
    κ::Float64
    μ::Float64
    β₁::β2Type
    β₂::β2Type
end




function f_chain(V::U1AlgebraFields,U::GaugeFields{U1,eType},dt,η,p::PolyakovChainModel{U1,βType}) where {eType,βType}
    for j in 1:U.NV
        M::U1Matrix = U[j]
        for i in vcat(j+1:U.NV, 1:j-1)
            M *= U[i]
        end

        MInv::U1Matrix = inv(U[j])
        for i in vcat(j-1:-1:1, U.NV:-1:j+1)
            MInv *= inv(U[i])
        end

        trT!(view(V.a,:,j),p.β₁*M + p.β₂*MInv)
    end
    
    trT!(V,U)
    muladd!(V,dt,η)
end

function f_chain(V::SU2AlgebraFields,U::GaugeFields{SU2,eType},dt,η,p::PolyakovChainModel{SU2,βType}) where {eType,βType}
    trT!(V,U)
    muladd!(V,dt*im*(p.β/2),η)
end

function f_chain(V::SU3AlgebraFields,U::GaugeFields{SU3,eType},dt,η,p::PolyakovChainModel{SU3,βType}) where {eType,βType}
    for j in 1:U.NV
        M::SU3Matrix = U[j]
        for i in vcat(j+1:U.NV, 1:j-1)
            M *= U[i]
        end

        MInv::SU3Matrix = inv(U[j])
        for i in vcat(j-1:-1:1, U.NV:-1:j+1)
            MInv *= inv(U[i])
        end

        trT!(view(V.a,:,j),p.β₁*M + p.β₂*MInv)
    end
    
    trT!(V,U)
    muladd!(V,dt,η)
end

function observable_chain(U::GaugeFields{SU2,eType},model::PolyakovChainModel{SU2,βType}) where {eType,βType}
    return (model.β/2) * tr(U)
end

function observable_chain(U::GaugeFields{SU3,eType},model::PolyakovChainModel{SU3,βType}) where {eType,βType}
    @unpack β, κ, μ = model
    return (β + κ*exp(μ)) * tr(U) + (β + κ*exp(-μ))*trInv(U)
end

"""
    Setting up the Polyakov chain model in a problem that can be used in the solve method

    sutype: U1, SU2, SU3, or SU{N}
    dt: Langevin-step size
    tspan: Simulate in langevin-time 0:dt:tspan
    β: \beta value of the model
        For SU(2) this is the only parameter you need to set,
        for SU(3) you need to give κ and μ
    Nlinks: Number of links in the chain
    NTr: Number of trajectories to run in parallel

    SU(2): S = β/2 Tr[sum_i U_i]
    SU(3): S = Tr[(β + κe^{μ})(sum_i U_i) + (β + κe^{-μ})(sum_i U_i)^{-1}]

"""
function PolyakovChainProblem(sutype::Type{SU{N}},dt,tspan,β,NLinks;NTr=1,κ=0.,μ=0.) where N
    if N == 2 || N == 1
        model = PolyakovChainModel{sutype,typeof(β),typeof(β)}(β,NLinks,0.,0.,0.,0.)
    else
        β₁ = β + κ*exp(μ)
        β₂ = β + κ*exp(-μ)
        model = PolyakovChainModel{sutype,typeof(β),typeof(β₁)}(β,NLinks,κ,μ,β₁,β₂)
    end
    
    u0 = GaugeSQM.IdentityGauges(N,NLinks)
    
    f(V,U,dt,η) = f_chain(V,U,dt,η,model)
    observable(U) = observable_chain(U,model)

    return Problem(u0,dt,tspan,model,NTr,f,observable)
end
