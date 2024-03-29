
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




function f_chain(V::U1AlgebraFields,U::GaugeFields_1D{U1,eType},dt,η,p::PolyakovChainModel{U1,βType}) where {eType,βType}
    for j in 1:U.NV
        M::U1Matrix = U[j]
        for i in vcat(j+1:U.NV, 1:j-1)
            M *= U[i]
        end

        MInv::U1Matrix = inv(U[j])
        for i in vcat(j-1:-1:1, U.NV:-1:j+1)
            MInv *= inv(U[i])
        end

        trT!(view(V.a,:,j),p.β₁*M - p.β₂*MInv)
    end
    
    muladd!(V,-im*dt,η)
end

function f_chain(V::SU2AlgebraFields,U::GaugeFields_1D{SU2,eType},dt,η,cache::Cache,p::PolyakovChainModel{SU2,βType}) where {eType,βType}
    trT!(V,U)
    muladd!(V,dt*im*(p.β/2),η)
end

function jac_chain(V,U::GaugeFields_1D{SU2,eType},dt,cache::Cache,p::PolyakovChainModel{SU2,βType}) where {eType,βType}
    @unpack Uinv = cache
    substitute!(Uinv,U)
    for i in 1:U.NV
        for a in 1:U.NC
            inxi = (i-1)*(U.NC^2-1) + a

            λUi!(Uinv,i,a)
            for j in 1:U.NV

                M::SU2Matrix = Uinv[j]
                for k in vcat(j+1:U.NV, 1:j-1)
                    M *= Uinv[k]
                end
                
                inxj = 1+(j-1)*U.NC
                trT!(view(V,inxi,inxj:inxj+(U.NC^2-1)-1),-dt*(p.β/2)*M)
            end
            Uinv[i] = U[i]
        end
    end         
end

function f_chain(V::SU3AlgebraFields,U::GaugeFields_1D{SU3,eType},dt,η,cache::Cache,p::PolyakovChainModel{SU3,βType}) where {eType,βType}
    
    @unpack Uinv = cache

    # Invert the GaugeFields
    inv!(Uinv,U)
    
    for j in 1:U.NV

        M::SU3Matrix = U[j]
        MInv::SU3Matrix = Uinv[j]
        for i in vcat(j+1:U.NV, 1:j-1)
            M *= U[i]
            MInv = Uinv[i]*MInv
        end


        trT!(view(V.a,:,j),p.β₁*M - p.β₂*MInv)
    end
    
    muladd!(V,im*dt,η)
end

function observable_chain(U::GaugeFields_1D{SU2,eType},model::PolyakovChainModel{SU2,βType}) where {eType,βType}
    return (model.β/2) * tr(U)
    #return unitarity_norm(U)
end

function observable_chain(U::GaugeFields_1D{SU3,eType},model::PolyakovChainModel{SU3,βType}) where {eType,βType}
    #@unpack β₁, β₂ = model
    return tr(U) #+ β₂*trInv(U)
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
    
    f(V,U,dt,η,cache) = f_chain(V,U,dt,η,cache,model)
    j(V,U,dt,cache) = jac_chain(V,U,dt,cache,model)
    observable(U) = observable_chain(U,model)

    return Problem(u0,dt,tspan,model,NTr,f,j,observable)
end

function PolyakovChainProblem(sutype::Type{SU{N}};dt=1e-3,tspan=10.,β=2.0,NLinks=1,NTr=1,κ=0.,μ=0.) where N
    PolyakovChainProblem(sutype,dt,tspan,β,NLinks;NTr=NTr,κ=κ,μ=μ)
end