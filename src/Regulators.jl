export Regulators
export GaugeCooling, GaugeCoolingAdaptive, NoGaugeCooling, NoDynamicStabilization
export GaugeCoolingUpdate!, adaptiveUpdateGCParameter

########################
abstract type AbstractGaugeCooling end
abstract type AbstractDynamicStabilization end

struct Regulators{GCType <: AbstractGaugeCooling, 
                  DSType <: AbstractDynamicStabilization}
    GC::GCType
    DS::DSType
end


## Dynamic Stabilization
struct NoDynamicStabilization <: AbstractDynamicStabilization end


## GaugeCooling
struct NoGaugeCooling <: AbstractGaugeCooling end
struct NoGaugeCoolingCache <: Cache end

get_cache(integrator,regs::Regulators{NoGaugeCooling,NoDynamicStabilization}) = NoGaugeCoolingCache()

struct GaugeCooling <: AbstractGaugeCooling
    α::Float64
    N::Float64
end

struct GaugeCoolingCache{T <: SUn,eType,lType} <: Cache
    U2::GaugeFields{T,eType}
    U3::GaugeFields{T,eType}
    U4::GaugeFields{T,eType}
    V1::LieAlgebraFields{T,lType}
    V2::LieAlgebraFields{T,lType}
end

function get_cache(integrator,regs::Regulators{GaugeCooling,DSType}) where {DSType}
    @unpack U = integrator

    U2 = similar(U)
    U3 = similar(U)
    U4 = similar(U)
    V1 = LieAlgebraFields(ComplexF64,U.NC,U.NV)
    V2 = LieAlgebraFields(ComplexF64,U.NC,U.NV)

    GaugeCoolingCache(U2,U3,U4,V1,V2)
end

struct GaugeCoolingAdaptive <: AbstractGaugeCooling
    α_eff::Float64
    D_0::Float64
    N::Integer
end

function get_cache(integrator,regs::Regulators{GaugeCoolingAdaptive,DSType}) where {DSType}
    @unpack U = integrator

    U2 = similar(U)
    U3 = similar(U)
    U4 = similar(U)
    V1 = LieAlgebraFields(ComplexF64,U.NC,U.NV)
    V2 = LieAlgebraFields(ComplexF64,U.NC,U.NV)

    GaugeCoolingCache(U2,U3,U4,V1,V2)
end

#"Polyakov loop"
#function GaugeCoolingUpdate!(integrator,cache::Cache)
#    @unpack U, Regulators = integrator
#
#    GaugeCoolingUpdate!(U,cache,Regulators.GC)
#end



"GC drift for SU(2)"
function GC_drift!(V::LieAlgebraFields{SU2,aType},V2::LieAlgebraFields{SU2,aType},U::GaugeFields{SU2,eType}) where {eType,aType}
    for j in 1:U.NV
        jm1 = j-1 == 0 ? U.NV : j-1
        trT!(view(V.a,:,j),U[j]*adjoint(U[j]) - adjoint(U[jm1])*U[jm1])
        
        jp1 = j+1 == U.NV + 1  ? 1 : j+1
        trT!(view(V2.a,:,j),U[jp1]*adjoint(U[jp1]) - adjoint(U[j])*U[j])
    end
end

"GC drift for SU(3)"
function GC_drift!(V::LieAlgebraFields{SU3,aType},V2::LieAlgebraFields{SU3,aType},U::GaugeFields{SU3,eType}) where {eType,aType}
    for j in 1:U.NV
        jm1 = j-1 == 0 ? U.NV : j-1
        trT!(view(V.a,:,j),U[j]*adjoint(U[j]) - adjoint(U[jm1])*U[jm1] - inv(adjoint(U[j]))*inv(U[j]) + inv(adjoint(U[jm1]))*inv(U[jm1]))
        
        jp1 = j+1 == U.NV + 1  ? 1 : j+1
        trT!(view(V2.a,:,j),U[jp1]*adjoint(U[jp1]) - adjoint(U[j])*U[j] - inv(adjoint(U[jp1]))*inv(U[jp1]) + inv(adjoint(U[j]))*inv(U[j]))
    end
end

"Adaptive Gauge Cooling"
function adaptiveUpdateGCParameter(GC::GaugeCoolingAdaptive,cache::GaugeCoolingCache,integrator)
    @unpack V1, V2 = cache 
    GC_drift!(V1,V2,integrator.U)

    G = 2*sum((x) -> sum(abs.(x)),V1) / integrator.U.NV
    return GC.α_eff / (G + GC.D_0)
end


function GaugeCoolingUpdate!(integrator,cache::NoGaugeCoolingCache) end

function GaugeCoolingUpdate!(integrator,cache::GaugeCoolingCache)
    
    @unpack U, dt, Regulators = integrator
    @unpack GC = Regulators
    
    for i in 1:GC.N
        if typeof(GC) == GaugeCoolingAdaptive
            α = dt*adaptiveUpdateGCParameter(GC,cache,integrator)
        elseif typeof(GC) == GaugeCooling
            α = dt*GC.α
        end
        GaugeCoolingUpdate!(U,cache; α=α)
    end
    substitute!(integrator.U,U)
end


function GaugeCoolingUpdate!(U::GaugeFields{SU{N},eType},cache::GaugeCoolingCache; α = 0.001) where {N,eType}
    
    @unpack U2, U3, U4, V1, V2 = cache

    GC_drift!(V1,V2,U)
    @. V1.a = -2*α*V1.a
    @. V2.a = 2*α*V2.a


    expA!(U2,V1)
    expA!(U3,V2)
    mul!(U4,U2,U)
    mul!(U,U4,U3)
end