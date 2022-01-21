
abstract type model end

"""
Testing docstring
"""
struct PolyakovChainModel{T <: SUn,βType} <: model
    β::βType
    N::Integer
end


function PolyakovChainProblem(sutype::Type{SU{N}},dt,tspan,β,NLinks;NTr=1) where N
    model = PolyakovChainModel{sutype,typeof(β)}(β,NLinks)
    u0 = GaugeSQM.IdentityGauges(N,NLinks)
    p = im*(β/2)
    return Problem(u0,dt,tspan,model,NTr,p)
end