
abstract type model end

"""
Testing docstring
"""
struct PolyakovChainModel{T <: SUn,βType} <: model
    β::βType
    N::Integer
end



function PolyakovChainProblem(sutype,dt,tspan,β,NLinks;NTr=1)
    model = PolyakovChainModel{sutype,typeof(β)}(β,NLinks)
    u0 = GaugeSQM.IdentityGauges(2,NLinks)
    return Setup(u0,dt,tspan,model,NTr)
end