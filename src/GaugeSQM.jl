module GaugeSQM

using LinearAlgebra
using NLsolve
using Random
import Parameters: @unpack

include("Manifolds.jl")
include("GaugeFields.jl")

include("Solver.jl")
include("plotHelperFuncitons.jl")
include("model.jl")
include("Integrator.jl")



export run_sim2



function run_sim(p::Setup,alg::Solver)

    Random.seed!(123)

    @unpack dt, tspan, model, NTr = p
    β = model.β

    nrSaves = length(0:dt:tspan)

    sol =            zeros(ComplexF64,NTr,nrSaves+1)

    Threads.@threads for tread in 1:NTr
    #for tread in 1:NTr
        U = copy(p.u0)
        U_tmp = similar(U)
        B = similar(U)

        V = LieAlgebraFields(ComplexF64,U.NC,U.NV)
        sol[tread,1] = (β/2) * tr(U)

        for (i,t) in enumerate(0:dt:tspan)
            R = im*(β / 2)
            
            η = sqrt(2*dt) .* randn(U.NC^2-1,U.NV)

            if typeof(alg) == gEM
                trT!(V,U)
                muladd!(V,dt*R,η)
                expiA!(B,V)
                mul!(U_tmp,B,U)

                substitute!(U,U_tmp)
            end
            sol[tread,i+1] = (β/2) * tr(U)
            

            #if isnan(sol[i]) return manifold_check[1:i-1], sol[1:i-1] end
    end
    return sol
end








end # module
