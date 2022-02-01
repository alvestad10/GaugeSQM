module GaugeSQM

using LinearAlgebra
using NLsolve
using Random
import Parameters: @unpack, @with_kw

include("Manifolds.jl")
include("GaugeFields.jl")

include("Solver.jl")
include("plotHelperFuncitons.jl")
include("Model.jl")
include("GaugeCooling.jl")
include("Integrator.jl")
include("Perform_step.jl")


export solve



function solve(opts::Problem,alg::Solver,regs::Regulators)

    @unpack U0, dt, tspan, model, NTr = opts
    
    Random.seed!(123)
    
    
    nrSaves = length(0:dt:tspan)
    sol = zeros(ComplexF64,NTr,nrSaves+1)
    
    Threads.@threads for tread in 1:NTr
    #for tread in 1:NTr
        
        integrator = Integrator(U0,dt,opts,alg,regs)
        algCache = get_cache(integrator,alg)
        regsCache = get_cache(integrator,regs)


        sol[tread,1] = (model.β/2) * tr(integrator.U)

        for (i,t) in enumerate(0:dt:tspan)
            
            num_of_basis = integrator.U.NC > 1 ? integrator.U.NC^2-1 : 1
            integrator.η = sqrt(2*dt) .* randn(num_of_basis,integrator.U.NV)
            perform_step!(integrator,algCache)
            
            GaugeCoolingUpdate!(integrator,regsCache)
            
            sol[tread,i+1] = (model.β/2) * tr(integrator.U)
            
        end
    end
    return sol
end








end # module
