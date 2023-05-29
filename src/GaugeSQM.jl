module GaugeSQM

using LinearAlgebra
using Statistics
using NLsolve
using Random
import Parameters: @unpack, @with_kw

include("Manifolds.jl")
include("GaugeFields.jl")

include("Solver.jl")
include("plotHelperFuncitons.jl")
include("Model.jl")
include("Regulators.jl")
include("Integrator.jl")
include("Perform_step.jl")


export solve



function solve(opts::Problem,alg::Solver,regs::Regulators; adaptive=false, ad_κ=5e-4, ad_p=2, saveat=0.001, cb=(integrator) -> nothing)

    @unpack U0, dt, tspan, model, NTr = opts
    
    Random.seed!(123)
    
    if saveat > dt
        nrSaves = length(0:saveat:tspan)
    else
        @warn "dt is larger that saveat, using dt instead of saveat for the save at interval"
        nrSaves = length(0:dt:tspan)
        saveat = dt
    end
    
    sol = zeros(ComplexF64,NTr,nrSaves+1)
    


    #Threads.@threads for tread in 1:NTr
    for tread in 1:NTr
        
        integrator = Integrator(U0,dt,opts,alg,regs)
        algCache = get_cache(integrator,alg)
        GCCache = get_cache(integrator,regs.GC)

        old_f = integrator.f
        
        
        if typeof(regs.DS) == DynamicStabilization
            DSCache = get_cache(integrator,regs.DS)
            
            function ff(V,U,dt,η,cache)
                @unpack M = DSCache
                old_f(V,U,dt,η,cache)
                
                get_DynamicStabilization!(M,integrator,DSCache,regs)

                @. V.a += im*dt*M.a
            end
            
            integrator = Integrator(integrator,ff)
        end 

        
        sol[tread,1] = opts.observable(integrator.U)

        Kmaxes = Float64[1.]
        Kmax = 1.
        dt0 = integrator.dt

        save_this_round = false
        prev_save = 0.
        i = 2
        t = 0.
        while i <= nrSaves + 1
        #for (i,t) in enumerate(0:dt:tspan)
            
            if adaptive
                if t < 1
                    append!(Kmaxes,calcKMax(integrator,algCache))
                    Kmax = mean(Kmaxes)
                end
                
                adaptive_stepsize_2!(integrator, algCache; ϵ̄=dt0, Kmax=Kmax)
                #adaptive_stepsize!(integrator, algCache; κ=ad_κ, p=ad_p)

                if prev_save + saveat <= t + integrator.dt + 5e-6
                    integrator.dt = prev_save + saveat - t
                    save_this_round = true
                    prev_save += saveat
                end
            else
                if prev_save + saveat <= t + integrator.dt
                    save_this_round = true
                    prev_save += saveat
                end
            end

            num_of_basis = integrator.U.NC > 1 ? integrator.U.NC^2-1 : 1
            integrator.η = sqrt(2*integrator.dt) .* randn(num_of_basis,integrator.U.NV)
            if !perform_step!(integrator,algCache)
                @warn "NLSolver failes to converge at t=$t, use smaller dt"
                break
            end
            
            GaugeCoolingUpdate!(integrator,GCCache)
            
            if save_this_round
                sol[tread,i] = opts.observable(integrator.U)
                save_this_round = false
                i += 1
            end

            cb(integrator)

            t += integrator.dt
        end
    end
    return sol
end








end # module
