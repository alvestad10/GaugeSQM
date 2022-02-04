



function perform_step!(integrator, cache::gEMCache)
    @unpack dt, U, alg, f, η = integrator
    @unpack V, U_tmp, U_tmp2 = cache    
    
    # Populate V with the model function 
    # it inncludes the whole exponential term, also the
    # the noise term, which should be changed when looking into
    # higher-order schemes
    f(V,U,dt,η)
    
    # Generic part of the algorithm not part of the problem    
    expiA!(U_tmp2,V)
    mul!(U_tmp,U_tmp2,U)

    substitute!(integrator.U,U_tmp)
end



function perform_step!(integrator, cache::gθEMCache)
    @unpack dt, U, alg, f, η = integrator
    @unpack U_tmp,U_tmp2, B, D, V, r0 = cache    
    
    # Calculating the initial value
    # for the nlsolver 
    f(D,U,dt,η)
    @. r0 = D.a 

    if (1-alg.θ) != 0.
        #@. B.a = D.a
        #muladd!(B,dt*p*(1-alg.θ),η)
        f(B,U,dt*(1-alg.θ),η)
    else
        @. B.a = η 
    end
    
    g!(F,x) = begin                
        
        # Convert to GaugeField
        @. V.a = x
        
        # Prepare to evaluate drift term
        expiA!(U_tmp2,V)
        mul!(U_tmp,U_tmp2,U)
        
        # Evaluate drift term according to
        # the model
        f(D,U_tmp,dt*alg.θ,B)

        # We want to minimize F ("Lie alebra vector")
        @. F = x - D.a 
    end

    # Perform the non-linear solve
    r = nlsolve(g!, r0, method = :anderson)#, autodiff = :forward)

    # Evaluate the expoenential 
    @. V.a = r.zero
    expiA!(U_tmp2,V)
    mul!(U_tmp,U_tmp2,U)

    substitute!(integrator.U,U_tmp)
end

function adaptive_stepsize_2!(integrator, cache::Cache; ϵ̄ = 1e-3, Kmax=2)

    @unpack U, opts = integrator
    @unpack model = opts

    dSMax = get_dSMax(U, cache, model)
    
    integrator.dt = ϵ̄ * Kmax/dSMax
end

function calcKMax(integrator, cache::Cache)
    @unpack U, opts = integrator
    @unpack model = opts

    return get_dSMax(U, cache, model)
end

function adaptive_stepsize!(integrator, cache::Cache; κ = 5e-4, p=2)

    @unpack U, opts = integrator
    @unpack model = opts

    dSMax = get_dSMax(U, cache, model)
    while !((1/p)*κ <= integrator.dt*dSMax <= p*κ)
        if integrator.dt*dSMax < (1/p)*κ
            integrator.dt = integrator.dt*p
        elseif integrator.dt*dSMax > p*κ
            integrator.dt = integrator.dt*(1/p)
        end
        
        if integrator.dt == Inf
            println(integrator.dt)
            break
        end
    end
end

function get_dSMax(U::GaugeFields_1D{SU3,eType}, cache::Cache, model::PolyakovChainModel{SU3,βType}) where {eType,βType}

    @unpack V = cache

    dSMax = 0.1

    for j in 1:U.NV

        M::SU3Matrix = U[j]
        MInv::SU3Matrix = inv(U[j])
        for i in vcat(j+1:U.NV, 1:j-1)
            M *= U[i]
            MInv = inv(U[i])*MInv
        end


        trT!(view(V.a,:,j),model.β₁*M - model.β₂*MInv)
        dS = maximum(abs.(view(V.a,:,j)))
        if dS > dSMax
            dSMax = dS
        end
    end
    return dSMax
end

function get_dSMax(U::GaugeFields_1D{SU2,eType},cache::Cache,model::PolyakovChainModel{SU2,βType}) where {eType,βType}
    
    @unpack V = cache
    
    trT!(V,U)
    return (model.β/2)*maximum(abs.(V.a))
end