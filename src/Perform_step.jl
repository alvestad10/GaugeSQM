



function perform_step!(integrator, cache::gEMCache)
    @unpack dt, U, alg, f, η = integrator
    @unpack V, U_tmp, U_tmp2 = cache    
    
    f(V,U,dt,η)
    
    # Related to the problem at the moment
    #trT!(V,U)
    
    # Generic part of the algorithm not part of the problem    
    #muladd!(V,dt*p,η)
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