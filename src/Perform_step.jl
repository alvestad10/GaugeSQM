



function perform_step!(integrator, cache::gEMCach)
    @unpack dt, U, alg, p, η = integrator
    @unpack V, U_tmp, U_tmp2 = cache    
    # Related to the problem at the moment
    trT!(V,U)
    
    # Generic part of the algorithm not part of the problem    
    muladd!(V,dt*p,η)
    expiA!(U_tmp2,V)
    mul!(U_tmp,U_tmp2,U)

    substitute!(integrator.U,U_tmp)
end



function perform_step!(integrator, cache::gθEMCach)
    @unpack dt, U, alg, p, η = integrator
    @unpack U_tmp,U_tmp2, B, D, V, r0 = cache    
    # Related to the problem at the moment
    trT!(B,U)
    @. D.a = B.a
    
    # Generic part of the algorithm not part of the problem    
    muladd!(D,dt*p,η)
    @. r0 = D.a 
    #@. r0[1:div(end,2),:] = real(D.a)
    #@. r0[div(end,2)+1:end,:] = imag(D.a)
    
    muladd!(B,dt*p*(1-alg.θ),η)
    
    
    g!(F,x) = begin                
        #D = ( dt * R .* alg.θ .* trT(_expi(x)*U) .+ B)
        
        @. V.a = x
        #V.a .= (@view x[1:div(end,2),:]) + im*(@view x[div(end,2)+1:end,:])
        
        expiA!(U_tmp2,V)
        mul!(U_tmp,U_tmp2,U)
        
        trT!(D,U_tmp)
        muladd!(D,dt*p*alg.θ,B)

        @. F = x - D.a 
        #@. F[1:div(end,2),:] = x[1:div(end,2),:] - real(D.a)
        #@. F[div(end,2)+1:end,:] = x[div(end,2)+1:end,:] - imag(D.a)
    end

    
    r = nlsolve(g!, r0, method = :anderson)#, autodiff = :forward)

    @. V.a = r.zero
    #@. V.a = (@view r.zero[1:div(end,2),:]) + im*(@view r.zero[div(end,2)+1:end,:])

    expiA!(U_tmp2,V)
    mul!(U_tmp,U_tmp2,U)

    substitute!(integrator.U,U_tmp)
end