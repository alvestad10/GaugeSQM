module GaugeSQM

using LinearAlgebra
using NLsolve
using Random

include("Manifolds.jl")
include("Solver.jl")
include("plotHelperFuncitons.jl")

struct PolyakovChainModel
    β::ComplexF64
end
struct Setup
    u0::SU2Matrix
    dt::Float64
    tspan::Float64
    model::PolyakovChainModel
end



struct Integrator{S<:Solver}
    opts::Setup
    alg::S
end

#function get_PolyakovLoopModel(β::Complex{Float64})
#    PolyakovChainModel(U::SU2Matrix,i::Int64) = im*(β / 2) * tr(T(i)*U)
#    return PolyakovChainModel
#end

function adaptive_stepsize(U,dt; κ = 1e-3, p=2)
    
    _dt = dt

    A = abs.((β / 2)*[ tr(T1*U), tr(T2*U), tr(T3*U)])
    dSMax = maximum(A)
    while !((1/p)*κ <= dt*dSMax <= p*κ)
        if dt*dSMax < (1/p)*κ
            dt = dt*p
        elseif dt*dSMax > p*κ
            dt = dt*(1/p)
        end
        
        if dt == Inf
            println(dt)
            return _dt
        end
    end
    return dt
end

function PolyakovChainProblem(u0,dt,tspan,β)
    model = PolyakovChainModel(β)
    return Setup(u0,dt,tspan,model)
end

function run_sim(p::Setup,alg::Solver)

    Random.seed!(123)

    dt = p.dt
    tspan = p.tspan
    β = p.model.β

    NTr = 1
    nrSaves = length(0:dt:tspan)

    #manifold_check = zeros(ComplexF64,NTr,nrSaves)
    sol =            zeros(ComplexF64,NTr,nrSaves)
    sol2 =            zeros(ComplexF64,NTr,nrSaves)
    x =            zeros(ComplexF64,NTr,nrSaves)
    y =            zeros(ComplexF64,NTr,nrSaves)
    z =            zeros(ComplexF64,NTr,nrSaves)

    II = Matrix(Diagonal(ones(2)))

    Threads.@threads for tread in 1:NTr

    U = copy(p.u0)
    
    A = zeros(ComplexF64,3)
    B = zeros(ComplexF64,3)
    D = zeros(ComplexF64,3)
    u0 = zeros(Float64,6)

    for (i,t) in enumerate(0:dt:tspan)

        R = im*(β / 2)
        #R = im*(β*conj(β) / 2)
        
        #η = zeros(3)#sqrt(2*dt) .* randn(3)
        η = sqrt(2*dt) .* randn(3)
        #η = sqrt(2*dt*conj(β)) .* randn(3)

        if typeof(alg) == gEM
            #A =  im *( ( dt * R * trT(U,1) + η[1])*T(1) + 
            #           ( dt * R * trT(U,2) + η[2])*T(2) + 
            #           ( dt * R * trT(U,3) + η[3])*T(3) )

            A .= (dt*R*trT(U) .+ η)

            U = _expi(A)*U
            
            # using Cay instad of Exp
            #A .= 0.5*(dt*R*trT(U) .+ η)

            #U = SU2Matrix(inv(II - hat(A).M))*SU2Matrix((II + hat(A).M))*U

            #U = exp(A)*U
        elseif typeof(alg) == gθEM
            A .= dt*R*trT(U) .+ η

            u0[1:3] .= real(A)
            u0[4:6] .= imag(A)

            B .= dt * R .* (1-alg.θ) .* trT(U) .+  η

            g!(F,x) = begin                
                D = ( dt * R .* alg.θ .* trT(_expi(x)*U) .+ B)

                F[1:div(end,2)] .= x[1:div(end,2)] .- real(D)
                F[div(end,2)+1:end] .= x[div(end,2)+1:end] .- imag(D)
            end

            r = nlsolve(g!, u0, method = :newton, autodiff = :forward)
            U = _expi_complex(r.zero)*U            
        
        elseif typeof(alg) == gHeuns
            A .= (dt*R*trT(U) .+ η)
            A .= (dt*R*trT(_expi(A)*U) .+ η)
            U = _expi(A)*U
        elseif typeof(alg) == gMidpoint
            A .= dt*R*trT(U) .+ η

            u0[1:3] .= real(A)
            u0[4:6] .= imag(A)


            _g!(F,x) = begin
                X = _expi(x)*U                
                
                D .= dt * R .* trT(X) .+ η

                F[1:div(end,2)] .= x[1:div(end,2)] .- real(D)
                F[div(end,2)+1:end] .= x[div(end,2)+1:end] .- imag(D)
            end

            r = nlsolve(_g!, u0, method = :newton,autodiff = :forward)

            U = _expi_complex(r.zero)*U

        elseif typeof(alg) == gKronfeld
            A .= (dt*R*trT(U) .+ η)
            A .= ((dt/2)*(1+2*dt/6)*R*(trT(_expi(A)*U) + trT(U)) .+ η)
            U = _expi(A)*U
        elseif typeof(alg) == gRK15

            ΔZ = 0.5*dt*(η .+ randn(3)*sqrt(2*dt/3))

            A(X) = dt*R*trT(X)

            H1 = 0
            H2 = (3/4)*A(H1)*dt + (3/2)*ΔZ/dt
            H3 = 0
            
            H̃1 = 0
            H̃2 = (1/9)*A(H1)*dt + (1/3)*sqrt(dt)
            H̃3 = (5/9)*A(H1)*dt + (1/3)*A(H2)*dt - (1/3)*sqrt(2*dt)+sqrt(2*dt)
            H̃4 = A(H1)*dt + (1/3)*A(H2)*dt + A(H3)*dt + sqrt(2*dt) - sqrt(2*dt) + sqrt(2*dt)
            
            Ω1 .= ((1/3)*A(H1) + (2/3)*A(H2))*dt +
                  ((13/4) - (9/4) - (9/4) + (9/4))*η
            _A .= ((dt/2)*(1+2*dt/6)*R*(trT(_expi(A)*U) + trT(U)) .+ η)
            U = _expi(A)*U
        end

        #manifold_check[i] = tr(U*adjoint(U))/2 - 1
        sol[tread,i] = (β/2) * tr(U)
        
        x[tread,i] = U.M[1,1]
        y[tread,i] = U.M[2,1]

        #if isnan(sol[i]) return manifold_check[1:i-1], sol[1:i-1] end
    end 
    end
    #return manifold_check, sol
    return sol, x, y
end





end # module
