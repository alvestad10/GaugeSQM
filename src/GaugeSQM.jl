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
    tspan::Int
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


    

    dt = p.dt
    tspan = p.tspan
    β = p.model.β

    NTr = 2
    nrSaves = length(0:dt:tspan)

    #manifold_check = zeros(ComplexF64,NTr,nrSaves)
    sol =            zeros(ComplexF64,NTr,nrSaves)

    Threads.@threads for tread in 1:NTr

    U = copy(p.u0)

    A = zeros(ComplexF64,3)
    B = zeros(ComplexF64,3)
    D = zeros(ComplexF64,3)
    u0 = zeros(Float64,6)

    for (i,t) in enumerate(0:dt:tspan)

        R = im*(β / 2)
        η = sqrt(2*dt) .* randn(3)

        if typeof(alg) == gEM
            #A =  im *( ( dt * R * trT(U,1) + η[1])*T(1) + 
            #           ( dt * R * trT(U,2) + η[2])*T(2) + 
            #           ( dt * R * trT(U,3) + η[3])*T(3) )

            A .= (dt*R*trT(U) .+ η)
            
            U = _expi(A)*U

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

            r = nlsolve(g!, u0, method = :newton,autodiff = :forward)

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

        elseif typeof(alg) == gRK
            A .= (dt*R*trT(U) .+ η)
            A .= ((dt/2)*(1+2*dt/6)*R*(trT(_expi(A)*U) + trT(U)) .+ η)
            U = _expi(A)*U
        end

        #manifold_check[i] = tr(U*adjoint(U))/2 - 1
        sol[tread,i] = (β/2) * tr(U)

        #if isnan(sol[i]) return manifold_check[1:i-1], sol[1:i-1] end
    end 
    end
    #return manifold_check, sol
    return sol
end





end # module
