abstract type Solver end
abstract type Cache end

export gEM, gθEM

"""
Geometric Euler-Maruyama scheme

Based in the N=1 approximation of the Magnus solution
"""
struct gEM <: Solver end
Base.show(io::IO, alg::gEM) = print(io, "gEM")

struct gEMCach{T <: SUn,eType,lType} <: Cache
    U_tmp::GaugeFields{T,eType}
    U_tmp2::GaugeFields{T,eType}
    V::LieAlgebraFields{T,lType}
end

function get_cach(integrator,alg::gEM)
    @unpack U = integrator

    U_tmp = similar(U)
    U_tmp2 = similar(U)
    V = LieAlgebraFields(ComplexF64,U.NC,U.NV)

    gEMCach(U_tmp,U_tmp2,V)    
end

"""
Implicit Geometric Euler-Maruyama scheme

Based in the N=1 approximation of the Magnus solution
where θ refferes to the impliciteness of the scheme
"""
struct gθEM <: Solver
    θ::Float64
end
Base.show(io::IO, alg::gθEM) = print(io, "gθEM: θ=", alg.θ)

struct gθEMCach{T <: SUn,eType,lType,r0Type} <: Cache
    U_tmp::GaugeFields{T,eType}
    U_tmp2::GaugeFields{T,eType}
    B::LieAlgebraFields{T,lType}
    D::LieAlgebraFields{T,lType}
    V::LieAlgebraFields{T,lType}
    r0::r0Type
end

function get_cach(integrator,alg::gθEM)
    @unpack U = integrator

    U_tmp = similar(U)
    U_tmp2 = similar(U)
    V = LieAlgebraFields(ComplexF64,U.NC,U.NV)
    B = LieAlgebraFields(ComplexF64,U.NC,U.NV)
    D = LieAlgebraFields(ComplexF64,U.NC,U.NV)
    
    r0 = D.a
    #r0 = zeros(Float64,2*(U.NC^2-1),U.NV)

    gθEMCach(U_tmp,U_tmp2,V,B,D,r0)    
end

"""
Implicit Geometric Euler-Maruyama scheme

Based in the N=1 approximation of the Magnus solution
where θ refferes to the impliciteness of the scheme
"""
struct gθEM2 <: Solver
    θ::Float64
end
Base.show(io::IO, alg::gθEM2) = print(io, "gθEM2: θ=", alg.θ)


"""
Geometric Heuns scheme

Based in the N=2 approximation of the Magnus solution
"""
struct gHeuns <: Solver end
Base.show(io::IO, alg::gHeuns) = print(io, "gHeuns")

"""
Geometric Implicit Midpoint scheme

Based in the N=2 approximation of the Magnus solution
"""
struct gMidpoint <: Solver end
Base.show(io::IO, alg::gMidpoint) = print(io, "gMidpoint")

"""
RK
"""
struct gRK <: Solver end






#################
# HELP FUNCTIONS
#################
@with_kw struct AdaptiveStepsize
    κ = 1e-3
    p = 2
end

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