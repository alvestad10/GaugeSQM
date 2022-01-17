abstract type Solver end

export gEM, gθEM

"""
Geometric Euler-Maruyama scheme

Based in the N=1 approximation of the Magnus solution
"""
struct gEM <: Solver end
Base.show(io::IO, alg::gEM) = print(io, "gEM")

"""
Implicit Geometric Euler-Maruyama scheme

Based in the N=1 approximation of the Magnus solution
where θ refferes to the impliciteness of the scheme
"""
struct gθEM <: Solver
    θ::Float64
end
Base.show(io::IO, alg::gθEM) = print(io, "gθEM: θ=", alg.θ)


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
