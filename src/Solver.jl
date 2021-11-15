abstract type Solver end

export gEM, gθEM

"""
Geometric Euler-Maruyama scheme

Based in the N=1 approximation of the Magnus solution
"""
struct gEM <: Solver end


"""
Implicit Geometric Euler-Maruyama scheme

Based in the N=1 approximation of the Magnus solution
where θ refferes to the impliciteness of the scheme
"""
struct gθEM <: Solver
    θ::Float64
end

"""
Geometric Heuns scheme

Based in the N=2 approximation of the Magnus solution
"""
struct gHeuns <: Solver end

"""
Geometric Implicit Midpoint scheme

Based in the N=2 approximation of the Magnus solution
"""
struct gMidpoint <: Solver end


"""
RK
"""
struct gRK <: Solver end
