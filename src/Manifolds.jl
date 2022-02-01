using StaticArrays

export SU2Element

#############################
abstract type Manifold end

abstract type UNElement end
abstract type SUNElement <: Manifold end

import Base: *, +
import Base: copy, exp
import LinearAlgebra: tr, adjoint, det, mul!

const SMatrix2 = SMatrix{2,2,Complex{Real}}
const SMatrix3 = SMatrix{3,3,Complex{Real}}

mutable struct SU2Element <: SUNElement
    M::SMatrix2
end

SU2Element(a,b,c,d) = SU2Element(SMatrix2(a,b,c,d))
SU2Element(a,b,c) = _expi([a,b,c])

mul!(A::T,B::T) where {T <: SUNElement} = A.M = A.M*B.M; nothing
*(A::T,B::T) where {T <: SUNElement} = T(A.M*B.M)
*(A::Number,B::T) where {T <: SUNElement} = T(A*B.M)
*(A::T,B::Number) where {T <: SUNElement} = T(A.M*B)
+(A::T,B::T) where {T <: SUNElement} = T(A.M + B.M)
+(A::Number,B::T) where {T <: SUNElement} = T(A + B.M)
+(A::T,B::Number) where {T <: SUNElement} = T(A.M + B)
copy(A::T) where {T <: SUNElement} = T(A.M)
exp(A::T) where {T <: SUNElement} = T(exp(A.M))
tr(A::T) where {T <: SUNElement} = tr(A.M)
adjoint(A::T) where {T <: SUNElement} = T(adjoint(A.M))
det(A::T) where {T <: SUNElement} = det(A.M)

#const SU2Matrix = SMatrix{2,2,Complex{Float64}}
#const su2CAlgebra = SVector{6,Float64}

II = SU2Element(1,0,0,1)

hat_complex(ω) = SU2element(ω[3] + im*ω[6],
                          (ω[1] + im*ω[4]) + im*(ω[2] + im*ω[5]),
                          (ω[1] + im*ω[4]) - im*(ω[2] + im*ω[5]), 
                          -ω[3] - im*ω[6])

hat(ω) = SU2Element(ω[3], ω[1] + im*ω[2], ω[1] - im*ω[2],  -ω[3])

function _expi(ω) 
    normA = sqrt(sum(ω.^2))
    B = hat(ω)
    return cos(normA)*II + im*sin(normA)*B*inv(normA)
end

function _expi_complex(ω) 
    normA = sqrt((ω[1]+im*ω[4])^2 + (ω[2]+im*ω[5])^2 + (ω[3]+im*ω[6])^2)
    B = hat_complex(ω)
    return cos(normA)*II + im*sin(normA)*B*inv(normA)
end

"""
    The representation of the su2 Lie Algebra
    which is the Pauli matrices

    Input
    i: Pauli matrix to get
"""
function T(i::Integer)
    if i == 1 return SU2Element(0., 1., 1., 0.) end
    if i == 2 return SU2Element(0., im, -im, 0.) end
    if i == 3 return SU2Element(1., 0., 0., -1.) end
end

"""
    TODO:
"""
function trT(X::T,i::Integer) where {T <: SUNElement}
    if i == 1 return X.M[2,1] + X.M[1,2]      end #SU2Matrix(0., 1., 1., 0.) end
    if i == 2 return im*(X.M[1,2] - X.M[2,1]) end #SU2Matrix(0., im, -im, 0.) end
    if i == 3 return X.M[1,1] - X.M[2,2]      end #SU2Matrix(1., 0., 0., -1.) end
end

"""
    TODO:
"""
function trT(X::T) where {T <: SUNElement}
    return [X.M[2,1] + X.M[1,2],
        im*(X.M[1,2] - X.M[2,1]),
            X.M[1,1] - X.M[2,2]]
end