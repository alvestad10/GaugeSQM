using StaticArrays

export SU2Matrix

#############################
abstract type SUNMatrix end

import Base: *, +
import Base: copy, exp
import LinearAlgebra: tr, adjoint

const SMatrix2 = SMatrix{2,2,Complex{Real}}
const SMatrix3 = SMatrix{3,3,Complex{Real}}

mutable struct SU2Matrix <: SUNMatrix
    M::SMatrix2
end

SU2Matrix(a,b,c,d) = SU2Matrix(SMatrix2(a,b,c,d))
*(A::T,B::T) where {T <: SUNMatrix} = T(A.M*B.M)
*(A::Number,B::T) where {T <: SUNMatrix} = T(A*B.M)
*(A::T,B::Number) where {T <: SUNMatrix} = T(A.M*B)
+(A::T,B::T) where {T <: SUNMatrix} = T(A.M + B.M)
+(A::Number,B::T) where {T <: SUNMatrix} = T(A + B.M)
+(A::T,B::Number) where {T <: SUNMatrix} = T(A.M + B)
copy(A::T) where {T <: SUNMatrix} = T(A.M)
exp(A::T) where {T <: SUNMatrix} = T(exp(A.M))
tr(A::T) where {T <: SUNMatrix} = tr(A.M)
adjoint(A::T) where {T <: SUNMatrix} = T(adjoint(A.M))

#const SU2Matrix = SMatrix{2,2,Complex{Float64}}
#const su2CAlgebra = SVector{6,Float64}

II = SU2Matrix(1,0,0,1)

hat_complex(ω) = SU2Matrix(ω[3] + im*ω[6],
                          (ω[1] + im*ω[4]) + im*(ω[2] + im*ω[5]),
                          (ω[1] + im*ω[4]) - im*(ω[2] + im*ω[5]), 
                          -ω[3] - im*ω[6])

hat(ω) = SU2Matrix(ω[3], ω[1] + im*ω[2], ω[1] - im*ω[2],  -ω[3])

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
    if i == 1 return SU2Matrix(0., 1., 1., 0.) end
    if i == 2 return SU2Matrix(0., im, -im, 0.) end
    if i == 3 return SU2Matrix(1., 0., 0., -1.) end
end

"""
    TODO:
"""
function trT(X::T,i::Integer) where {T <: SUNMatrix}
    if i == 1 return X.M[2,1] + X.M[1,2]      end #SU2Matrix(0., 1., 1., 0.) end
    if i == 2 return im*(X.M[1,2] - X.M[2,1]) end #SU2Matrix(0., im, -im, 0.) end
    if i == 3 return X.M[1,1] - X.M[2,2]      end #SU2Matrix(1., 0., 0., -1.) end
end

"""
    TODO:
"""
function trT(X::T) where {T <: SUNMatrix}
    return [X.M[2,1] + X.M[1,2],
        im*(X.M[1,2] - X.M[2,1]),
            X.M[1,1] - X.M[2,2]]
end