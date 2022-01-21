using StaticArrays

export U1, SU2, SU3
export GaugeFields
export SUMatrix, SU2Matrix

const SUMatrix{N,L} = SMatrix{N,N,ComplexF64,L}
const SU2Matrix = SUMatrix{2,4}

### Setting up the manifold ####
abstract type SUn end
abstract type SU{N} <: SUn end


const U1 = SU{1}
const SU2 = SU{2}
const SU3 = SU{3}


###### 1D Gauge field #####
struct GaugeFields{T <: SUn, eType} 
    g::Array{eType}
    NC::Int64
    NV::Int64

    function GaugeFields(NC,NV) 
        sutype = SU{NC}
        g = zeros(SUMatrix{NC,NC^2},NV)
        eType = eltype(g)
        return new{sutype,eType}(g,NC,NV)
    end
end

const U1GaugeFields  = GaugeFields{U1}
const SU2GaugeFields  = GaugeFields{SU2}
const SU3GaugeFields  = GaugeFields{SU3}
const SUNGaugeFields{N}  = GaugeFields{SU{N}}


function Base.setindex!(x::GaugeFields{T,eType},v,i) where {T <: SUn, eType} 
    x.g[i] = v
end

function Base.getindex(x::GaugeFields,i)
    return x.g[i]
end

function Base.similar(x::GaugeFields) 
    return GaugeFields(x.NC,x.NV)
end

function Base.copy(x::GaugeFields)
    S = similar(x)
    S.g .= copy(x.g)
    return S
end

function substitute!(a::GaugeFields{SU{NC},eType},b::GaugeFields{SU{NC},eType}) where {NC,eType}
    a.g .= b.g
    return
end

function clear!(a::GaugeFields{T,eType}) where {T <: SUn, eType}
    @. a.g = zero(eType)
    return 
end


identity(::SUMatrix{N,L}) where {N,L} = SUMatrix{N,L}(Diagonal(ones(N)))
identity(::Type{SUMatrix{N,L}}) where {N,L} = SUMatrix{N,L}(Diagonal(ones(N)))
identity(::GaugeFields{SU{N},eType}) where {N,eType} = SUMatrix{N,N^2}(Diagonal(ones(N)))

function identity!(a::GaugeFields{T,eType}) where {T <: SUn, eType}
    @. a.g = identity(eType)
    return 
end

function IdentityGauges(NC,NV)
    U = GaugeFields(NC,NV)
    identity!(U)
    return U
end




function RandomGauges(NC,NV)
    U = GaugeFields(NC,NV)
    A = LieAlgebraFields(Float64,NC,NV)
    
    gauss_distribution_lie!(A)
    expiA!(U,A)
            
    return U
end

function LinearAlgebra.mul!(c::GaugeFields{T,eType},a::GaugeFields{T,eType}, b::GaugeFields{T,eType}) where {T <: SUn, eType}
    @. c.g = a.g * b.g
    return
end










#####################
# Lie Algebra Fields
#####################

export LieAlgebraFields

abstract type AbstractLieAlgebraFields end

struct LieAlgebraFields{T <: SUn,fType <: Number}
    a::Array{fType,2}
    NC::Int64
    NV::Int64
    NumOfBasis::Int64
    #generators::Generator

    function LieAlgebraFields(fType,NC,NV)
        sutype = SU{NC}
        NumofBasis = NC^2-1
        #generators = Generator(NC)
        return new{sutype,fType}(zeros(fType,NumofBasis,NV),NC,NV,NumofBasis)
    end
end

const U1AlgebraFields{T} = LieAlgebraFields{U1,T}
const SU2AlgebraFields{T} = LieAlgebraFields{SU2,T}
const SU3AlgebraFields{T} = LieAlgebraFields{SU3,T}

function Base.setindex!(x::LieAlgebraFields,v,i...) 
    x.a[i...] = v
end

function Base.getindex(x::LieAlgebraFields,i...)
    return x.a[i...]
end

function Base.similar(x::LieAlgebraFields)
    return LieAlgebraFields(eltype(x.a),x.NC,x.NV)
end

#function GaugeFields.substitute!(x::LieAlgebraFields,pwork)
#    n1,n2 = size(x.a)
    #println(size(pwork))
#    x.a[:,:] = reshape(pwork,(n1,n2))
#end

function Base.display(a::LieAlgebraFields)
    n1,n2 = size(a.a)
    for i2=1:n2
        for i1=1:n1
            println("$i1 $i2 ",a.a[i1,i2])# = α*a.a[i1,i2,i3,i4,i5]
        end
    end
    return 
end

function clear!(a::LieAlgebraFields)
    @. a.a = 0
    return 
end

function LinearAlgebra.mul!(c::LieAlgebraFields,α::Number,a::LieAlgebraFields)
    @. c.a = α*a.a
    return
end

function Base.:*(x::T,y::T) where T <: LieAlgebraFields
    return sum(x.a .* y.a)
end

function add!(a::LieAlgebraFields,α,b::LieAlgebraFields)
    @. a.a += α*b.a
    return
end

function muladd!(a::LieAlgebraFields,α,η)
    @. a.a = α*a.a + η
    return
end

function gauss_distribution_lie!(p::LieAlgebraFields)
    @inbounds @simd for i=1:p.NV
        p.a[1:div(p.NumOfBasis,2),i] = convert(typeof(p.a),randn(div(p.NumOfBasis,2),1))
    end
    return
end


const pi23 = 2pi/3
const tinyvalue =1e-100

function expiA!(v::SU2GaugeFields, u::SU2AlgebraFields)   
    NV=u.NV
    @inbounds @simd for i=1:NV
        
        #u1 = (u[1,i] + im*u[4,i])/2
        #u2 = (u[2,i] + im*u[5,i])/2
        #u3 = (u[3,i] + im*u[6,i])/2

        u1 = u[1,i]#/2
        u2 = u[2,i]#/2
        u3 = u[3,i]#/2

        R = sqrt(u1^2+u2^2+u3^2) +  tinyvalue
        sR = sin(R)/R
        cR = cos(R)
        a1 = u1*sR
        a2 = u2*sR
        a3 = u3*sR

        v[i] = SU2Matrix(cR + im*a3,
                         im*a1 - a2,
                         im*a1 + a2,
                         cR - im*a3)
        #v[i][1,1] = cos(R) + im*a3
        #v[i][1,2] = im*a1 + a2
        #v[i][2,1] = im*a1 - a2
        #v[i][2,2] = cos(R) - im*a3
    end
end






const sr3 = sqrt(3)
const sr3i = 1/sr3
const sr3i2 = 2*sr3i

function Gauge2Lie!(c::SU3AlgebraFields,x::SU3GaugeFields)
    NV = x.NV
    
    for i=1:NV
        x11 = x[i][1,1]
        x12 = x[i][1,2]
        x13 = x[i][1,3]
        x21 = x[i][2,1]
        x22 = x[i][2,2]
        x23 = x[i][2,3]
        x31 = x[i][3,1]
        x32 = x[i][3,2]
        x33 = x[i][3,3]

        c[1,i] = ( imag(x12) + imag(x21) )
        c[2,i] = ( real(x12) - real(x21) )
        c[3,i] = ( imag(x11) - imag(x22) )
        c[4,i] = ( imag(x13) + imag(x31) )
        c[5,i] = ( real(x13) - real(x31) )
        
        c[6,i] = ( imag(x23) + imag(x32) )
        c[7,i] = ( real(x23) - real(x32) )
        c[8,i] = sr3i *
                ( imag(x11) + imag(x22) -
                        2*imag(x33) )
    end
end


function Gauge2Lie!(c::SU2AlgebraFields,x::SU2GaugeFields)
    NV = x.NV

    for i=1:NV

        x11 = x[i][1,1]
        x12 = x[i][1,2]
        x21 = x[i][2,1]
        x22 = x[i][2,2]

        c[1,i] = (imag(x12)+imag(x21))
        c[2,i] = (real(x12)-real(x21))
        c[3,i] = (imag(x11)-imag(x22))

    end
end


function Gauge2Lie!(c::U1AlgebraFields,x::U1GaugeFields)
    NV = c.NV

    for i=1:NV
        c[1,i] = real(x[i][1,1])
    end
end



####################################################################


function trT(X::SU2Matrix)
    return 0.5 .* [X[2,1] + X[1,2],
        im*(X[1,2] - X[2,1]),
            X[1,1] - X[2,2]]
end

function trT!(V,X::SU2Matrix)
    V[1] = (X[2,1] + X[1,2])
    V[2] = im*(X[1,2] - X[2,1])
    V[3] = (X[1,1] - X[2,2])
end

function trT!(V::LieAlgebraFields{SU{N},aType}, X::GaugeFields{SU{N},eType}) where {N, eType, aType}
    for j in 1:X.NV
        M::SUMatrix{N} = X[j]
        for i in j+1:X.NV
            M *= X[i]
        end
        
        for i in 1:j-1
            M *= X[i]
        end
        trT!(view(V.a,:,j),M)
    end

    return V 
end

function tr(X::GaugeFields{SU{N},eType}) where {N, eType}
    M = X[1]
    for i in 2:X.NV
        M *= X[i]
    end
    return tr(M)
end

#=















function Base.setindex!(x::GaugeFields,v,i1,i2,i3,i4,i5,i6) 
    x.g[i1,i2,i3 .+ x.NDW,i4 + x.NDW,i5 + x.NDW,i6 + x.NDW] = v
end

function Base.getindex(x::GaugeFields,i1,i2,i3,i4,i5,i6)
    return x.g[i1,i2,i3 .+ x.NDW,i4 .+ x.NDW,i5 .+ x.NDW,i6 .+ x.NDW]
end

function clear!(a::GaugeFields)
    @. a.g = 0
    return 
end




function IdentityGauges(NC,NX,NY,NZ,NT,NDW)
    U = GaugeFields(NC,NDW,NX,NY,NZ,NT)
    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    @simd for ic=1:NC
                        U[ic,ic,ix,iy,iz,it] = 1 
                    end
                end
            end
        end
    end
    return U
end


function RandomGauges(NC,NX,NY,NZ,NT,NDW)
    U = GaugeFields(NC,NDW,NX,NY,NZ,NT)
    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    for j=1:NC
                        @simd for i=1:NC
                            U[i,j,ix,iy,iz,it] = rand()-0.5 + im*(rand()-0.5)
                        end
                    end
                end
            end
        end
    end
    normalize!(U)
    return U
end


function Base.similar(x::GaugeFields) 
    return GaugeFields(x.NC,x.NDW,x.NX,x.NY,x.NZ,x.NT)
end

function Base.similar(x::GaugeFields_1d) 
    return GaugeFields_1d(x.NC,x.NV) 
end

function Base.similar(x::Array{T,1}) where T <: GaugeFields
    xout = Array{T,1}(undef,4)
    for μ=1:4
        xout[μ] = similar(x[μ])
    end
    return xout
end

function Base.copyto!(Unew::GaugeFields,U::GaugeFields)
    Unew = deepcopy(U)
    return
end


function LinearAlgebra.tr(a::GaugeFields{SU{NC}}) where NC
    NX=a.NX
    NY=a.NY
    NZ=a.NZ
    NT=a.NT


    s = 0
    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    @simd for k=1:NC
                        s += a[k,k,ix,iy,iz,it]
                    end
                end
            end
        end
    end
    return s

end


function LinearAlgebra.tr(a::GaugeFields_1d{SU{NC}}) where NC
    NV=a.NV

    s = 0

    for i=1:NV

        @simd for k=1:NC
            s += a[k,k,i]
        end

    end
    return s

end






MProd(A::Matrix,B::T) where {T <: SUNElement} = A*B.M
arrMult(V::Vector{T}) where {T <: SUNElement} = T(reduce(MProd,V;init=[1 0 ; 0 1]))
arrMult(V::Vector{T}) where {T <: SMatrix} = reduce(*,V;init=[1 0 ; 0 1])

















##### 4D Gauge field #####
struct GaugeFields{T <: SUn}
    g::Array{ComplexF64,6}
    NX::Int64
    NY::Int64
    NZ::Int64
    NT::Int64
    NC::Int64
    NDW::Int64
    NV::Int64


    function GaugeFields(NC,NDW,NX,NY,NZ,NT)
        sutype = SU{NC}
            
        NV = NX*NY*NZ*NT
        g = zeros(ComplexF64,NC,NC,NX+2NDW,NY+2NDW,NZ+2NDW,NT+2NDW)
        return new{sutype}(g,NX,NY,NZ,NT,NC,NDW,NV)
    end

end

const U1GaugeFields  = GaugeFields{U1}
const SU2GaugeFields  = GaugeFields{SU2}
const SU3GaugeFields  = GaugeFields{SU3}
const SUNGaugeFields{N}  = GaugeFields{SU{N}}

struct Adjoint_GaugeFields{T <: SUn} 
    parent::GaugeFields{T}
end

const Adjoint_U1GaugeFields  = Adjoint_GaugeFields{U1}
const Adjoint_SU2GaugeFields  = Adjoint_GaugeFields{SU2}
const Adjoint_SU3GaugeFields  = Adjoint_GaugeFields{SU3}

function Base.adjoint(x::GaugeFields{T}) where T <: SUn
    Adjoint_GaugeFields{T}(x)
end

=#