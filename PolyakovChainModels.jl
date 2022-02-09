using GaugeSQM
using Plots
using Statistics

#prob = GaugeSQM.PolyakovChainProblem(SU2,1e-3,1.,0.5*(1. + 5*sqrt(3)*im),5;NTr=10)
prob = GaugeSQM.PolyakovChainProblem(SU3,1e-3,1.,2.,30;NTr=10,κ=0.1,μ=1.)

noRegs = Regulators(NoGaugeCooling(), NoDynamicStabilization())
regs_GC = Regulators(GaugeCooling(1.0,1), NoDynamicStabilization())
regs_GCAD = Regulators(GaugeCoolingAdaptive(1.0,0.8,10), NoDynamicStabilization())

#regs_I = Regulators(GaugeCooling(0.03,5), NoDynamicStabilization())
#regs_I_Ad = Regulators(GaugeCoolingAdaptive(1.0,0.8,1), NoDynamicStabilization())

DSRegs = Regulators(GaugeCoolingAdaptive(1.0,0.8,1), DynamicStabilization(1.0))

saveat = 0.01
@time sol_E = solve(prob,gEM(),noRegs)
@time sol_E_DS = solve(prob,gEM(),DSRegs,adaptive=true,saveat=saveat)
@time sol_E_GC = solve(prob,gEM(),regs_GC; adaptive=true, saveat=saveat)
@time sol_E_GCAD = solve(prob,gEM(),regs_GCAD; adaptive=false, saveat=saveat)

@time sol_I = solve(prob,gθEM(1.0),noRegs)
@time sol_I_DS = solve(prob,gθEM(1.0),DSRegs,adaptive=true,saveat=saveat)
@time sol_I_GC = solve(prob,gθEM(1.0),regs_GC; adaptive=true,saveat=saveat)
@time sol_I_GCAD = solve(prob,gθEM(1.0),regs_GCAD; adaptive=false, saveat=saveat)

dts = Float64[]
### Why adaptive stepsize is worse for implicit scheme?
cb = function (integrator)
    append!(dts,integrator.dt)    
end
dts_I = dts
dts_E = dts

begin
fig = plot(yaxis=:log,ylim=[7e-5,1e-3])
plot!(fig,dts_E)
plot!(fig,dts_I)
end

termTime = 1 
sol_E = reshape(sol_E[:,floor(Int64,1/saveat)*termTime:end],:)
sol_E_DS = reshape(sol_E_DS[:,floor(Int64,1/saveat)*termTime:end],:)
sol_E_GC = reshape(sol_E_GC[:,floor(Int64,1/saveat)*termTime:end],:)#[1:10:end]
sol_E_GCAD = reshape(sol_E_GCAD[:,floor(Int64,1/saveat)*termTime:end],:)#[1:10:end]
sol_E = sol_E[abs.(real(sol_E)) .< 1e6]
sol_E = sol_E[abs.(imag(sol_E)) .< 1e6]
sol_E_GC = sol_E_GC[abs.(real(sol_E_GC)) .< 1e10]
sol_E_GC = sol_E_GC[abs.(imag(sol_E_GC)) .< 1e10]
sol_E_GCAD = sol_E_GCAD[abs.(real(sol_E_GCAD)) .< 1e6]
sol_E_GCAD = sol_E_GCAD[abs.(imag(sol_E_GCAD)) .< 1e6]


sol_I = reshape(sol_I[:,floor(Int64,1/saveat)*termTime:end],:)
sol_I_DS = reshape(sol_I_DS[:,floor(Int64,1/saveat)*termTime:end],:)
sol_I_GC = reshape(sol_I_GC[:,floor(Int64,1/saveat)*termTime:end],:)
sol_I_GCAD = reshape(sol_I_GCAD[:,floor(Int64,1/saveat)*termTime:end],:)

begin
ylim = [-10,10] 
plot(real(sol_I_GC_AD),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_I_GC_AD)),ylim=[1e-3,10],yaxis=:log,label="Implicit Re GC_AD")
plot!(real(sol_E_GCAD),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_E_GCAD)),ylim=[1e-3,10],yaxis=:log,label="Explicit Re GC_AD")
plot!(real(sol_I_GC),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_I_GC)),ylim=[1e-3,10],yaxis=:log,label="Implicit Re GC")
plot!(real(sol_E_GC),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_E_GC)),ylim=[1e-3,10],yaxis=:log,label="Explicit Re GC")
#plot!(real(sol_I),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_I)),ylim=[1e-3,10],yaxis=:log,label="Implicit Re")
#plot!(real(sol_E),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_E)),ylim=[1e-3,10],yaxis=:log,label="Explicit Re")

plot!(imag(sol_I_GC_AD),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_I_GC_AD)),ylim=[1e-3,10],yaxis=:log,label="Implicit Im GC_AD")
plot!(imag(sol_E_GCAD),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_E_GCAD)),ylim=[1e-3,10],yaxis=:log,label="Explicit Im GC_AD")
plot!(imag(sol_I_GC),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_I_GC)),ylim=[1e-3,10],yaxis=:log,label="Implicit Im GC")
plot!(imag(sol_E_GC),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_E_GC)),ylim=[1e-3,10],yaxis=:log,label="Explicit Im GC")
#plot!(imag(sol_E),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_E)),ylim=[1e-3,10],yaxis=:log,label="Explicit Im")
#plot!(imag(sol_I),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_I)),ylim=[1e-3,10],yaxis=:log,label="Implicit Im")
end

begin
ylim = [-15,15] 
plot(real(sol_I_GC_AD),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_I_GC_AD)),ylim=[1e-3,10],yaxis=:log,label="Implicit Re GC_AD")
plot!(imag(sol_I_GC_AD),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_I_GC_AD)),ylim=[1e-3,10],yaxis=:log,label="Implicit Im GC_AD")
end

begin
    ylim = [-5,5] 
    plot(real(sol_I_GC),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_I_GC_AD)),ylim=[1e-3,10],yaxis=:log,label="Implicit Re GC_AD")
    plot!(imag(sol_I_GC),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_I_GC_AD)),ylim=[1e-3,10],yaxis=:log,label="Implicit Im GC_AD")
    end

begin
    ylim = [-5,5] 
    plot(real(sol_E_GCAD),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_E_GCAD)),ylim=[1e-3,10],yaxis=:log,label="Explicit Re GC_AD")
    plot!(imag(sol_E_GCAD),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_E_GCAD)),ylim=[1e-3,10],yaxis=:log,label="Explicit Im GC_AD")
end


begin
    ylim = [-5,5] 
    
    plot(real(sol_I),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_I)),ylim=[1e-3,10],yaxis=:log,label="Implicit Re")
    plot!(imag(sol_I),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_I)),ylim=[1e-3,10],yaxis=:log,label="Implicit Im")
    #plot!(real(sol_E),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_E)),ylim=[1e-3,10],yaxis=:log,label="Explicit Re")
    #plot!(imag(sol_E),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_E)),ylim=[1e-3,10],yaxis=:log,label="Explicit Im")
end

begin
    ylim = [-50,50] 
    
    plot(real(sol_I_GC),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_I_GC)),ylim=[1e-3,10],yaxis=:log,label="Implicit Re GC")
    plot!(imag(sol_I_GC),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_I_GC)),ylim=[1e-3,10],yaxis=:log,label="Implicit Im GC")
    plot!(real(sol_E_GC),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_E_GC)),ylim=[1e-3,10],yaxis=:log,label="Explicit Re GC")
    plot!(imag(sol_E_GC),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_E_GC)),ylim=[1e-3,10],yaxis=:log,label="Explicit Im GC")
end

begin
    ylim = [-20,20] 
    
    plot(real(sol_I_DS),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_I_DS)),ylim=[1e-3,10],yaxis=:log,label="Implicit Re DS")
    plot!(imag(sol_I_DS),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_I_DS)),ylim=[1e-3,10],yaxis=:log,label="Implicit Im DS")
    plot!(real(sol_E_DS),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_E_DS)),ylim=[1e-3,10],yaxis=:log,label="Explicit Re DS")
    plot!(imag(sol_E_DS),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_E_DS)),ylim=[1e-3,10],yaxis=:log,label="Explicit Im DS")
end

begin
    ylim = [-5,5] 
    
    plot(real(sol_E_DS),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_E_DS)),ylim=[1e-3,10],yaxis=:log,label="Explicit Re DS")
    plot!(imag(sol_E_DS),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_E_DS)),ylim=[1e-3,10],yaxis=:log,label="Explicit Im DS")
    plot!(real(sol_E_GC),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_E_GC)),ylim=[1e-3,10],yaxis=:log,label="Explicit Re GC")
    plot!(imag(sol_E_GC),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_E_GC)),ylim=[1e-3,10],yaxis=:log,label="Explicit Im GC")
    plot!(real(sol_E_GCAD),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_E_GCAD)),ylim=[1e-3,10],yaxis=:log,label="Explicit Re GCAD")
    plot!(imag(sol_E_GCAD),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_E_GCAD)),ylim=[1e-3,10],yaxis=:log,label="Explicit Im GCAD")
    plot!(real(sol_I_DS),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_I_DS)),ylim=[1e-3,10],yaxis=:log,label="Implicit Re DS")
    plot!(imag(sol_I_DS),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_I_DS)),ylim=[1e-3,10],yaxis=:log,label="Implicit Im DS")
    plot!(real(sol_I_GC),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_I_GC)),ylim=[1e-3,10],yaxis=:log,label="Implicit Re GC")
    plot!(imag(sol_I_GC),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_I_GC)),ylim=[1e-3,10],yaxis=:log,label="Implicit Im GC")
    plot!(real(sol_I_GCAD),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(real(sol_I_GCAD)),ylim=[1e-3,10],yaxis=:log,label="Implicit Re GC")
    plot!(imag(sol_I_GCAD),xlim=ylim,seriestype=:stephist, norm=true,bins=GaugeSQM.get_nr_bins(imag(sol_I_GCAD)),ylim=[1e-3,10],yaxis=:log,label="Implicit Im GC")
end


mean(mean(sol_E_GC,dims=2))
mean(sol_I_GC_AD)

sol

begin
fig = plot(ylim=[-10,10])
plot!(fig,real(sol_E_GCAD[1300:2000]), label="Re Explicit")
plot!(fig,real(sol_I_GC_AD[1300:2000]), label="Re Implicit")
end

begin
fig = plot(ylim=[-100,100])
plot!(fig,imag(sol_E_GCAD), label="Im Explicit")
plot!(fig,imag(sol_I_GC_AD), label="Im Implicit")
end

sol_E
sol_E = [s for s in sol_E if !isnan(s)]



using Measurements

begin
exact = 2.09
N16PC_E = Measurement{Float64}[]
N16PC_I = Measurement{Float64}[]
prob = GaugeSQM.PolyakovChainProblem(SU3,1e-4,30.,2.,10;NTr=20,κ=0.1,μ=1.)
αs = collect(0.03:0.2:1.53)
for α in αs
    println("Starting on ", α)
    regs = Regulators(GaugeCooling(α,5), NoDynamicStabilization())
    
    saveat = 0.005
    termTime = 1
    @time sol_E_GC = solve(prob,gEM(),regs; adaptive=true, saveat=saveat)
    M_E_GC = mean(sol_E_GC[:,floor(Int64,1/prob.dt)*termTime:end],dims=2)
    append!(N16PC_E,real(mean(M_E_GC)) ± real(std(M_E_GC))/sqrt(length(M_E_GC)))
    #append!(N16PC_E,real(mean(sol_E_GC)) ± real(std(sol_E_GC))/sqrt(length(sol_E_GC)))

    @time sol_I_GC = solve(prob,gθEM(1.0),regs; adaptive=true, saveat=saveat)
    M_I_GC = mean(sol_I_GC[:,floor(Int64,1/prob.dt)*termTime:end],dims=2)
    append!(N16PC_I,real(mean(M_I_GC)) ± real(std(M_I_GC))/sqrt(length(M_I_GC)))
    #append!(N16PC_I,real(mean(sol_I_GC)) ± real(std(sol_I_GC))/sqrt(length(sol_I_GC)))
end
end

begin
    fig = plot()
    hline!(fig,[exact])
    scatter!(fig,N16PC_E)
    scatter!(fig,N16PC_I)
end


















A = LieAlgebraFields(ComplexF64,2,10)

sum(norm,A)

for i in A
    @show i
end

GaugeSQM.clear!(A)

Aview = view(A.a,:,1)

A[1:2,1] = [2.,1.]

A

randn(size(A.a))

GaugeSQM.mul!(B,2.,A)
A*B
GaugeSQM.gauss_distribution_lie!(A)
A

B.a .= convert(typeof(A.a),randn(size(A.a)))



B
B = similar(A)


g = zeros(SUMatrix{2},10)
[SUMatrix{2}(0,0,0,0) for i in 1:10]


SU2Matrix(0.,0.,0.,0. + im) *
SU2Matrix(0.,0.,0.,0. + im)

copy(GF.g)

GF1 = GaugeFields(2,10)
GF2 = GaugeSQM.IdentityGauges(2,1)
GF3 = GaugeSQM.RandomGauges(2,10)

for i in 1:1000000
    GaugeCoolingUpdate!(GF3,A,similar(A),similar(GF3),similar(GF3),similar(GF3))
end

GF3


sum(tr,GF3)

using LinearAlgebra
mul!(GF1,GF2,GF3)
GF1
GF3
GaugeSQM.trT(GF)


(1/2)*GaugeSQM.tr(GF2[1])


GaugeSQM.substitute!(GF1,GF2)

GF1

GF2[1] = SUMatrix{2}(2000,0,0,0)
GF2


M = GF[1]

G = M

AF[1]
GF[1]


AF[1] = SUMatrix{2}(2,1,1,1)

GF

copy(GF)

GaugeSQM.Gauge2Lie!(A,GF)
A

GF[1][1,1] = 1

GF[5] = SUMatrix{GF.NC}(1. + im,2.,3.,4.)
GF = GaugeSQM.IdentityGauges(2,10)
GF = GaugeSQM.RandomGauges(2,10)



eltype(GF.g)(1,2,3,4)

GF
GaugeSQM.clear!(GF)


f(G::GaugeFields{S,T}) where {S <: SU2, T} = T


f(GF)








@. GF.g = GaugeSQM.identity(GF.g[1])

typeof(zero(eltype(GF.g)))
typeof(identity(GF.g[1]))
typeof(GF.g)

eltype(GF.g)
typeof(identity(GF.g[1]))

GF

V = [SUMatrix{GF.NC}(1. + im,2.,3.,4.) for i in 1:3]

eltype(V)
SUMatrix{2}(1. + im,2.,3.,4.) 
zeros(SUMatrix{2},10)
[SUMatrix{2}(zeros(2^2)) for i in 1:10]







SUMatrix2{2}(0,0,0,0)


const T{N,L} = GaugeSQM.SMatrix{N,N,Complex{Real},L}

T{2,4}(1.,1.,1.,1.)