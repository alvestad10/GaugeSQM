

# GaugeSQM.jl
<!---
[![DOI](https://zenodo.org/badge/552342814.svg)](https://zenodo.org/badge/latestdoi/552342814)
-->

To run the code you need a version of Julia installed, then you can make separate script, and run it as `julia --project=. script.jl` (add -t X if more than one trajectory is run as this will parallelize the rund, X is the number of threads) or line by line using the Julia vscode extension. Before you run the code follow the Instantite section below to setup the necessary packages.

## Instantiate

To initialize the project run these comments inside the Julia REPL (From inside the project directory)
```julia
    import Pkg
    Pkg.activate(".")
    Pkg.instantiate()
```
For more information see: https://docs.julialang.org/en/v1/stdlib/Pkg/

Now all dependencies should be downloaded and the code is ready to be run.


#  Overview of code

Define a model by, e.g., for the Polyakov chain model; GaugeSQM.PolyakovChainModel(), where the first input is the group SU(N), which can be SU2 or SU3. Then make the regulator object that can contain a GaugeCooling type and a SynamicalStabilization type. Then solve system using the solve() function which takes as input the problem and the regulator object. 

The ouput of this function is the action observable. TODO: Make it easier to return more observables, or all the configurations. (It is possible to change this to returning the unitarity norm, by changing the observable funciton in the "src/Model.jl" file)

TODO: Fix the line where we need to reshape the whole thing before plotting when using more than one trajectory.

Under is some example scripts:

## Example 1
First example is to compare the fixed step-size explicit scheme to the implicit scheme

```julia
using GaugeSQM
using Plots

prob = GaugeSQM.PolyakovChainProblem(SU2;dt=1e-3,tspan=30.,β = 0.5*(1. + sqrt(3)*im), NLinks = 30, NTr=1)

noRegs = Regulators(NoGaugeCooling(), NoDynamicStabilization())

termTime = 2; saveat=1e-2
sol_E = solve(prob,gEM(),noRegs,saveat=saveat)
sol_I = solve(prob,gθEM(1.0),noRegs,saveat=saveat)
sol_E = reshape(sol_E[:,floor(Int64,1/saveat)*termTime:end],:)
sol_I = reshape(sol_I[:,floor(Int64,1/saveat)*termTime:end],:)

# We do a check if any trajectories have diverges or failed, then we get nan, 0 or a large number.
begin
    _sol_I = sol_I[isnan.(sol_I) .!= 1 .&& abs.(sol_I) .< 100 .&& abs.(sol_I) .!= 0.]
    _sol_E = sol_E[isnan.(sol_E) .!= 1 .&& abs.(sol_E) .< 100 .&& abs.(sol_E) .!= 0.]
    NI = length(_sol_I)/length(sol_I)
    NE = length(_sol_E)/length(sol_E)
    @info "Explicit: $NE ratio of events expected, Implicit: $NI ratio of events expected" 
    
    fig = plot(;xlim=[-5,5],ylim=[1e-3,10],yaxis=:log)
    plot!(fig,real(_sol_I),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(real(_sol_I)),label="Implicit Re")
    plot!(fig,imag(_sol_I),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(imag(_sol_I)),label="Implicit Im")
    plot!(fig,real(_sol_E),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(real(_sol_E)),label="Explicit Re")
    plot!(fig,imag(_sol_E),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(imag(_sol_E)),label="Explicit Im")
end
```


## Example 2

In this example we add adaptive stepsize and compare with and without gauge cooling, where the implciit scheme is without gauge cooling and the explicit is with gauge cooling;

```julia
using GaugeSQM
using Plots

prob = GaugeSQM.PolyakovChainProblem(SU2;dt=1e-3,tspan=30.,β = 0.5*(1. + sqrt(3)*im), NLinks = 30, NTr=10)

noRegs = Regulators(NoGaugeCooling(), NoDynamicStabilization())
regs_GC = Regulators(GaugeCooling(1.0,10), NoDynamicStabilization())

termTime = 2; saveat=1e-2
@time sol_E = solve(prob,gEM(),regs_GC,saveat=saveat)
@time sol_I_GC = solve(prob,gθEM(1.0),regs_GC,saveat=saveat)
@time sol_I = solve(prob,gθEM(1.0),noRegs,saveat=saveat)
sol_E = reshape(sol_E[:,floor(Int64,1/saveat)*termTime:end],:)
sol_I = reshape(sol_I[:,floor(Int64,1/saveat)*termTime:end],:)
sol_I_GC = reshape(sol_I_GC[:,floor(Int64,1/saveat)*termTime:end],:)

# We do a check if any trajectories have diverges or failed, then we get nan, 0 or a large number.
begin
    _sol_I = sol_I[isnan.(sol_I) .!= 1 .&& abs.(sol_I) .< 100 .&& abs.(sol_I) .!= 0.]
    _sol_I_GC = sol_I_GC[isnan.(sol_I_GC) .!= 1 .&& abs.(sol_I_GC) .< 100 .&& abs.(sol_I_GC) .!= 0.]
    _sol_E = sol_E[isnan.(sol_E) .!= 1 .&& abs.(sol_E) .< 100 .&& abs.(sol_E) .!= 0.]
    NI = length(_sol_I)/length(sol_I)
    NI_GC = length(_sol_I_GC)/length(sol_I_GC)
    NE = length(_sol_E)/length(sol_E)
    @info "Explicit with GC: $NE ratio of events expected, Implicit: $NI ratio of events expected, Implicit with GC: $NI_GC ratio of events expected" 
    
    fig = plot(;xlim=[-5,5],ylim=[1e-3,10],yaxis=:log)
    plot!(fig,real(_sol_I),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(real(_sol_I)),label="Implicit Re")
    plot!(fig,imag(_sol_I),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(imag(_sol_I)),label="Implicit Im")
    plot!(fig,real(_sol_I_GC),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(real(_sol_I_GC)),label="Implicit GC Re")
    plot!(fig,imag(_sol_I_GC),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(imag(_sol_I_GC)),label="Implicit GC Im")
    plot!(fig,real(_sol_E),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(real(_sol_E)),label="Explicit GC Re")
    plot!(fig,imag(_sol_E),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(imag(_sol_E)),label="Explicit GC Im")
end
```


## Example 2

In this example we add adaptive stepsize and compare with and without gauge cooling, where the implicit scheme is without gauge cooling and the explicit is with gauge cooling;

TODO: Fix bug where DS failes

```julia
using GaugeSQM
using Plots

prob = GaugeSQM.PolyakovChainProblem(SU2;dt=1e-3,tspan=30.,β = 0.5*(1. + sqrt(3)*im), NLinks = 30, NTr=10)

noRegs = Regulators(NoGaugeCooling(), NoDynamicStabilization())
regs_GC = Regulators(GaugeCooling(1.0,10), NoDynamicStabilization())
regs_GCAD = Regulators(GaugeCoolingAdaptive(1.0,0.8,10), NoDynamicStabilization())
regs_GC_DS = Regulators(GaugeCooling(1.0,1), DynamicStabilization(1.0))
regs_DS = Regulators(NoGaugeCooling(), DynamicStabilization(1.0))


termTime = 2; saveat=1e-2
@time sol_I = solve(prob,gθEM(1.0),noRegs,saveat=saveat)
@time sol_E_GC = solve(prob,gEM(),regs_GC,saveat=saveat)
@time sol_E_GCAD = solve(prob,gEM(),regs_GCAD,saveat=saveat)
@time sol_E_GC_DS = solve(prob,gEM(),regs_GC_DS,saveat=saveat)
@time sol_I_DS = solve(prob,gθEM(1.0),regs_DS,saveat=saveat)
sol_I = reshape(sol_I[:,floor(Int64,1/saveat)*termTime:end],:)
sol_E_GC = reshape(sol_E_GC[:,floor(Int64,1/saveat)*termTime:end],:)
sol_E_GCAD = reshape(sol_E_GCAD[:,floor(Int64,1/saveat)*termTime:end],:)
sol_E_GC_DS = reshape(sol_E_GC_DS[:,floor(Int64,1/saveat)*termTime:end],:)
sol_I_DS = reshape(sol_I_DS[:,floor(Int64,1/saveat)*termTime:end],:)

# We do a check if any trajectories have diverges or failed, then we get nan, 0 or a large number.
begin
    _sol_I = sol_I[isnan.(sol_I) .!= 1 .&& abs.(sol_I) .< 100 .&& abs.(sol_I) .!= 0.]
    _sol_E_GC = sol_E_GC[isnan.(sol_E_GC) .!= 1 .&& abs.(sol_E_GC) .< 100 .&& abs.(sol_E_GC) .!= 0.]
    _sol_E_GCAD = sol_E_GCAD[isnan.(sol_E_GCAD) .!= 1 .&& abs.(sol_E_GCAD) .< 100 .&& abs.(sol_E_GCAD) .!= 0.]
    _sol_E_GC_DS = sol_E_GC_DS[isnan.(sol_E_GC_DS) .!= 1 .&& abs.(sol_E_GC_DS) .< 100 .&& abs.(sol_E_GC_DS) .!= 0.]
    _sol_I_DS = sol_I_DS[isnan.(sol_I_DS) .!= 1 .&& abs.(sol_I_DS) .< 100 .&& abs.(sol_I_DS) .!= 0.]
    NI = length(_sol_I)/length(sol_I)
    NE_GC = length(_sol_E_GC)/length(sol_E_GC)
    NE_GCAD = length(_sol_E_GCAD)/length(sol_E_GCAD)
    NE_GC_DS = length(_sol_E_GC_DS)/length(sol_E_GC_DS)
    NE_DS = length(_sol_I_DS)/length(sol_I_DS)
    @info "Implicit: $NI ratio of events expected, \n Explicit with GC: $NE_GC ratio of events expected, \n Explicit with adaptive GC: $NE_GCAD ratio of events expected, \n Explicit with adaptive GC and DS: $NE_GC_DS ratio of events expected, \n Implicit with DS: $NE_DS ratio of events expected" 
    
    fig = plot(;xlim=[-5,5],ylim=[1e-3,10],yaxis=:log)
    plot!(fig,real(_sol_I),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(real(_sol_I)),label="Implicit Re")
    plot!(fig,imag(_sol_I),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(imag(_sol_I)),label="Implicit Im")
    plot!(fig,real(_sol_E_GC),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(real(_sol_E_GC)),label="Explicit GC Re")
    plot!(fig,imag(_sol_E_GC),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(imag(_sol_E_GC)),label="Explicit GC Im")
    plot!(fig,real(_sol_E_GCAD),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(real(_sol_E_GCAD)),label="Explicit GCAD Re")
    plot!(fig,imag(_sol_E_GCAD),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(imag(_sol_E_GCAD)),label="Explicit GCAD Im")
    plot!(fig,real(_sol_E_GC_DS),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(real(_sol_E_GC_DS)),label="Explicit GCAD DS Re")
    plot!(fig,imag(_sol_E_GC_DS),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(imag(_sol_E_GC_DS)),label="Explicit GCAD DS Im")
    plot!(fig,real(_sol_I_DS),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(real(_sol_I_DS)),label="Implciit DS Re")
    plot!(fig,imag(_sol_I_DS),seriestype=:stephist, norm=true, bins=GaugeSQM.get_nr_bins(imag(_sol_I_DS)),label="Implicit DS Im")
end
```