
using Plots, ColorTypes
using CSV
using LaTeXStrings
using DataFrames, Statistics

orange = RGB(230/255,159/255,0/255)  
vermillion = RGB(213/255, 94/255, 0/255) 
bluishgreen = RGB(0/255, 158/255, 115/255)  

Wudata = CSV.File("/home/cornell/BEMJulia/BEM.jl/paper/Plots/heuristics_dx.csv") |> DataFrame
Deldata = CSV.File("/home/cornell/BEMJulia/BEM.jl/paper/Plots/heuristics_dx_DELhommeau.csv") |> DataFrame


Wudata.Aerror = abs.(Wudata.A12_grad_r .- Deldata.A12_grad_r)./Deldata.A12_grad_r
Wudata.Berror = abs.(Wudata.B12_grad_r .- Deldata.B12_grad_r)./Deldata.B12_grad_r
plot(Wudata.dxs, Wudata.Aerror, label=L"Error A(w)", marker=:circle, lw=2, color=bluishgreen,legend=:topright)
plot!(Wudata.dxs, Wudata.Berror, label=L"Error B(w)", marker=:square, lw=2, color=vermillion)
hline!([0], label="", lw=2, linestyle=:dash, color=orange)
hline!([-0.01], label="", lw=2, linestyle=:dash, color=orange)
hline!([0.01], label="", lw=2, linestyle=:dash, color=orange)

# Labels and Title
xlabel!("separation distance (x)")
ylabel!("Sensitivity values (normalized)")
title!("Comparison of Wu vs Delhommeau sensitivities")

# Save the plot
savefig("/home/cornell/BEMJulia/BEM.jl/paper/Plots/heuristics_Wu_DEL.pdf")