
using Plots, ColorTypes
using CSV
using LaTeXStrings
using DataFrames, Statistics

orange = RGB(230/255,159/255,0/255)  
vermillion = RGB(213/255, 94/255, 0/255) 
bluishgreen = RGB(0/255, 158/255, 115/255)  
default(size = (800, 600))
default(
    guidefont = font(15),   # Axis labels
    tickfont = font(12),    # Tick labels
    legendfont = font(10)   # Legend text
)

data = CSV.File("/home/cornell/BEMJulia/MarineHydro.jl/paper/Plots/11_heuristics_dx_DELhommeau.csv") |> DataFrame

plot(data.dx, data.A11_grad_r, label=L"\frac{\partial A_{11}(\omega = 1.03\,\mathrm{rad/s})}{\partial r_2}", marker=:circle, lw=2, color=bluishgreen,legend=:right)
plot!(data.dx, data.B11_grad_r, label=L"\frac{\partial B_{11}(\omega = 1.03\,\mathrm{rad/s})}{\partial r_2}", marker=:square, lw=2, color=vermillion)
vline!([5], label="PWA heuristics[Singh (2013)]", lw=2, linestyle=:dash, color=orange)

# Labels and Title
xlabel!("separation distance (x) [m]")
ylabel!("Sensitivity values")
title!("Sensitivity of coeffficients of a sphere with nearby sphere's dimension ")

# Save the plot
savefig("/home/cornell/BEMJulia/MarineHydro.jl/paper/Plots/11_heuristics_dx.pdf")


#switch for damping and added mass here and corresponding labels , columns below
data = CSV.File("/home/cornell/BEMJulia/MarineHydro.jl/paper/Plots/added_mass_data_dimensionless.csv") |> DataFrame

# Unique values for dimensionless parameters
dx_r_ratios = unique(data.dx_r_ratio)
kr_values = unique(data.kr)

dx_values = unique(data.dx_r_ratio)
kr_values =  unique(data.kr)

grad_r_matrix = reshape(data.grad_r, length(dx_values), length(kr_values))

# Normalize function
function normalize(matrix)
    min_val = minimum(matrix)
    max_val = maximum(matrix)
    println("Maximum: $max_val")
    println("Minimum: $min_val")
    return (matrix .- min_val) ./ (max_val - min_val)
end

# Normalize gradient matrices
grad_r_matrix = normalize(grad_r_matrix)


# Plot heatmap using dimensionless parameters
p1 = heatmap(
    dx_r_ratios, kr_values, grad_r_matrix,
    ylabel="Kr",
    xlabel="x/r",
    color=:magma,
    xguidefontsize=14,
    yguidefontsize=14,
    titlefontsize=16,
    tickfontsize=8,
    colorbar_title="∂A/∂r" #switch
)

savefig("/home/cornell/BEMJulia/MarineHydro.jl/paper/Plots/added_mass_dimensionless_grad_dr.pdf")

