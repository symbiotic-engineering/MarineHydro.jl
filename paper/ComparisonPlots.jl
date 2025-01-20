## COEFFICIENTS CHECKS - DO NOT DELETE even if commented out; MAKE TEST
using LinearAlgebra,Statistics,Plots
using BEM
using ColorTypes
include("MeshGradients_singlebody.jl")
include("HulmeData.jl")
include("BemProgram.jl")

# Define colorblind-friendly colors using hex values
orange = RGB(230/255,159/255,0/255)  
vermillion = RGB(213/255, 94/255, 0/255) 
bluishgreen = RGB(0/255, 158/255, 115/255) 

default(
    guidefont = font(15),   # Axis labels
    tickfont = font(12),    # Tick labels
    legendfont = font(13)   # Legend text
)

# ffix resoluution first in meshGradients.jl - for gradient calculations they need to be decreased as it takes a long time with this naive implementation.
# Combine colors in a plot or visualization
using Plots
heave = [0,0,1]
surge = [1,0,0]
# #compare with analytical and BEMjulia
rho = 1023
omegas = sqrt.(K_heave .* 9.8)

results_damping =  [damping_program(1.0,wn,heave)/wn for wn in omegas ] ./ ((2/3)*pi*1*1023);

results_amass =  [added_mass_program(1.0,wn,heave) for wn in omegas ] ./ ((2/3)*pi*1*1023);

# # # Plotting results with respect to wavenumbers
plot(K_heave,  results_damping, xlabel = "kr", ylabel = "B(kr)",
      label = "new solver", marker = :circle, linewidth = 2,linecolor = vermillion,markercolor = vermillion)

plot!(K_heave,B_heave, xlabel = "kr",
      label = "analytical", marker = :circle, linewidth = 2,linecolor = bluishgreen,markercolor = bluishgreen, linestyle=:dot )

savefig("/home/cornell/BEMJulia/MarineHydro.jl/paper/Plots/HemisphereCoefficientComparison_damping.pdf")

plot(K_heave,  results_amass, xlabel = "kr",ylabel = "A(kr)",
     label = "new solver", marker = :diamond, linewidth = 2,linecolor = vermillion,markercolor = vermillion)

plot!(K_heave,A_heave, xlabel = "kr", 
     label = "analytical", marker = :diamond, linewidth = 2,linecolor = bluishgreen,markercolor = bluishgreen , linestyle=:dot )

savefig("/home/cornell/BEMJulia/MarineHydro.jl/paper/Plots/HemisphereCoefficientComparison_amass.pdf")
  


# # ##### SUURGE ####

omegas = sqrt.(K_surge .* 9.8)
results_damping =  [damping_program(1.0,wn,surge)/wn for wn in omegas ] ./ ((2/3)*pi*1*1023);
results_amass =  [added_mass_program(1.0,wn,surge) for wn in omegas ] ./ ((2/3)*pi*1*1023);

# # Plotting results with respect to wavenumbers
plot(K_surge,  results_damping, xlabel = "kr", ylabel = "B(kr)",
     label = "new solver", marker = :circle, linewidth = 2,linecolor = vermillion,markercolor = vermillion)

plot!(K_surge, B_surge, xlabel = "kr",
     label = "analytical", marker = :circle, linewidth = 2,linecolor = bluishgreen,markercolor = bluishgreen, linestyle=:dot )

savefig("/home/cornell/BEMJulia/MarineHydro.jl/paper/Plots/Surge_HemisphereCoefficientComparison_damping.pdf")

plot(K_surge,  results_amass, xlabel = "kr",ylabel = "A(kr)",
     label = "new solver", marker = :diamond, linewidth = 2,linecolor = vermillion,markercolor = vermillion)

plot!(K_surge, A_surge, xlabel = "kr", 
     label = "analytical", marker = :diamond, linewidth = 2,linecolor = bluishgreen,markercolor = bluishgreen , linestyle=:dot )

savefig("/home/cornell/BEMJulia/MarineHydro.jl/paper/Plots/Surge_HemisphereCoefficientComparison_amass.pdf")




# #             ##### Diffraction results comparing with Capytaine ####

using PyCall
using Plots

# Import Capytaine
cpt = pyimport("capytaine")
r = 1.0
# Create the mesh and body
cptmesh = cpt.mesh_sphere(radius=r, center=(0, 0, 0), resolution=(14, 14)).immersed_part()
cptbody = cpt.FloatingBody(cptmesh, name="sphere")
cptbody.add_translation_dof(name="Heave")

# Define the wave numbers and corresponding frequencies
K_heave_diff = [0.1,0.5, 1.0, 2.0,3.0, 5.0,7.0,8.0, 10.0]
omegas = sqrt.(K_heave_diff .* 9.8)

# Create the test matrix
xr = pyimport("xarray")
test_matrix = xr.Dataset(coords=Dict("omega" => omegas, "wave_direction" => [0.0]))
results = cpt.BEMSolver().fill_dataset(test_matrix, cptbody, method="direct")
Froude_heave =  vec(results["Froude_Krylov_force"].values)
Diff_heave =  vec(results["diffraction_force"].values)

# Display the results
Diff_heave = Diff_heave ./ ((2/3)*pi*1*1023)
Froude_heave = Froude_heave ./ ((2/3)*pi*1*1023)

# Non-dimensionalize the results
non_dimensional_const = ((2/3) * pi * r * 1023)
julia_diff_omega = [diffraction_program(w, r, heave) for w in omegas] ./ non_dimensional_const

# Plot the results
plot(K_heave_diff, abs.(Diff_heave), marker=:circle, label="capytaine", color=bluishgreen, linestyle=:dash)
plot!(K_heave_diff, abs.(julia_diff_omega), marker=:square, label="new solver", color=vermillion, linestyle=:solid)

xlabel!("kr")
ylabel!("Diffraction Force")
plot!(legend=:topright, grid=true)

savefig("/home/cornell/BEMJulia/MarineHydro.jl/paper/Plots/diffraction_comparison.pdf")

## Frroude Krylov

function froudeKrylovProgram(r,ω,dof)
    mesh = differentiableMesh(r) #fd
    force = FroudeKrylovForce(mesh::Mesh, ω,dof)
    return force
end

non_dimensional_const = ((2/3) * pi * r * 1023)
julia_froude_omega = [froudeKrylovProgram(r, w,heave) for w in omegas] ./ non_dimensional_const
# Plot the results
plot(K_heave_diff, abs.(Froude_heave), marker=:circle, label="capytaine", color=bluishgreen, linestyle=:dash)
plot!(K_heave_diff, abs.(julia_froude_omega), marker=:square, label="new solver", color=vermillion, linestyle=:solid)

xlabel!("kr")
ylabel!("Froude Krylov Force")
plot!(legend=:topright, grid=true)
savefig("/home/cornell/BEMJulia/MarineHydro.jl/paper/Plots/froudeKrylov_comparison.pdf")




###### =====comparing the value and gradients for sphere with no frequency =====#
omega = 1.0

# using Printf

# Constants
const π = 3.141592653589793 
rho = 1023  

radii = [1,2,3,4,5]

# #  reference added mass when no wave frequency
function reference_added_mass(radius)
    return (2 / 3) * π * rho * radius^3
end

# # Function to calculate the differential of the added mass with respect to radius
function reference_diff_added_mass(radius)
    return 2 * π * rho * radius^2
end


results_amass_solver =  [added_mass_program(r,omega,surge) for r in radii ] #./  ((2/3)*pi*1*1023);
results_amass_solver =  [reference_added_mass(r) for r in radii ]

analy_grady = [reference_diff_added_mass(r) ./  ((2/3)*pi*r^3*rho) for r in radii];


# _gradients = Array{Any}(undef, length(radii))

# # for (i, r) in enumerate(radii)
# #         grads = Zygote.gradient(X -> rankine_program(X[1], X[2]), [r, omega, surge])[1][1]
# #         _gradients[i] = grads # completelly submerged not just half floating free surface at 0
# #     end

_gradients = [6555.47838568984, 26221.913542847837, 58999.30547136783, 104887.65417050153, 163886.9596426078]

_gradients = _gradients ./ ((2/3).*pi.*(radii.^3).*rho)

##_grad_coarse =  [6529.024840896878, 26116.099363590965, 58761.22356808913, 104464.3974542766, 163225.6210222983] ./ ((2/3).*pi.*radii.*1000);
_grad_fine =  [6552.7934287808785, 26211.17371512438, 58975.14085905893, 104844.69486049563, 163819.83571945396] ./ ((2/3).*pi.*(radii.^3).*rho);

#plot(radii, _grad_coarse, xlabel="r [m]", ylabel="∂A_∂ω", label="AD (coarse mesh)",marker = "*")
plot(radii, _gradients , xlabel="r [m]", ylabel="∂A_∂r", label="AD",marker = "*",color = vermillion)
plot!(radii,analy_grady  , xlabel="r [m]", label="analytical",marker = "*" , color = bluishgreen)
savefig("/home/cornell/BEMJulia/MarineHydro.jl/paper/Plots/analy_ad_AMass_surge.pdf")