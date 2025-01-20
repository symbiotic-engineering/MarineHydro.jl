using PyCall, MarineHydro, Zygote, Plots, ColorTypes

orange = RGB(230/255,159/255,0/255)  
vermillion = RGB(213/255, 94/255, 0/255) 
bluishgreen = RGB(0/255, 158/255, 115/255) 

include("BemProgram.jl")
include("MeshGradients_singlebody.jl")

cpt = pyimport("capytaine")
radius = 1.0
resolution = (8, 8) 
cptmesh = cpt.mesh_sphere(name="sphere", radius=radius, center=(0, 0, 0), resolution=resolution)
cptmesh.keep_immersed_part(inplace=true)
mesh = Mesh(cptmesh)
omega = 1.03
dof = [0.0,0.0,1.0]



function added_mass_program(radius,omega,dof)  
    mesh = differentiableMesh(radius) #fd
    A = calculate_radiation_forces(mesh,dof,omega)[1]
    return A
end


green_functionsWu = (
    Rankine(),
    RankineReflected(),
    GFWu(),
)

green_functionsDel = (
    Rankine(),
    RankineReflected(),
    ExactGuevelDelhommeau(),
)

function compute_gradients(radius,dof,k_range)
    A(w) = added_mass_program(radius,w,dof)
    gradients = Float64[]
    amass = Float64[]
    for k in k_range
        omega = sqrt(k * 9.8)
        grad, = Zygote.gradient(A, omega)
        am = added_mass_program(radius,omega,dof)
        push!(gradients, grad)
        push!(amass, am)
    end
   return amass, gradients
end


k_range = collect(range(0.01, stop=10, step=0.4))
AmassWu, A_w_gradwu = compute_gradients(radius,dof,k_range)
AmassDel, A_w_graddel = compute_gradients(radius,dof,k_range)


plot(k_range, AmassWu, label="Added mass [Wu]", lw=2, linestyle=:dash, color=:vermillion)
plot!(k_range, AmassDel, label="Added mass [Delhommeau]", lw=2, linestyle=:solid, color=:bluishgreen)
xlabel!("k")
ylabel!("Gradient")
title!("Comparison of Gradients: Wu vs Delhommeau")
legend()
savefig("/home/cornell/BEMJulia/MarineHydro.jl/paper/Plots/added_mass_del_wu.pdf")

plot(k_range, A_w_gradwu, label="A_w_gradwu", lw=2, linestyle=:dash, color=vermillion)
plot!(k_range, A_w_graddel, label="A_w_graddel", lw=2, linestyle=:solid, color=bluishgreen)
xlabel!("k")
ylabel!("Gradient")
title!("Comparison of Gradients: Wu vs Delhommeau")
legend()
savefig("/home/cornell/BEMJulia/MarineHydro.jl/paper/Plots/gradw_added_mass_del_wu.pdf")
