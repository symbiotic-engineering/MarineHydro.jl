using Plots
using CSV
using DataFrames
using BEM
include("/home/cornell/BEMJulia/BEM.jl/paper/meshGradients_pair.jl")
 #takes a while depending on this - some faces_max_radius gives weird answer
#check meshes.jl to change faces_max_radius.
radius_range = [1,2,3,4,5]
heave = [0,0,1] #heave
surge = [1,0,0]

plot(size=(800, 600))



#impact of one body on the other added mass coeffficients
function added_mass_off_diagonal(radius,omega ,dx1)  
    mesh = differentiableMeshPairs(radius, dx1)  
    total_nfaces = mesh.nfaces 
    element_is_in_sphere_1(j) = get_center_x(j, mesh) <= radius
    element_is_in_sphere_2(j) = get_center_x(j, mesh) > radius
    sphere_1_heave_normal = [element_is_in_sphere_1(j) ? mesh.normals[j,:]' * [0.0, 0.0, 1.0] : 0.0 for j in 1:total_nfaces]
    sphere_2_heave_normal = [element_is_in_sphere_2(j) ? mesh.normals[j,:]' * [0.0, 0.0, 1.0] : 0.0 for j in 1:total_nfaces]
    k = omega^2 / 9.81  # Wave number
    S, D = assemble_matrices((Rankine(), RankineReflected(), GFWu()), mesh, k)
    potential = BEM.solve(D, S, -1im * omega * sphere_1_heave_normal)
    pressure = 1im * 1000 * omega * potential
    # force = -sum(pressure .* sphere_2_heave_normal .* mesh.areas)
    # A12 = real(force) / omega^2
    force = -sum(pressure .* sphere_1_heave_normal .* mesh.areas)
    A11 = real(force) / omega^2
    return A11
end

function damping_off_diagonal(radius,omega ,dx1)  
    mesh = differentiableMeshPairs(radius, dx1)  
    total_nfaces = mesh.nfaces 
    element_is_in_sphere_1(j) = get_center_x(j, mesh) <= radius
    element_is_in_sphere_2(j) = get_center_x(j, mesh) > radius
    sphere_1_heave_normal = [element_is_in_sphere_1(j) ? mesh.normals[j,:]' * [0.0, 0.0, 1.0] : 0.0 for j in 1:total_nfaces]
    sphere_2_heave_normal = [element_is_in_sphere_2(j) ? mesh.normals[j,:]' * [0.0, 0.0, 1.0] : 0.0 for j in 1:total_nfaces]
    k = omega^2 / 9.81  # Wave number
    S, D = assemble_matrices((Rankine(), RankineReflected(), GFWu()), mesh, k) # Assemble matrices tuple error -- use default
    potential = BEM.solve(D, S, -1im * omega * sphere_1_heave_normal)
    pressure = 1im * 1000 * omega * potential
    force = -sum(pressure .* sphere_1_heave_normal .* mesh.areas)
    B11 = imag(force) / omega
    return B11
end


# Set parameters
g = 9.8 
heave = [0, 0, 1]  # Heave
dx_r_ratios = collect(range(1.5, stop=8.0, step=0.5))
kr_values = collect(range(0.5, stop=10.0, step=1.0))  # k*r dimensionless parameter
r = 1.0   
omega = 1.03
#check with heuristics that the  - plot sensitivity of added mass and damping with dx and see if they go to zero or low.
data = DataFrame(A11_grad_r=Float64[], B11_grad_r =Float64[], dx = Float64[])
for dx in dx_r_ratios
    @show dx
    A11_grad_r, = Zygote.gradient(x -> added_mass_off_diagonal(r,omega ,x), dx) #at fix r = 1.0
    B11_grad_r, =  Zygote.gradient(x -> damping_off_diagonal(r,omega ,x), dx)
    push!(data, (A11_grad_r, B11_grad_r, dx))
end

CSV.write("/home/cornell/BEMJulia/BEM.jl/paper/Plots/11_heuristics_dx_DELhommeau.csv", data)


# data = DataFrame(dx_r_ratio=Float64[], kr=Float64[], grad_r = Float64[])

# for dx_r in dx_r_ratios
#     for kr in kr_values
#         dx1 = dx_r * r
#         @show dx1
#         omega = sqrt(kr * g / r)
#         @show omega
#         grad_r, = Zygote.gradient(x -> added_mass_off_diagonal(x, omega, dx1), r)
#         @show grad_r
#         push!(data, (dx_r, kr, grad_r))
#     end
# end

# CSV.write("/home/cornell/BEMJulia/BEM.jl/paper/Plots/added_mass_data_dimensionless.csv", data)

# println("damping data")
# data = DataFrame(dx_r_ratio=Float64[], kr=Float64[], grad_r = Float64[])

# for dx_r in dx_r_ratios
#     for kr in kr_values
#         dx1 = dx_r * r
#         omega = sqrt(kr * g / r)
#         grad_r, = Zygote.gradient(x -> damping_off_diagonal(x,omega ,dx1), r)
#         @show grad_r
#         push!(data, (dx_r, kr, grad_r))
#     end
# end

# CSV.write("/home/cornell/BEMJulia/BEM.jl/paper/Plots/damping_data_dimensionless.csv", data)
