using MarineHydro
using PyCall
using Test
using DimensionalData
using Plots
using ColorTypes
cpt = pyimport("capytaine")



# omegas = range(1.5, 2.25, length=50)
omegas = range(1.5, 2.25, length=15)
betas = [0.0]
DOFs = ["Surge"]

radius = 12.5
draught = 37.5

# cpt mesh and lid
cpt_mesh = cpt.mesh_vertical_cylinder(
    length=draught*2, radius=radius, center=(0, 0, 0), faces_max_radius=3.5,
    ).immersed_part()
cpt_lid_mesh = cpt_mesh.generate_lid(z=cpt_mesh.lowest_lid_position(omega_max=maximum(omegas)))

# mh mesh and lid
mh_mesh = Mesh(cpt_mesh)
mh_lid_mesh = Mesh(cpt_lid_mesh)

# cpt body 1 (no lid)
cpt_body_1 = cpt.FloatingBody(mesh=cpt_mesh)
cpt_body_1.add_translation_dof(name=DOFs[1])

# cpt body 2 (with lid)
cpt_body_2 = cpt.FloatingBody(mesh=cpt_mesh, lid_mesh=cpt_lid_mesh)
cpt_body_2.add_translation_dof(name=DOFs[1])

# mh body 1 (no lid)
mh_body_1 = FloatingBody(mh_mesh, DOFs, "mh_body_1")

# mh body 2 (with lid)
mh_body_2 = FloatingBody(mh_mesh, mh_lid_mesh, DOFs, "mh_body_2")

# Solve settings
method = "direct"
if method == "direct"
    direct = true
else
    direct = false
end
gf = "ExactGuevelDelhommeau"

# Solve using cpt 
solver = cpt.BEMSolver()
xr = pyimport("xarray")
test_matrix = xr.Dataset(coords=Dict("omega" => omegas, "wave_direction" => betas, "radiating_dof" => DOFs))

cpt_results_1 = cpt.BEMSolver().fill_dataset(test_matrix, cpt_body_1, method="direct")  
cpt_results_2 = cpt.BEMSolver().fill_dataset(test_matrix, cpt_body_2, method="direct") 

# Solve using mh
parameters = (wave_frequencies=omegas, 
        wave_directions=betas,
        radiating_dofs=Symbol.(DOFs),
        influenced_dofs=Symbol.(DOFs))

mh_results_1 = compute_and_label_hydrodynamic_coefficients(parameters, mh_body_1)
mh_results_2 = compute_and_label_hydrodynamic_coefficients(parameters, mh_body_2)

cpt_group = (cpt_results_1, cpt_results_2)
mh_group  = (mh_results_1, mh_results_2)

cpt_results = cpt_results_2
mh_results = mh_results_2

# orange = RGB(230/255,159/255,0/255)  
# vermillion = RGB(213/255, 94/255, 0/255) 
# bluishgreen = RGB(0/255, 158/255, 115/255) 


# plot(omegas, abs.(vec(cpt_results_1.excitation_force.values)), 
#      xlabel = "Frequency (rad/s)", ylabel = "Excitation Force Magnitude [N]",
#      label = "CPT without lid", marker = :square, markersize = 6, linewidth = 6,
#      linecolor = bluishgreen, markercolor = bluishgreen, linestyle = :dot)
# plot!(omegas, abs.(vec(cpt_results_2.excitation_force.values)), 
#       label = "CPT with lid", marker = :utriangle, markersize = 6, linewidth = 6,
#       linecolor = bluishgreen, markercolor = bluishgreen)
# plot!(omegas, abs.(vec(mh_results_1.excitation_force.data)), 
#       label = "MH without lid", marker = :diamond, markersize = 5, linewidth = 3,
#       linecolor = vermillion, markercolor = vermillion, linestyle = :dot)
# plot!(omegas, abs.(vec(mh_results_2.excitation_force.data)), 
#       label = "MH with lid", marker = :dtriangle, markersize = 5, linewidth = 3,
#       linecolor = vermillion, markercolor = vermillion)
# savefig("Excitation_force_for_lid_method.png")

# plot(omegas, abs.(vec(cpt_results_1.added_mass.values)), 
#      xlabel = "Frequency (rad/s)", ylabel = "Added Mass [kg]",
#      label = "CPT without lid", marker = :square, markersize = 6, linewidth = 6,
#      linecolor = bluishgreen, markercolor = bluishgreen, linestyle = :dot)
# plot!(omegas, abs.(vec(cpt_results_2.added_mass.values)), 
#       label = "CPT with lid", marker = :utriangle, markersize = 6, linewidth = 6,
#       linecolor = bluishgreen, markercolor = bluishgreen)
# plot!(omegas, abs.(vec(mh_results_1.added_mass.data)), 
#       label = "MH without lid", marker = :diamond, markersize = 5, linewidth = 3,
#       linecolor = vermillion, markercolor = vermillion, linestyle = :dot)
# plot!(omegas, abs.(vec(mh_results_2.added_mass.data)), 
#       label = "MH with lid", marker = :dtriangle, markersize = 5, linewidth = 3,
#       linecolor = vermillion, markercolor = vermillion)
# savefig("Added_mass_for_lid_method.png")

# plot(omegas, abs.(vec(cpt_results_1.radiation_damping.values)), 
#      xlabel = "Frequency (rad/s)", ylabel = "Radiation Damping [N s / m]",
#      label = "CPT without lid", marker = :square, markersize = 6, linewidth = 6,
#      linecolor = bluishgreen, markercolor = bluishgreen, linestyle = :dot)
# plot!(omegas, abs.(vec(cpt_results_2.radiation_damping.values)), 
#       label = "CPT with lid", marker = :utriangle, markersize = 6, linewidth = 6,
#       linecolor = bluishgreen, markercolor = bluishgreen)
# plot!(omegas, abs.(vec(mh_results_1.radiation_damping.data)), 
#       label = "MH without lid", marker = :diamond, markersize = 5, linewidth = 3,
#       linecolor = vermillion, markercolor = vermillion, linestyle = :dot)
# plot!(omegas, abs.(vec(mh_results_2.radiation_damping.data)), 
#       label = "MH with lid", marker = :dtriangle, markersize = 5, linewidth = 3,
#       linecolor = vermillion, markercolor = vermillion)
# savefig("Radiation_damping_for_lid_method.png")

@testset "Comparison of hydrodynamic coefficients with Capytaine using lid method" begin
    # Get Capytaine outputs
    A_cpt = cpt_results.added_mass
    B_cpt = cpt_results.radiation_damping
    F_FK_cpt = cpt_results.Froude_Krylov_force 
    F_D_cpt = cpt_results.diffraction_force
    F_ex_cpt = cpt_results.excitation_force

    # Get MarineHydro outputs
    A_mh = mh_results.added_mass
    B_mh = mh_results.radiation_damping
    F_FK_mh = mh_results.Froude_Krylov_force
    F_D_mh = mh_results.diffraction_force
    F_ex_mh = mh_results.excitation_force

    for omega in omegas
        for influenced_dof in DOFs
            for radiating_dof in DOFs
                @testset "Omega: $omega, influenced_dof: $influenced_dof, radiating_dof: $radiating_dof" begin
                    # Test added mass
                    a_cpt = A_cpt.sel(omega=omega, radiating_dof=radiating_dof, influenced_dof=influenced_dof).values[]
                    a_mh = A_mh[influenced_dofs = At(Symbol(influenced_dof)),radiating_dofs = At(Symbol(radiating_dof)), wave_frequencies = At(omega), forward_speeds = At(0)]
                    @test  a_cpt ≈ a_mh atol=1e-4 rtol = 3e-1
                    # Test radiation damping
                    b_cpt = B_cpt.sel(omega=omega, radiating_dof=radiating_dof, influenced_dof=influenced_dof).values[]
                    b_mh = B_mh[influenced_dofs = At(Symbol(influenced_dof)),radiating_dofs = At(Symbol(radiating_dof)), wave_frequencies = At(omega), forward_speeds = At(0)]
                    @test  b_cpt ≈ b_mh atol=1e-4 rtol = 3e-1
                end                          
            end
            for beta in betas
                @testset "Omega: $omega, influenced_dof: $influenced_dof, beta: $beta" begin
                    # Test FK force
                    f_FK_cpt = F_FK_cpt.sel(omega=omega, influenced_dof=influenced_dof, wave_direction=beta).values[]
                    f_FK_mh = F_FK_mh[influenced_dofs = At(Symbol(influenced_dof)), wave_frequencies = At(omega), wave_directions = At(beta), forward_speeds = At(0)]
                    @test real(f_FK_cpt) ≈ real(f_FK_mh) atol=1e-4 rtol = 3e-1
                    @test imag(f_FK_cpt) ≈ imag(f_FK_mh) atol=1e-4 rtol = 3e-1
                    # Test diffraction force
                    f_D_cpt = F_D_cpt.sel(omega=omega, influenced_dof=influenced_dof, wave_direction=beta).values[]
                    f_D_mh = F_D_mh[influenced_dofs = At(Symbol(influenced_dof)), wave_frequencies = At(omega), wave_directions = At(beta), forward_speeds = At(0)]
                    @test real(f_D_cpt) ≈ real(f_D_mh) atol=1e-4 rtol = 3e-1
                    @test imag(f_D_cpt) ≈ imag(f_D_mh) atol=1e-4 rtol = 3e-1
                    # Test excitation force
                    f_ex_cpt = F_ex_cpt.sel(omega=omega, influenced_dof=influenced_dof, wave_direction=beta).values[]
                    f_ex_mh = F_ex_mh[influenced_dofs = At(Symbol(influenced_dof)), wave_frequencies = At(omega), wave_directions = At(beta), forward_speeds = At(0)]
                    @test real(f_ex_cpt) ≈ real(f_ex_mh) atol=1e-4 rtol = 3e-1
                    @test imag(f_ex_cpt) ≈ imag(f_ex_mh) atol=1e-4 rtol = 3e-1
                end 
            end           
        end        
    end
end

