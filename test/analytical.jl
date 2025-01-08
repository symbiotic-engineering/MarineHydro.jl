using BEM
using Test
using PyCall


cpt = pyimport("capytaine")
radius = 1.0 #fixed
resolution = (10, 10)
cptmesh = cpt.mesh_sphere(name="sphere", radius=radius, center=(0, 0, 0), resolution=resolution) 
cptmesh.keep_immersed_part(inplace=true)
mesh = Mesh(cptmesh)

@testset "Comparison with analytical added mass and radiation damping using exact Delhommeau  rtol=1e-1" begin
    K_heave = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]

    A_heave = [0.8764, 0.8627, 0.7938, 0.7157, 0.6452, 0.5861, 0.5381, 0.4999, 0.4698, 0.4464, 0.4284, 0.4047, 0.3924, 0.3871, 0.3864, 0.3884, 0.3988, 0.4111, 0.4322, 0.4471, 0.4574, 0.4647, 0.4700, 0.4740, 0.4771]

    B_heave = [0.1036, 0.1816, 0.2793, 0.3254, 0.3410, 0.3391, 0.3271, 0.3098, 0.2899, 0.2691, 0.2484, 0.2096, 0.1756, 0.1469, 0.1229, 0.1031, 0.0674, 0.0452, 0.0219, 0.0116, 0.0066, 0.0040, 0.0026, 0.0017, 0.0012]

    rho = 1000
    ω = sqrt.(K_heave .* 9.8)
    heave = [0.0, 0.0, 1.0]
    results_vec = [calculate_radiation_forces(mesh,heave,omega) for omega in ω]
    A_bem = [results[1]./((2/3)*pi*radius*rho)  for results in results_vec]
    @test A_bem ≈ A_heave  rtol=1e-1
    B_bem = [results[2]./((2/3)*pi*radius*rho)  for results in results_vec]
    ##dimensional for heave and damping(requires /wn)
    B_bem = B_bem./ω
    @test B_bem ≈ B_heave  rtol=1e-1
end


@testset "Comparison with analytical added mass and radiation damping using Wu gf  rtol=1e-1" begin
    #doi:10.1017/S0022112082001980
    K_heave = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]

    A_heave = [0.8764, 0.8627, 0.7938, 0.7157, 0.6452, 0.5861, 0.5381, 0.4999, 0.4698, 0.4464, 0.4284, 0.4047, 0.3924, 0.3871, 0.3864, 0.3884, 0.3988, 0.4111, 0.4322, 0.4471, 0.4574, 0.4647, 0.4700, 0.4740, 0.4771]

    B_heave = [0.1036, 0.1816, 0.2793, 0.3254, 0.3410, 0.3391, 0.3271, 0.3098, 0.2899, 0.2691, 0.2484, 0.2096, 0.1756, 0.1469, 0.1229, 0.1031, 0.0674, 0.0452, 0.0219, 0.0116, 0.0066, 0.0040, 0.0026, 0.0017, 0.0012]

    rho = 1025
    ω = sqrt.(K_heave .* 9.8)
    heave = [0.0, 0.0, 1.0]
    results_vec = [calculate_radiation_forces(mesh,heave,omega) for omega in ω]
    A_bem = [results[1]./((2/3)*pi*radius*rho)  for results in results_vec]
    @test A_bem ≈ A_heave  rtol=1e-1
    B_bem = [results[2]./((2/3)*pi*radius*rho)  for results in results_vec]
    ##dimensional for heave and damping(requires /wn)
    B_bem = B_bem./ω
    @test B_bem ≈ B_heave  rtol=1e-1

end

@testset "Comparison of Excitation Forces with Capytaine (rtol = 1e-1) " begin

   ### diffraction forces analytical by https://www.sciencedirect.com/science/article/pii/S002980182202813X
    #buut that did not match 
    #from capytaine 
    cpt = pyimport("capytaine")
    r = 1.0
    # Create the mesh and body
    cptmesh = cpt.mesh_sphere(radius=r, center=(0, 0, 0), resolution=(14, 14)).immersed_part()
    cptbody = cpt.FloatingBody(cptmesh, name="sphere")
    cptbody.add_translation_dof(name="Heave")

    # Define the wave numbers and corresponding frequencies
    K_heave_diff = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    omegas = sqrt.(K_heave_diff .* 9.81)
    non_dimensional_const = 2/3 * pi * r^3 * 1000

    # Create the test matrix
    xr = pyimport("xarray")
    test_matrix = xr.Dataset(coords=Dict("omega" => omegas, "wave_direction" => [0.0]))
    results = cpt.BEMSolver().fill_dataset(test_matrix, cptbody, method="direct")
    Froude_heave =  vec(results["Froude_Krylov_force"].values)./ non_dimensional_const
    Diff_heave =   vec(results["diffraction_force"].values) ./ non_dimensional_const
    mesh = Mesh(cptmesh)
    g = 9.81
    dof = [0,0,1]
    julia_diff_omega = [DiffractionForce(mesh,w,dof) for w in omegas] ./ non_dimensional_const
    julia_fr_force = [FroudeKrylovForce(mesh,w,dof) for w in omegas] ./ non_dimensional_const
    @test real.(julia_diff_omega) ≈ real.(Diff_heave) rtol=1e-1
    @test imag.(julia_diff_omega) ≈ imag.(Diff_heave) rtol=1e-1
    @test real.(julia_fr_force) ≈ real.(Froude_heave) rtol=1e-1
    @test imag.(julia_fr_force) ≈ imag.(Froude_heave) atol=1e-1

end 