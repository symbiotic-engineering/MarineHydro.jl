using MarineHydro
using Test
using PyCall

cpt = pyimport("capytaine")
radius = 1.0
resolution = (5, 5)
cpt_mesh_sphere = cpt.mesh_sphere(radius=radius, center=(0, 0, 0), resolution=resolution, name="floating sphere").immersed_part()
cpt_mesh_two_spheres = (cpt_mesh_sphere + cpt_mesh_sphere.translated_x(5.0)).copy(name="floating two spheres")

@testset "Comparison with Capytaine for $(cptmesh.name)" for cptmesh in [cpt_mesh_sphere, cpt_mesh_two_spheres]

    mesh = Mesh(cptmesh)

    @testset "Matrix shape" begin
        wavenumber = 1.0
        ω = √(wavenumber*9.8)
        dof = [0,0,1]
        green_functions = (
            Rankine(),
            RankineReflected(),
            #=GFWu()=#
        )
        S, D = assemble_matrices(green_functions, mesh, 1.0);
        S_, K = assemble_matrices(green_functions, mesh, 1.0; direct=false);
        @test size(S) == size(S_) == size(D) == size(K)
        @test S ≈ S_
        @test !(D ≈ K)
    end

    @testset "Comparison with Capytaine" begin
        # As long as we use Capytaine for the mesh, we might as well use it also for validation.
        # When we get rid of the dependency to Capytaine, this test might need to be restructured.

        @testset "Pure Rankine, direct=$direct" for direct in [true, false]
            green_functions = (Rankine(),)
            S, D = assemble_matrices(green_functions, mesh, 1.0, direct=direct)
            capy_S, capy_D = cpt.Delhommeau().evaluate(cptmesh, cptmesh, free_surface=Inf, water_depth=Inf, wavenumber=1.0, adjoint_double_layer=!direct)
            @test S ≈ capy_S  atol=1e-6
            @test D ≈ capy_D  atol=1e-6
        end

        @testset "Rankine + reflected, direct=$direct" for direct in [true, false]
            green_functions = (Rankine(), RankineReflected())
            S, D = assemble_matrices(green_functions, mesh, 1.0, direct=direct)
            capy_S, capy_D = cpt.Delhommeau().evaluate(cptmesh, cptmesh, free_surface=0.0, water_depth=Inf, wavenumber=0.0, adjoint_double_layer=!direct)
            @test S ≈ capy_S  atol=1e-6
            @test D ≈ capy_D  atol=1e-6
        end

        @testset "Rankine - reflected, direct=$direct" for direct in [true, false]
            green_functions = (Rankine(), RankineReflectedNegative())
            S, D = assemble_matrices(green_functions, mesh, 1.0, direct=direct)
            capy_S, capy_D = cpt.Delhommeau().evaluate(cptmesh, cptmesh, free_surface=0.0, water_depth=Inf, wavenumber=Inf, adjoint_double_layer=!direct)
            @test S ≈ capy_S  atol=1e-6
            @test D ≈ capy_D  atol=1e-6
        end

        @testset "Full Green function with Guével-Delhommeau, direct=$direct, k=$k" for direct in [true, false], k in [1.0, 2.0]
            green_functions = (Rankine(), RankineReflected(), ExactGuevelDelhommeau(),)
            S, D = assemble_matrices(green_functions, mesh, k, direct=direct)
            capy_S, capy_D = cpt.Delhommeau(tabulation_nr=0).evaluate(cptmesh, cptmesh, free_surface=0.0, water_depth=Inf, wavenumber=k, adjoint_double_layer=!direct)
            @test real.(S) ≈ real.(capy_S)  atol=1e-3
            @test imag.(S) ≈ imag.(capy_S)  atol=1e-6
            @test real.(D) ≈ real.(capy_D)  atol=1e-4
            @test imag.(D) ≈ imag.(capy_D)  atol=1e-6
        end

        @testset "Full Green function with Wu, direct=$direct, k=$k" for direct in [true, false], k in [1.0, 2.0]
            green_functions = (Rankine(), RankineReflected(), GFWu())
            S, D = assemble_matrices(green_functions, mesh, k, direct=direct)
            capy_S, capy_D = cpt.Delhommeau(tabulation_nr=0).evaluate(cptmesh, cptmesh, free_surface=0.0, water_depth=Inf, wavenumber=k, adjoint_double_layer=!direct)
            @test real.(S) ≈ real.(capy_S)  atol=1e-2
            @test imag.(S) ≈ imag.(capy_S)  atol=1e-6
            @test real.(D) ≈ real.(capy_D)  atol=1e-2
            @test imag.(D) ≈ imag.(capy_D)  atol=1e-6
        end

        @testset "Diffraction and Excitation methods, boundary condtion, airywaves test" begin 
        cpt = pyimport("capytaine")
        r = 1.0
        ω = 1.03
        # Create the mesh and body
        cptmesh = cpt.mesh_sphere(radius=r, center=(0, 0, 0), resolution=(14, 14)).immersed_part()
        cptbody = cpt.FloatingBody(cptmesh, name="sphere")
        cptbody.add_translation_dof(name="Heave")
        pb = cpt.DiffractionProblem(body = cptbody,omega = ω,wave_direction= 0 )
        diffProbBC = pb.boundary_condition
        mesh  = Mesh(cptmesh)
        juliaBC = AiryBC(mesh,ω)
        @test diffProbBC ≈ juliaBC atol=1e-3 rtol = 1e-3

        capyairy = cpt.bem.airy_waves.airy_waves_pressure(cptmesh.faces_centers, pb)
        juliaairy = airy_waves_pressure(mesh.centers, ω)
        @test abs.(capyairy) ≈ abs.(juliaairy) atol=1e-3 rtol = 1e-3

        capyairyV =cpt.bem.airy_waves.airy_waves_velocity(cptmesh.faces_centers,pb)
        juliaairyV = airy_waves_velocity(mesh.centers,ω)
        @test capyairyV ≈ juliaairyV atol=1e-4 rtol = 1e-4

        capyairyPot =cpt.bem.airy_waves.airy_waves_potential(cptmesh.faces_centers,pb)
        juliaairyPot = airy_waves_potential(mesh.centers, ω)
        @test capyairyPot ≈ juliaairyPot atol=1e-4 rtol = 1e-4

        end 

    end

end
