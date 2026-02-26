using Test
using Zygote
using MarineHydro
using PyCall
using LinearAlgebra


@testset "Integrate Pressure " begin
    cpt = pyimport("capytaine")
    radius = 1.0
    resolution = (4, 4) #smaller data just to test
    cptmesh = cpt.mesh_sphere(name="sphere", radius=radius, center=(0, 0, 0), resolution=resolution)
    cptmesh.keep_immersed_part(inplace=true)
    mesh = Mesh(cptmesh)
    @testset "Re(Forces) and Im(Forces) test" begin
        pressure = rand(size(mesh.areas))   # Random pressure values
        dof = [0.0,0.0,1.0]
        # Differentiability check
        real_grad_zygote, = Zygote.gradient(mesh -> real.(integrate_pressure(mesh, pressure, dof)), mesh) #with respect to mesh
        @test real_grad_zygote !== nothing

        imag_grad_zygote, = Zygote.gradient(mesh -> imag.(integrate_pressure(mesh, pressure, dof)), mesh) #with respect to mesh
        @test imag_grad_zygote !== nothing

        grad_areas = real_grad_zygote.areas #integration depends on normal and area only
        @test grad_areas !== nothing
        grad_normals = real_grad_zygote.normals
        @test typeof(grad_normals) == Matrix{Float64}
        @test grad_normals !== nothing

        # Test the gradient shape (should be the same shape as mesh panels)
        @test size(grad_areas) == size(mesh.areas)
        @test size(grad_normals) == size(mesh.normals)
    end
end

@testset "Gradient accuracy check with Finite diff [w.r.t omega]" begin
    using FiniteDifferences
    cpt = pyimport("capytaine")
    radius = 1.0
    resolution = (4, 4) #smaller data just to test
    cptmesh = cpt.mesh_sphere(name="sphere", radius=radius, center=(0, 0, 0), resolution=resolution).keep_immersed_part()
    mesh = Mesh(cptmesh)
    omega = 1.12
    dof = [0.0,0.0,1.0]
    green_functions = (
        Rankine(),
        RankineReflected(),
        GFWu(),
    )

    function check_added_mass(ω,mesh,dof)
        A = calculate_radiation_forces(mesh,dof,ω)[1]
        return A
    end

    @testset "A(w) partials with Wu " begin
        A(w) = check_added_mass(w, mesh, dof)
        A_w_grad, = Zygote.gradient(A,omega)
        #checking accuracy with finite FiniteDifferences
    # Central difference with order 5
        fd_grad = FiniteDifferences.central_fdm(5, 1)(A, omega)
        @test A_w_grad !== nothing
        @test typeof(A_w_grad) == Float64
        @test A_w_grad ≈  fd_grad atol=1e-6 rtol=1e-6
    end

    @testset "A(w) with w using exact Delhommeau" begin
        green_functions = (
            Rankine(),
            RankineReflected(),
            ExactGuevelDelhommeau(),
        )
        A(w) = check_added_mass(w, mesh, dof)
        A_w_grad1, = Zygote.gradient(A,omega)
        #checking accuracy with FiniteDifferences
    # Central difference with order 5
        fd_grad1 = FiniteDifferences.central_fdm(5, 1)(A, omega)
        @test A_w_grad1 !== nothing
        @test typeof(A_w_grad1) == Float64
        @test A_w_grad1 ≈  fd_grad1 atol=1e-6 rtol=1e-6
    end
end

@testset "Gradient accuracy check with Finite diff [w.r.t omega] for MDOF cylnder" begin
    # Description of problem
    omegas = [1.0, 1.5] # frequencies [rad/s]
    beta = 0 # incident wave angle [rad]
    t_DOFs = ["Heave"] # translational DOFs
    r_DOFs = ["Pitch"] # rotational DOFs
    DOFs = [t_DOFs; r_DOFs] # all DOFs

    # Create Mesh object
    radius = 1.5  
    center = (0.0,0.0,0.0) 
    len = 2.5
    faces_max_radius = 0.7
    cptmesh = cpt.meshes.predefined.mesh_horizontal_cylinder(
                radius=radius,
                center=center, 
                length=len, 
                faces_max_radius = faces_max_radius
                ).keep_immersed_part(inplace=true)

    # Get MarineHydro values
    mesh = Mesh(cptmesh)
    rigid_dof_list = DOFs
    rotation_center = [1.0, 1.0, 0.0] # off set for nonzero off-diagoinal elements
    fb = FloatingBody(mesh, rigid_dof_list, rotation_center)

    # Radiation solve functions
    function A_and_B_vec(w)
        added_mass_dict, damping_dict = calculate_radiation_forces(fb, w)
        A_vals = [real(added_mass_dict[(w, i, r)]) for i in DOFs, r in DOFs]
        B_vals = [real(damping_dict[(w, i, r)]) for i in DOFs, r in DOFs]
        return vcat(vec(A_vals), vec(B_vals)) # Note this is [A_11,A_12,A_21...,B_22]
    end
    function Jacobian_of_rad_problem(Omega)
        # This inclusion of imaginary inputs, then taking the real value was required to get jacobian to work
        j = Zygote.jacobian(Omega + 0im) do w
            A_and_B_vec(real(w)) 
        end
        return vec(real.(j)[1])
    end

    # Incident solve function
    function F_FK_vec(w)
        F_FK_dict = FroudeKrylovForce(fb, w)
        F_FK_vals = [F_FK_dict[(w, i)] for i in DOFs]
        return vcat(real.(vec(F_FK_vals)),imag.(vec(F_FK_vals))) # Note this is [real(F_FK_1),real(F_FK_2),...,imag(F_FK_1),imag(F_FK_2),... ]
    end
    function Jacobian_of_inc_problem(Omega)
        # This inclusion of imaginary inputs, then taking the real value was required to get jacobian to work
        j = Zygote.jacobian(Omega + 0im) do w
            F_FK_vec(real(w)) 
        end
        return vec(real.(j)[1])
    end

    # Diffraction solve function
    function F_D_vec(w)
        F_D_dict = DiffractionForce(fb, w)
        F_D_vals = [F_D_dict[(w, i)] for i in DOFs]
        return vcat(real.(vec(F_D_vals)),imag.(vec(F_D_vals))) # Note this is [real(F_D_1),real(F_D_2),...,imag(F_D_1),imag(F_D_2),... ]
    end
    function Jacobian_of_diff_problem(Omega)
        # This inclusion of imaginary inputs, then taking the real value was required to get jacobian to work
        j = Zygote.jacobian(Omega + 0im) do w
            F_D_vec(real(w)) 
        end
        return vec(real.(j)[1])
    end

    for omega in omegas
        @testset "Verify sensitivity of added mass and damping wrt Omega: $omega" begin
            A_and_B_FD = FiniteDifferences.central_fdm(5, 1)(A_and_B_vec, omega)
            A_and_B_AD = Jacobian_of_rad_problem(omega)
            @test A_and_B_AD !== nothing
            @test typeof(A_and_B_AD) == Vector{Float64}
            @test A_and_B_AD ≈ A_and_B_FD atol=1e-6 rtol=1e-6
        end
        @testset "Verify sensitivity of Froude Krylov force wrt Omega: $omega" begin
            F_FK_FD = FiniteDifferences.central_fdm(5, 1)(F_FK_vec, omega)
            F_FK_AD = Jacobian_of_inc_problem(omega)
            @test F_FK_AD !== nothing
            @test typeof(F_FK_AD) == Vector{Float64}
            @test F_FK_AD ≈ F_FK_FD atol=1e-6 rtol=1e-6
        end
        @testset "Verify sensitivity of diffraction force wrt Omega: $omega" begin
            F_D_FD = FiniteDifferences.central_fdm(5, 1)(F_D_vec, omega)
            F_D_AD = Jacobian_of_diff_problem(omega)
            @test F_D_AD !== nothing
            @test typeof(F_D_AD) == Vector{Float64}
            @test F_D_AD ≈ F_D_FD atol=1e-6 rtol=1e-6
        end
    end   

end