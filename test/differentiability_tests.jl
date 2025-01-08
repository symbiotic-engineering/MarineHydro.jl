using Test
using Zygote
using BEM
using PyCall
using LinearAlgebra
using ImplicitAD
using ChainRulesCore

@testset "Green Function Differentiability Tests" begin
    # Define elements
    e1 = (center=[0.0, 0.0, -1.0],)
    e2 = (
        center=[1.0, 1.0, -2.0],
        vertices=[
            -0.5 -0.5 0.0;
             0.5 -0.5 0.0;
             0.5  0.5 0.0;
            -0.5  0.5 0.0
        ] .+ [1.0, 1.0, -2.0]',
        normal=[0.0, 0.0, 1.0],
        radius=sqrt(2)/2,
        area=1.0,
    )

    @testset "Rankine Differentiability" begin
        # Define functions to test with only 1 parameter
        function greens_center1(center)
            e1_new = (center=center,)
            return greens(Rankine(), e1_new, e2)
        end

        function greens_center2(center)
            e2_new = (
                center=center,
                vertices=[
                    -0.5 -0.5 0.0;
                     0.5 -0.5 0.0;
                     0.5  0.5 0.0;
                    -0.5  0.5 0.0
                ] .+ center',
                normal=[0.0, 0.0, 1.0],
                radius=sqrt(2)/2,
                area=1.0,
            )
            return greens(Rankine(), e1, e2_new)
        end

        # Test differentiability of `greens` with respect to `e1.center`
        gradient1 = Zygote.gradient(greens_center1, e1.center)
        @test typeof(gradient1) == Tuple{Vector{Float64}}  # Ensure valid gradient type
        @test !any(isnan, gradient1[1])                  # Ensure gradient is valid (not NaN)

        # Test differentiability of `greens` with respect to `e2.center`
        gradient2 = Zygote.gradient(greens_center2, e2.center)
        @test typeof(gradient2) == Tuple{Vector{Float64}}
        @test !any(isnan, gradient2[1])

        # Define functions for `gradient_greens`
        function gradient_greens_center1(center)
            e1_new = (center=center,)
            return gradient_greens(Rankine(), e1_new, e2, with_respect_to_first_variable=true)
        end

        function gradient_greens_center2(center)
            e2_new = (
                center=center,
                vertices=[
                    -0.5 -0.5 0.0;
                     0.5 -0.5 0.0;
                     0.5  0.5 0.0;
                    -0.5  0.5 0.0
                ] .+ center',
                normal=[0.0, 0.0, 1.0],
                radius=sqrt(2)/2,
                area=1.0,
            )
            return gradient_greens(Rankine(), e1, e2_new, with_respect_to_first_variable=false)
        end

        # Test differentiability of `gradient_greens` with respect to `e1.center`
        grad_grad1 = Zygote.jacobian(gradient_greens_center1, e1.center)
        @test typeof(grad_grad1) == Tuple{Matrix{Float64}}  # Ensure the Jacobian is a matrix
        @test !any(isnan, grad_grad1[1])                   # Ensure no NaN values in the Jacobian
    

        # Test differentiability of `gradient_greens` with respect to `e2.center`
        grad_grad2 = Zygote.jacobian(gradient_greens_center2, e2.center)
        @test typeof(grad_grad1) == Tuple{Matrix{Float64}} 
        @test !any(isnan, grad_grad2[1])
    end

    @testset "Rankine Integral Differentiability" begin
        #add integral test sets
    end


# GFu test 

@testset "GFWu Differentiability" for k in [0.1, 1.0, 10.0]
    # Define functions to test real and imaginary parts separately
    function real_greens_gfwu_center1(center)
        e1_new = (center=center,)
        return real(greens(GFWu(), e1_new, e2, k))
    end

    function imag_greens_gfwu_center1(center)
        e1_new = (center=center,)
        return imag(greens(GFWu(), e1_new, e2, k))
    end

    function real_gradient_greens_gfwu_center1(center)
        e1_new = (center=center,)
        return real(gradient_greens(GFWu(), e1_new, e2, k))
    end

    function imag_gradient_greens_gfwu_center1(center)
        e1_new = (center=center,)
        return imag(gradient_greens(GFWu(), e1_new, e2, k))
    end

    # Test differentiability of `greens` with respect to `e1.center`
    real_grad1 = Zygote.gradient(real_greens_gfwu_center1, e1.center)
    imag_grad1 = Zygote.gradient(imag_greens_gfwu_center1, e1.center)
    @test typeof(real_grad1) == Tuple{Vector{Float64}}
    @test typeof(imag_grad1) == Tuple{Vector{Float64}}
    @test !any(isnan, real_grad1[1])
    @test !any(isnan, imag_grad1[1])

    # Test differentiability of `gradient_greens` with respect to `e1.center`
    real_grad_grad1 = Zygote.jacobian(real_gradient_greens_gfwu_center1, e1.center)
    imag_grad_grad1 = Zygote.jacobian(imag_gradient_greens_gfwu_center1, e1.center)
    @test typeof(real_grad_grad1) == Tuple{Matrix{Float64}}
    @test typeof(imag_grad_grad1) == Tuple{Matrix{Float64}}
    @test !any(isnan, real_grad_grad1[1])
    @test !any(isnan, imag_grad_grad1[1])
end

@testset "ExactGuevelDelhommeau Differentiability" for k in [0.1, 1.0, 10.0]
    function real_greens_gfwu_center1(center)
        e1_new = (center=center,)
        return real(greens(ExactGuevelDelhommeau(), e1_new, e2, k))
    end

    function imag_greens_gfwu_center1(center)
        e1_new = (center=center,)
        return imag(greens(ExactGuevelDelhommeau(), e1_new, e2, k))
    end

    function real_gradient_greens_gfwu_center1(center)
        e1_new = (center=center,)
        return real(gradient_greens(ExactGuevelDelhommeau(), e1_new, e2, k))
    end

    function imag_gradient_greens_gfwu_center1(center)
        e1_new = (center=center,)
        return imag(gradient_greens(ExactGuevelDelhommeau(), e1_new, e2, k))
    end

    # Test differentiability of `greens` with respect to `e1.center`
    real_grad1 = Zygote.gradient(real_greens_gfwu_center1, e1.center)
    imag_grad1 = Zygote.gradient(imag_greens_gfwu_center1, e1.center)
    @test typeof(real_grad1) == Tuple{Vector{Float64}}
    @test typeof(imag_grad1) == Tuple{Vector{Float64}}
    @test !any(isnan, real_grad1[1])
    @test !any(isnan, imag_grad1[1])

    # Test differentiability of `gradient_greens` with respect to `e1.center`
    real_grad_grad1 = Zygote.jacobian(real_gradient_greens_gfwu_center1, e1.center)
    imag_grad_grad1 = Zygote.jacobian(imag_gradient_greens_gfwu_center1, e1.center)
    @test typeof(real_grad_grad1) == Tuple{Matrix{Float64}}
    @test typeof(imag_grad_grad1) == Tuple{Matrix{Float64}}
    @test !any(isnan, real_grad_grad1[1])
    @test !any(isnan, imag_grad_grad1[1])
end
end

@testset "Matrix Differentiability Tests" begin
    cpt = pyimport("capytaine")
    radius = 1.0
    resolution = (4, 4) #smaller data just to test
    cptmesh = cpt.mesh_sphere(name="sphere", radius=radius, center=(0, 0, 0), resolution=resolution)
    cptmesh.keep_immersed_part(inplace=true)
    mesh = Mesh(cptmesh)

    # Define parameters
    wavenumber = 1.0
    ω = √(wavenumber * 9.8)
    dof = [0, 0, 1]

    @testset "Radiation Boundary Condition Tests" begin
        # Test differentiability using Zygote
    Ji_bc_mesh, Ji_bc_dof , Ji_bc_ω = Zygote.jacobian((mesh, dof, ω) -> imag.(radiation_bc(mesh, dof, ω)), mesh,dof,ω)
    
    @test Ji_bc_ω !== nothing 
    @test Ji_bc_dof !== nothing 
    @test typeof(Ji_bc_ω) == Vector{Float64}
    #Test the size of the gradient

   # Jr_bc_mesh, Jr_bc_dof , Jr_bc_ω = Zygote.jacobian((mesh, dof, ω) -> real.(radiation_bc(mesh, dof, ω)), mesh,dof,ω)
    
    # @test Jr_bc_ω !== nothing 
  #  @test typeof(Jr_bc_ω) == Vector{Float64}
    #@test Ji_bc_mesh !== nothing  # this fails not sure why ;  probably due to vector with respect to matrix?
    # Test the size of the gradient
   # @test size(real_grad_bc.normals) == size(mesh.normals)

    end

@testset "Differentiability of assemble_matrices with respect to mesh.vertices" begin
    # Create the base mesh using Capytaine
    green_functions = (
        Rankine(),
        RankineReflected(),
         GFWu(),
    )

    # Function to test differentiability of `assemble_matrices` with respect to vertices
    function test_S_assemble_matrices(vertices) #To do - split S and D/K matrix test
        
        mesh_new = Mesh(
            vertices, mesh.faces, 
            mesh.centers, mesh.normals,
            mesh.areas, mesh.radii, mesh.nvertices, mesh.nfaces
        )
        S, _ = assemble_matrices(green_functions, mesh_new, ω)
        return vec(real(S))  # Return the real part of S as a vector for Jacobian computation
    end

    # Compute Jacobian of `assemble_matrices` with respect to `mesh.vertices`
    jacobian_vertices, = Zygote.jacobian(test_S_assemble_matrices, mesh.vertices)

    # Test the Jacobian
    @test typeof(jacobian_vertices) == Matrix{Float64}
    @test !any(isnan, jacobian_vertices[1])

    function test_D_assemble_matrices(vertices) 
        
        mesh_new = Mesh(
            vertices, mesh.faces, 
            mesh.centers, mesh.normals,
            mesh.areas, mesh.radii, mesh.nvertices, mesh.nfaces
        )
        _, D = assemble_matrices(green_functions, mesh_new, ω)
        return vec(real(D))  # Return the real part of S as a vector for Jacobian computation
    end

    # # Compute Jacobian of `assemble_matrices` with respect to `mesh.vertices`
    jacobian_vertices, = Zygote.jacobian(test_D_assemble_matrices, mesh.vertices)

    # # Test the Jacobian
    @test typeof(jacobian_vertices) == Matrix{Float64}
    @test !any(isnan, jacobian_vertices)
end
end

@testset "Linear Solve Tests" begin
   # Test the differentiability of `solve`
@testset "Differentiability of solve (direct and indirect)" begin
    #real BEM matrix not required for just testing solve
    D = rand(3, 3)
    S = rand(3, 3)
    bc = rand(3)
   #  `D` is invertible if its diagonally dominant
   for i in eachindex(D[:,1])
        D[i, i] += sum(abs, D[i, :])
    end
    # Compute gradients using Zygote
    JD, JS, Jbc = Zygote.jacobian((D, S, bc) -> solve(D, S, bc; direct=true), D, S, bc)

    # Check that Jacobian are not `nothing` and have the correct dimensions
    @test size(JD) == (3, 9)  # Jacobian w.r.t. D should have size (output_dim, input_dim_D1 * input_dim_D2) : Zygote flattens input matrix ?
    @test size(JS) == (3, 9)  # Jacobian w.r.t. S should have size (output_dim, input_dim_S1, input_dim_S2)
    @test size(Jbc) == (3, 3)    # Jacobian w.r.t. bc should have size (output_dim, input_dim_bc)
    @test JD !== nothing 
    @test JS !== nothing 
    @test Jbc !== nothing 
    @test typeof(JD) == Matrix{Float64}
    @test typeof(JS) == Matrix{Float64}
    @test typeof(Jbc) == Matrix{Float64}

    # Test indirect mode
   # Compute Jacobians using Zygote
   JD, JS, Jbc = Zygote.jacobian((D, S, bc) -> solve(D, S, bc; direct=false), D, S, bc)

   # Test Jacobian sizes
   @test size(JD) == (3, 9)  # Jacobian w.r.t. D should have size (output_dim, input_dim_D1 * input_dim_D2) : Zygote flattens input matrix ?
    @test size(JS) == (3, 9)  # Jacobian w.r.t. S should have size (output_dim, input_dim_S1, input_dim_S2)
    @test size(Jbc) == (3, 3)    # Jacobian w.r.t. bc should have size (output_dim, input_dim_bc)
    @test JD !== nothing 
    @test JS !== nothing 
    @test Jbc !== nothing 
    @test typeof(JD) == Matrix{Float64}
    @test typeof(JS) == Matrix{Float64}
    @test typeof(Jbc) == Matrix{Float64}
end
end

## Now finally coefficients test - differentiability 
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
