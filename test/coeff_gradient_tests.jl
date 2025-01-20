using Test
using Zygote
using BEM
using PyCall
using LinearAlgebra
using ImplicitAD
using ChainRulesCore


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