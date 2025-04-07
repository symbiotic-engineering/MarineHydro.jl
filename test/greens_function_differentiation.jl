
using Test
using Zygote
using MarineHydro
using PyCall
using LinearAlgebra

@testset "Greens Function Differentiability Tests" begin
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
