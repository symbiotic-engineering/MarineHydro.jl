
module BEM

export greens, gradient_greens, integral, integral_gradient, solve,calculate_radiation_forces
export Mesh, element, integrate_pressure, DiffractionForce,assemble_matrix_wu
export FroudeKrylovForce , AiryBC , airy_waves_pressure ,airy_waves_velocity,airy_waves_potential
export assemble_matrices, radiation_bc, solve_hydro_coefficients 
export GreensFunction, Rankine, RankineReflected, RankineReflectedNegative, WaveComponent
export GFWu, ExactGuevelDelhommeau


using Base: Base  # getindex
using LinearAlgebra: cross, diagm, dot, norm
using SpecialFunctions: besselj0, besselj1, bessely0, bessely1
using LinearAlgebra
using .MathConstants: γ
using ImplicitAD


abstract type GreensFunction end

function greens end
function gradient_greens end
function integral end
function integral_gradient end

const ∑ = sum
const τ̅ = 2π
#τ = 2.0unit_π
const PointNonDim = AbstractVector{<:Real}

include("Meshes.jl") #capytaine one
include("green_functions/rankine.jl")
include("green_functions/rankine_reflected.jl")
include("green_functions/rankine_reflected_negative.jl")
include("green_functions/struve.jl")
include("green_functions/wu_.jl")
include("green_functions/exact_Guevel_Delhommeau.jl")



function assemble_matrices(green_functions::NTuple{N, GreensFunction} where N, mesh::Mesh, wavenumber; direct::Bool=true)
    """
    assemble_matrices(green_functions::NTuple{N, GreensFunction} where N, mesh::Mesh, wavenumber; direct::Bool=true)

        Assembles the influence matrices based on the tuple of provided Green's functions, mesh, and wavenumber.

        # Arguments
        - `green_functions::NTuple{N, GreensFunction} where N`: A tuple of Green's functions (rankine, reflecteRankine,Wave). Look for the methods they need to have.
        - `mesh::Mesh`: Floating body mesh with panel information such as vertices, faces, normals, areas etc.
        - `wavenumber`: Incoming ocean wavenumber
        - `direct::Bool=true`: A flag to specify whether to use direct BEM vs Indirect BEM.

        # Returns
        - A tuple of assembled matrices. S and (D or K) depending on the flag.
    """
    # Use comprehensions to build S and D matrices
    S = @inbounds [sum(-1/2τ̅ * Complex(integral(gf, element(mesh, i), element(mesh, j), wavenumber)) for gf in green_functions) for i in 1:mesh.nfaces, j in 1:mesh.nfaces]
    
    D = @inbounds [begin
            element_i = element(mesh, i)
            element_j = element(mesh, j)

            # Select the normal based on direct flag
            normal = direct ? element_j.normal : element_i.normal

            sum(-1/2τ̅ * Complex(normal' * integral_gradient(gf, element_i, element_j, wavenumber; with_respect_to_first_variable=!direct)) for gf in green_functions)
        end for i in 1:mesh.nfaces, j in 1:mesh.nfaces]

    # Add diagonal elements to D
    @inbounds begin
        @views begin
            D1 = D .+ Diagonal(0.5 .* I(mesh.nfaces))
        end
    end

    return S, D1
end

function assemble_matrix_wu(mesh::Mesh, wavenumber; direct::Bool=true)
    """
    assemble_matrices(green_functions::NTuple{N, GreensFunction} where N, mesh::Mesh, wavenumber; direct::Bool=true)

        Assembles the influence matrices based on the tuple of provided Green's functions, mesh, and wavenumber.

        # Arguments
        - `green_functions::NTuple{N, GreensFunction} where N`: A tuple of Green's functions (rankine, reflecteRankine,Wave). Look for the methods they need to have.
        - `mesh::Mesh`: Floating body mesh with panel information such as vertices, faces, normals, areas etc.
        - `wavenumber`: Incoming ocean wavenumber
        - `direct::Bool=true`: A flag to specify whether to use direct BEM vs Indirect BEM.

        # Returns
        - A tuple of assembled matrices. S and (D or K) depending on the flag.
    """
    # Use comprehensions to build S and D matrices
    green_functions = [Rankine(), RankineReflected(), GFWu()]
    S = @inbounds [sum(-1/2τ̅ * Complex(integral(gf, element(mesh, i), element(mesh, j), wavenumber)) for gf in green_functions) for i in 1:mesh.nfaces, j in 1:mesh.nfaces]
    
    D = @inbounds [begin
            element_i = element(mesh, i)
            element_j = element(mesh, j)

            # Select the normal based on direct flag
            normal = direct ? element_j.normal : element_i.normal

            sum(-1/2τ̅ * Complex(normal' * integral_gradient(gf, element_i, element_j, wavenumber; with_respect_to_first_variable=!direct)) for gf in green_functions)
        end for i in 1:mesh.nfaces, j in 1:mesh.nfaces]

    # Add diagonal elements to D
    @inbounds begin
        @views begin
            D1 = D .+ Diagonal(0.5 .* I(mesh.nfaces))
        end
    end

    return S, D1
end



function radiation_bc(mesh::Mesh, dof, omega)
    """
        radiation_bc(mesh::Mesh, dof, omega)

    Calculates the radiation boundary conditions for floating bodies at each panel.

    # Arguments
    - `mesh::Mesh`: The mesh of the floating body.
    - `dof`: The degrees of freedom (assumed same for each panel).
    - `omega`: The frequency of the incident ocean wave ~~~.

    # Returns
    - The (Neumann) radiation boundary condition values for each panel.
"""
    return -1im .* omega .* sum(mesh.normals .* dof', dims=2)
end

function solve(D, S, bc; direct::Bool=true)
    if direct
        ϕ = implicit_linear(D,S*bc)
    else
        K = D
        sources = implicit_linear(K,bc)
        ϕ = S * sources
    end
    return ϕ
end

function integrate_pressure(mesh::Mesh, pressure, dof)
    normal_dof_amp = -sum(transpose(dof) .* mesh.normals, dims=2)
    forces = sum(pressure .* normal_dof_amp .* mesh.areas)
  return forces
  end


  
function calculate_radiation_forces(mesh::Mesh, dof, omega)
    k = omega^2 / 9.8
    S, D = assemble_matrix_wu(mesh, k)
    bc = radiation_bc(mesh, dof, omega)
    potential = solve(D, S, bc)
    pressure = 1im * 1023 * omega * potential
    forces = integrate_pressure(mesh, pressure, dof)
    return [real(forces)/omega^2, imag(forces)/omega]
end

################################ Diffraction and Excitation methods #########################################
function airy_waves_potential(points, omega )
    g = 9.81
    wavenumber = omega^2/g
    x, y, z = points[:, 1], points[:, 2], points[:, 3]
    beta = 0.0
    wbar = x .* cos(beta) .+ y .* sin(beta)
    cih = exp.(wavenumber .* z)
    phi = -1im*g/omega .* cih .* exp.(1im * wavenumber * wbar)
    return phi
end

function airy_waves_velocity(points, omega, water_depth = Inf, beta = 0)
    """Compute the fluid velocity for Airy waves at a given point (or array of points)."""
    g = 9.8
    k = omega^2/g

    x, y, z = points[:, 1], points[:, 2], points[:, 3]

    wbar = x .* cos(beta) .+ y .* sin(beta)
    cih = exp.(k .* z)
    sih = exp.(k .* z)

    v = g * k / omega .* exp.(1im * k .* wbar) .* 
        hcat(cos(beta) .* cih, sin(beta) .* cih, -1im .* sih)
    return v
end


#boundary conditions from airy wave for solving diffraction problem

function AiryBC(mesh,omega)
    """Boundary condition for diffraction problem : the velocity on the floating body is the velocity of Airy wave field."""
    bcs = -sum(airy_waves_velocity(mesh.centers,omega) .* mesh.normals, dims = 2)
return bcs
end



function airy_waves_pressure(points, omega)
    """Compute the pressure for Airy waves."""
    rho = 1000
    return 1im .* omega .* rho .* airy_waves_potential(points, omega)
end



function FroudeKrylovForce(mesh::Mesh, ω,dof)
    """Compute the Froude-Krylov force."""
    pressure =  airy_waves_pressure(mesh.centers,  ω)
    return  integrate_pressure(mesh::Mesh, pressure, dof) 
end

function DiffractionForce(mesh::Mesh,ω,dof)
    green_functions = (
        Rankine(),
        RankineReflected(),
        GFWu(),
    )
    k = ω^2 / 9.8
    S, D = assemble_matrices(green_functions, mesh, k)
    bc = AiryBC(mesh, ω)
    potential = solve(D, S, bc)
    forces = diffraction_force(potential,mesh, ω,dof)
    return forces
end

# F1 = -sum(dPressure .* sphere_1_heave_normal .* mesh.areas)
# normal_dof_amp = -sum(transpose(dof) .* mesh.normals, dims=2)
#forces = sum(pressure .* normal_dof_amp .* mesh.areas)

function diffraction_force(potential,mesh, omega,dof)
    rho = 1000
    pressure = 1im*rho* potential * omega 
    forces = integrate_pressure(mesh,pressure,dof) 
    return forces  
  end

end
