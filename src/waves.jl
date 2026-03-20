function integrate_pressure(floatingbody::FloatingBody, pressure)
    mesh = floatingbody.mesh

    # generator
    force_pairs = (dof_symbol => begin
        dof_mat = floatingbody.dofs[dof_symbol]
        normal_dof_amp_on_face = -sum(dof_mat .* mesh.normals, dims=2)
        sum(pressure .* normal_dof_amp_on_face .* mesh.areas) # output
    end for dof_symbol in keys(floatingbody.dofs))

    # convert Pair to NamedTuple using ; and splat
    forces = (; force_pairs...)
    return forces
end

################################ Radiation methods #########################################

function radiation_bc(mesh::Mesh, dof_mat::Matrix{Float64}, omega::Real)
    """
    Calculates the radiation boundary conditions for floating bodies at each panel.

    # Arguments
    - `floatingbody::FloatingBody`: The floating body
    - `omega`: The frequency of the incident ocean wave ~~~.

    # Returns
    - The (Neumann) radiation boundary condition values for each panel.
"""
    bc =  -1im .* omega .* sum(mesh.normals .* dof_mat, dims=2)
    return bc
end

# Version of compute_bc for a radiation problem
function compute_bc(problem::RadiationProblem)
    return radiation_bc(problem.floatingbody.mesh,
    problem.floatingbody.dofs[problem.radiating_dof],
    problem.omega)
end

################################ Diffraction and Excitation methods #########################################
function airy_waves_potential(points, omega, beta=0)
    wavenumber = omega^2/SETTINGS.g
    x, y, z = points[:, 1], points[:, 2], points[:, 3]
    wbar = x .* cos(beta) .+ y .* sin(beta)
    cih = exp.(wavenumber .* z)
    phi = -1im*SETTINGS.g/omega .* cih .* exp.(1im * wavenumber * wbar)
    return phi
end

function airy_waves_velocity(points, omega, beta=0, water_depth = Inf)
    """Compute the fluid velocity for Airy waves at a given point (or array of points)."""
    k = omega^2/SETTINGS.g

    x, y, z = points[:, 1], points[:, 2], points[:, 3]

    wbar = x .* cos(beta) .+ y .* sin(beta)
    cih = exp.(k .* z)
    sih = exp.(k .* z)

    v = SETTINGS.g * k / omega .* exp.(1im * k .* wbar) .* 
        hcat(cos(beta) .* cih, sin(beta) .* cih, -1im .* sih)
    return v
end


#boundary conditions from airy wave for solving diffraction problem
function AiryBC(mesh,omega,beta=0)
    """Boundary condition for diffraction problem : the velocity on the floating body is the velocity of Airy wave field."""
    bcs = -sum(airy_waves_velocity(mesh.centers,omega,beta) .* mesh.normals, dims = 2)
    return bcs
end

# Version of compute_bc for a diffraction problem
function compute_bc(problem::DiffractionProblem)
    return AiryBC(problem.floatingbody.mesh, problem.omega, problem.beta)
end


function airy_waves_pressure(points, omega, beta=0)
    """Compute the pressure for Airy waves."""  

    return 1im .* omega .* SETTINGS.rho .* airy_waves_potential(points, omega, beta)
end


function FroudeKrylovForce(floatingbody::FloatingBody, ω, beta=0)
    """Compute the Froude-Krylov force."""

    mesh = floatingbody.mesh
    pressure =  airy_waves_pressure(mesh.centers,  ω, beta)
    forces = integrate_pressure(floatingbody, pressure) 
    return forces 
end
