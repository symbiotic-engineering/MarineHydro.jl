

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
