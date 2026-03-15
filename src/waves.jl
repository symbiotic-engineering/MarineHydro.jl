# Old version
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


function radiation_bc(floatingbody::FloatingBody, omega)
    """
    Calculates the radiation boundary conditions for floating bodies at each panel.

    # Arguments
    - `floatingbody::FloatingBody`: The floating body
    - `omega`: The frequency of the incident ocean wave ~~~.

    # Returns
    - The (Neumann) radiation boundary condition values for each panel.
"""
    bc = Dict()
    for dof_name in keys(floatingbody.dofs)
        dof_mat = floatingbody.dofs[dof_name]
        normals_mat = floatingbody.mesh.normals
        bc[dof_name] = -1im .* omega .* sum(normals_mat .* dof_mat, dims=2)
    end
    return bc
end


function integrate_pressure(floatingbody::FloatingBody, pressure)
  mesh = floatingbody.mesh
  forces = Dict()
  for dof_name in keys(floatingbody.dofs)
    dof_mat = floatingbody.dofs[dof_name]
    normal_dof_amp_on_face = -sum(dof_mat .* mesh.normals, dims=2)
    forces[dof_name] = sum(pressure .* normal_dof_amp_on_face .* mesh.areas)
  end
  return forces
end

# Old version
function integrate_pressure(mesh::Mesh, pressure, dof)
    normal_dof_amp = -sum(transpose(dof) .* mesh.normals, dims=2)
    forces = sum(pressure .* normal_dof_amp .* mesh.areas)
  return forces
end
  
function calculate_radiation_forces(floatingbody::FloatingBody, omega)
    Added_mass = Dict{Tuple{Any, String, String}, Any}()
    Radiation_damping = Dict{Tuple{Any, String, String}, Any}()
    k = omega^2 / SETTINGS.g
    mesh = floatingbody.mesh
    S, D = assemble_matrix_wu(mesh, k)
    rad_bcs = radiation_bc(floatingbody, omega) 
    for rad_dof in keys(floatingbody.dofs) # radiating dofs
        rad_bc = rad_bcs[rad_dof]
        potential = solve(D, S, rad_bc)
        pressure = 1im * SETTINGS.rho * omega * potential
        forces = integrate_pressure(floatingbody, pressure)
        for inf_dof in keys(forces) # influenced dofs
            force = forces[inf_dof]
            Added_mass[(omega, inf_dof, rad_dof)] = real(force)/omega^2
            Radiation_damping[(omega, inf_dof, rad_dof)] = imag(force)/omega 
        end
    end  
    return Added_mass, Radiation_damping
end

# Old version 
function calculate_radiation_forces(mesh::Mesh, dof, omega)
    k = omega^2 / SETTINGS.g
    S, D = assemble_matrix_wu(mesh, k)
    bc = radiation_bc(mesh, dof, omega)
    potential = solve(D, S, bc)
    pressure = 1im * SETTINGS.rho * omega * potential
    forces = integrate_pressure(mesh, pressure, dof)
    return [real(forces)/omega^2, imag(forces)/omega]
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



function airy_waves_pressure(points, omega, beta=0)
    """Compute the pressure for Airy waves."""
   

    return 1im .* omega .* SETTINGS.rho .* airy_waves_potential(points, omega, beta)
end

function FroudeKrylovForce(floatingbody::FloatingBody, ω, beta=0)
    """Compute the Froude-Krylov force."""
    F_FK = Dict{Tuple{Any, String}, Any}()
    mesh = floatingbody.mesh
    pressure =  airy_waves_pressure(mesh.centers,  ω, beta)
    forces = integrate_pressure(floatingbody::FloatingBody, pressure) # this is a Dict with keys associated with influenced dofs
    for inf_dof in keys(forces) # influenced dofs
            force = forces[inf_dof]
            F_FK[(ω, inf_dof)] = force
    end
    return F_FK  
end

# Old version
function FroudeKrylovForce(mesh::Mesh, ω,dof,beta=0)
    """Compute the Froude-Krylov force."""
    pressure =  airy_waves_pressure(mesh.centers,  ω, beta)
    return  integrate_pressure(mesh::Mesh, pressure, dof) 
end


function DiffractionForce(floatingbody::FloatingBody,ω,beta=0)
    F_D = Dict{Tuple{Any, String}, Any}()
    mesh = floatingbody.mesh
    green_functions = (
        Rankine(),
        RankineReflected(),
        GFWu(),
    )
    k = ω^2 / SETTINGS.g
    S, D = assemble_matrices(green_functions, mesh, k)
    bc = AiryBC(mesh, ω, beta)
    potential = solve(D, S, bc)
    forces = diffraction_force(floatingbody, potential,ω)
    for inf_dof in keys(forces) # influenced dofs
            force = forces[inf_dof]
            F_D[(ω, inf_dof)] = force
    end
    return F_D
end



# Old version
function DiffractionForce(mesh::Mesh,ω,dof,beta=0)
    green_functions = (
        Rankine(),
        RankineReflected(),
        GFWu(),
    )
    k = ω^2 / SETTINGS.g
    S, D = assemble_matrices(green_functions, mesh, k)
    bc = AiryBC(mesh, ω, beta)
    potential = solve(D, S, bc)
    forces = diffraction_force(potential,mesh,ω,dof)
    return forces
end

# F1 = -sum(dPressure .* sphere_1_heave_normal .* mesh.areas)
# normal_dof_amp = -sum(transpose(dof) .* mesh.normals, dims=2)
#forces = sum(pressure .* normal_dof_amp .* mesh.areas)



function diffraction_force(floatingbody::FloatingBody,potential,omega)
    pressure = 1im*SETTINGS.rho* potential * omega 
    forces = integrate_pressure(floatingbody::FloatingBody, pressure) # this is a Dict with keys associated with influenced dofs
    return forces  
end


# Old version
function diffraction_force(potential::Matrix{ComplexF64},mesh::Mesh,omega,dof)
    pressure = 1im*SETTINGS.rho* potential * omega 
    forces = integrate_pressure(mesh,pressure,dof) 
    return forces  
end
