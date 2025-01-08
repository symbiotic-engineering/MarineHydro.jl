

function added_mass_program(radius,omega,dof)  
    mesh = differentiableMesh(radius) #fd
    A = calculate_radiation_forces(mesh,dof,omega)[1]
    return A
end

function damping_program(radius,omega,dof)  
    mesh = differentiableMesh(radius) #fd
    B = calculate_radiation_forces(mesh,dof,omega)[2]
    return B
end

function diffraction_program(ω,radius,dof)  
    mesh = differentiableMesh(radius) #fd
    force = DiffractionForce(mesh,ω,dof)
    return force
end

function omega_added_mass_bem_program(omega,radius=1,dof = [0,0,1])  
    mesh = differentiableMesh(radius) #fd
    A = calculate_radiation_forces(mesh,dof,omega)[1]
    return A
end

function omega_damping_bem_program(omega,radius=1,dof = [0,0,1])  
    mesh = differentiableMesh(radius) #fd
    B = calculate_radiation_forces(mesh,dof,omega)[2]
    return B
end