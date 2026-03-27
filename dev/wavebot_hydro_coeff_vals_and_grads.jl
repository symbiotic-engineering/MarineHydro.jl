using ChainRulesCore
using FiniteDifferences
using PyCall
using Zygote
using MarineHydro
using Plots
using Base.Threads

cpt = pyimport("capytaine")
np = pyimport("numpy")


# define wavebot profile
function profile_points(r1,r2,d1,d2)

    # nonhorizontal portion
    function radial_profile(z,r1,r2,d1,d2)
        if z>=-d1
            r = r1
        elseif z>=-d2 && z<=-d1
            m = (r1-r2)/(d2-d1)
            b = r1 + (m*d1)
            r = m * z + b
        end
        return r        
    end
    z_vals1 = range(-d2, stop=0.0, length=10)
    r_vals1 = radial_profile.(z_vals1,r1,r2,d1,d2)

    # horizontal portion
    num_points_h = 3
    z_vals2 = - ones(num_points_h) * d2
    r_vals2 = range(r2, stop=0.0, length=num_points_h)

    # concatinate
    r_vals = vcat(r_vals2,r_vals1)
    z_vals = vcat(z_vals2,z_vals1)

    return r_vals, z_vals
end


# Plot profile 
# set arbitrary vals
r1 = 4.0
r2 = 2.0
d1 = 2.0
d2 = 5.0
(r_vals, z_vals) = profile_points(r1,r2,d1,d2)
plot(r_vals, z_vals,
     xlabel="Radial Dimension",
     ylabel="Vertical Dimension",  
     title="Wavebot Profile Geometry",
     lw=2)


# Generate Capytaine Mesh
function cptMesh(r1,r2,d1,d2)
    # r1 = 4.0
    # r2 = 2.0
    # d1 = 2.0
    # d2 = 5.0
    (r_vals, z_vals) = profile_points(r1,r2,d1,d2)

    # 3D points
    x_vals = r_vals
    y_vals = zeros(length(r_vals)) # A column of zeros
    z_vals = z_vals
    pts= np.array(hcat(x_vals, y_vals, z_vals))

    # make cpt mesh
    cptmesh = cpt.AxialSymmetricMesh.from_profile(pts,nphi=5)
    return cptmesh.keep_immersed_part(inplace=true)
end

############################## Zygote info ####################################
# Functions for generating input to MarineHydro Mesh struct 
function get_mesh_sizes(r1,r2,d1,d2) 
    cptmesh = cptMesh(r1,r2,d1,d2)
    return cptmesh.nb_faces, cptmesh.nb_vertices
end
function get_vertices(r1,r2,d1,d2)
    cptmesh = cptMesh(r1,r2,d1,d2)
    return reduce(vcat,cptmesh.vertices)
end
function get_centers(r1,r2,d1,d2)
    cptmesh = cptMesh(r1,r2,d1,d2)
    return reduce(vcat,cptmesh.faces_centers)
end
function get_normals(r1,r2,d1,d2)
    cptmesh = cptMesh(r1,r2,d1,d2)
    return reduce(vcat,cptmesh.faces_normals)
end
function get_radii(r1,r2,d1,d2) 
    cptmesh = cptMesh(r1,r2,d1,d2)
    return reduce(vcat,cptmesh.faces_radiuses)
end
function get_areas(r1,r2,d1,d2) 
    cptmesh = cptMesh(r1,r2,d1,d2)
    return reduce(vcat,cptmesh.faces_areas)
end
function get_faces(r1,r2,d1,d2) 
    cptmesh = cptMesh(r1,r2,d1,d2)
    return reduce(vcat,cptmesh.faces)
end

function ∂J_arg_fd(f, args, i; h=1e-5)
    args_plus = collect(args)
    args_minus = collect(args)
    
    args_plus[i] += h
    args_minus[i] -= h
    
    return (f(args_plus...) .- f(args_minus...)) ./ (2*h)
end

function ChainRulesCore.rrule(::typeof(get_vertices), r1, r2, d1, d2)
    y = get_vertices(r1, r2, d1, d2)
    args = (r1, r2, d1, d2)
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function

        # Calculate partials for each input
        dv_dr1 = ∂J_arg_fd(get_vertices, args, 1)
        dv_dr2 = ∂J_arg_fd(get_vertices, args, 2)
        dv_dd1 = ∂J_arg_fd(get_vertices, args, 3)
        dv_dd2 = ∂J_arg_fd(get_vertices, args, 4)

        # Calculate pullbacks
        dr1 = dv_dr1' * dy
        dr2 = dv_dr2' * dy
        dd1 = dv_dd1' * dy
        dd2 = dv_dd2' * dy
        
        return (df, dr1, dr2, dd1, dd2)  
    end    
    return y, f_pullback
end

function ChainRulesCore.rrule(::typeof(get_centers), r1, r2, d1, d2)
    y = get_centers(r1, r2, d1, d2)
    args = (r1, r2, d1, d2)
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function

        # Calculate partials for each input
        dv_dr1 = ∂J_arg_fd(get_centers, args, 1)
        dv_dr2 = ∂J_arg_fd(get_centers, args, 2)
        dv_dd1 = ∂J_arg_fd(get_centers, args, 3)
        dv_dd2 = ∂J_arg_fd(get_centers, args, 4)

        # Calculate pullbacks
        dr1 = dv_dr1' * dy
        dr2 = dv_dr2' * dy
        dd1 = dv_dd1' * dy
        dd2 = dv_dd2' * dy
        
        return (df, dr1, dr2, dd1, dd2)  
    end    
    return y, f_pullback
end

function ChainRulesCore.rrule(::typeof(get_normals), r1, r2, d1, d2)
    y = get_normals(r1, r2, d1, d2)
    args = (r1, r2, d1, d2)
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function

        # Calculate partials for each input
        dv_dr1 = ∂J_arg_fd(get_normals, args, 1)
        dv_dr2 = ∂J_arg_fd(get_normals, args, 2)
        dv_dd1 = ∂J_arg_fd(get_normals, args, 3)
        dv_dd2 = ∂J_arg_fd(get_normals, args, 4)

        # Calculate pullbacks
        dr1 = dv_dr1' * dy
        dr2 = dv_dr2' * dy
        dd1 = dv_dd1' * dy
        dd2 = dv_dd2' * dy
        
        return (df, dr1, dr2, dd1, dd2)  
    end    
    return y, f_pullback
end

function ChainRulesCore.rrule(::typeof(get_radii), r1, r2, d1, d2)
    y = get_radii(r1, r2, d1, d2)
    args = (r1, r2, d1, d2)
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function

        # Calculate partials for each input
        dv_dr1 = ∂J_arg_fd(get_radii, args, 1)
        dv_dr2 = ∂J_arg_fd(get_radii, args, 2)
        dv_dd1 = ∂J_arg_fd(get_radii, args, 3)
        dv_dd2 = ∂J_arg_fd(get_radii, args, 4)

        # Calculate pullbacks
        dr1 = dv_dr1' * dy
        dr2 = dv_dr2' * dy
        dd1 = dv_dd1' * dy
        dd2 = dv_dd2' * dy
        
        return (df, dr1, dr2, dd1, dd2)  
    end    
    return y, f_pullback
end

function ChainRulesCore.rrule(::typeof(get_areas), r1, r2, d1, d2)
    y = get_areas(r1, r2, d1, d2)
    args = (r1, r2, d1, d2)
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function

        # Calculate partials for each input
        dv_dr1 = ∂J_arg_fd(get_areas, args, 1)
        dv_dr2 = ∂J_arg_fd(get_areas, args, 2)
        dv_dd1 = ∂J_arg_fd(get_areas, args, 3)
        dv_dd2 = ∂J_arg_fd(get_areas, args, 4)

        # Calculate pullbacks
        dr1 = dv_dr1' * dy
        dr2 = dv_dr2' * dy
        dd1 = dv_dd1' * dy
        dd2 = dv_dd2' * dy
        
        return (df, dr1, dr2, dd1, dd2)  
    end    
    return y, f_pullback
end

function ChainRulesCore.rrule(::typeof(get_faces), r1, r2, d1, d2)
    y = get_faces(r1, r2, d1, d2)
    return y, _ -> (NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent())
end

function differentiableMesh(r1, r2, d1, d2)
    nf,nv = Zygote.@ignore get_mesh_sizes(r1, r2, d1, d2) #hack..need to know how many first to reshape
    vertices = reshape(get_vertices(r1, r2, d1, d2),(nv,3))
    faces =  reshape(get_faces(r1, r2, d1, d2),(nf,4))
    centers = reshape(get_centers(r1, r2, d1, d2),(nf,3))
    normals = reshape(get_normals(r1, r2, d1, d2),(nf,3))
    areas = get_areas(r1, r2, d1, d2)
    radii = get_radii(r1, r2, d1, d2)
    nvertices = size(vertices,1)
    nfaces = size(centers,1)
    return Mesh(vertices,faces,centers,normals,areas,radii,nvertices,nfaces)
end

############################## Functions for hydro coeff calcs #################
# Hydrostatics
function hydrostatic_results(r1, r2, d1, d2)
    rho_w = 1025 # density of fluid [kg/m^3]
    g = 9.81 # acceleration due to gravity [m/s^2]
    K = rho_w * g * pi * r1^2 # hydrostatic stiffness [kg/s^2]
    M = rho_w * (pi * r1^2 * d1 + (pi/3) * (d2-d1) * (r1^2 + r1*r2 + r2^2)) # mass of body [kg]
    return [M, K]
end
function hyd_sta_prob(r1, r2, d1, d2)
    j_hydrostatic_results_AD_raw = Zygote.withjacobian(hydrostatic_results,r1, r2, d1, d2) # first row is real pert and second row is imag part
    mass_val = j_hydrostatic_results_AD_raw.val[1]
    mass_grad = hcat(j_hydrostatic_results_AD_raw.grad ...)[1,:]
    hydrostatic_stiffness_val = j_hydrostatic_results_AD_raw.val[2]
    hydrostatic_stiffness_grad = hcat(j_hydrostatic_results_AD_raw.grad ...)[2,:]
    return mass_val, mass_grad, hydrostatic_stiffness_val, hydrostatic_stiffness_grad  
end

# Radiation
function radiation_program(r1, r2, d1, d2, omega, dof) 
    mesh = differentiableMesh(r1, r2, d1, d2) #fd
    A, B = calculate_radiation_forces(mesh,dof,omega)
    return [A, B]
end
function rad_prob(r1, r2, d1, d2, omega, dof)
    radiation_results(r1, r2, d1, d2) = radiation_program(r1, r2, d1, d2, omega, dof)
    j_radiation_results_AD_raw = Zygote.withjacobian(radiation_results,r1, r2, d1, d2) # first row is added mass and second row is damping
    added_mass_val = j_radiation_results_AD_raw.val[1]
    added_mass_grad = hcat(j_radiation_results_AD_raw.grad ...)[1,:]
    damping_val = j_radiation_results_AD_raw.val[2]
    damping_grad = hcat(j_radiation_results_AD_raw.grad ...)[2,:]
    return added_mass_val, added_mass_grad, damping_val, damping_grad 
end

# Diffraction
function diffraction_program(r1, r2, d1, d2, omega, dof) 
    mesh = differentiableMesh(r1, r2, d1, d2) #fd
    F_D = DiffractionForce(mesh,omega,dof)
    F_FK = FroudeKrylovForce(mesh,omega,dof)
    F_ex = F_D + F_FK
    return [real(F_ex),imag(F_ex)]
end
function diff_prob(r1, r2, d1, d2, omega, dof)
    diffraction_results(r1, r2, d1, d2)  = diffraction_program(r1, r2, d1, d2, omega, dof)
    j_diffraction_results_AD_raw = Zygote.withjacobian(diffraction_results,r1, r2, d1, d2) # first row is real pert and second row is imag part
    real_excitation_force_val = j_diffraction_results_AD_raw.val[1]
    real_excitation_force_grad = hcat(j_diffraction_results_AD_raw.grad ...)[1,:]
    imag_excitation_force_val = j_diffraction_results_AD_raw.val[2]
    imag_excitation_force_grad = hcat(j_diffraction_results_AD_raw.grad ...)[2,:]
    return real_excitation_force_val, real_excitation_force_grad, imag_excitation_force_val, imag_excitation_force_grad 
end

# Everything
function compute_all_hydrodynamic_coefficinets_vals_and_grads(r1, r2, d1, d2, omegas)

    # Hydrostatic problem (not frequency dependent)
    mass_val, mass_grad, hydrostatic_stiffness_val, hydrostatic_stiffness_grad  = hyd_sta_prob(r1, r2, d1, d2)
    dof = [0.0,0.0,1.0]

    # Radiation and diffraction problems
    n_omega = length(omegas)
    results = Vector{Any}(undef, n_omega)

    for i in 1:n_omega
    ω = omegas[i]    
    am_v, am_g, dmp_v, dmp_g = rad_prob(r1, r2, d1, d2, ω, dof)
    re_v, re_g, im_v, im_g   = diff_prob(r1, r2, d1, d2, ω, dof)

    results[i] = (
        omega = ω,
        added_mass = (val = am_v, grad = am_g),
        damping    = (val = dmp_v, grad = dmp_g),
        force_real = (val = re_v, grad = re_g),
        force_imag = (val = im_v, grad = im_g)
    )
    end
end

################################## Trying out the function ####################
omegas = range(start=0.005,stop=0.3,length=3)
@timev results = compute_all_hydrodynamic_coefficinets_vals_and_grads(r1, r2, d1, d2, omegas)
# @timev results = compute_all_hydrodynamic_coefficinets_vals_and_grads(r1, r2, d1, d2, omegas)
display(results)