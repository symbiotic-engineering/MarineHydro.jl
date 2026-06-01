

# Solve single problem (one frequency and one radiating dof or wave direction)
function solve_problem(problem::LinearPotentialFlowProblem; direct::Bool=true, gf::String="Wu")

    # Apply boundary conditions 
    bc_on_hull = vec(compute_bc(problem)) 

    # Specify omega and wavenumber based on forward speed 
    if problem.forward_speed==0
        omega = problem.omega
        wavenumber = problem.wavenumber
    else
        omega = problem.encountered_omega
        wavenumber = problem.encountered_wavenumber # use encountered_wavenumber in gfs
    end

    # Select a Greens function
    if gf=="Wu"  
        selected_GF = GFWu()
    elseif gf=="ExactGuevelDelhommeau"
        selected_GF = ExactGuevelDelhommeau()
    end

    # combine mesh and lid_mesh and define bc accordingly
    if isnothing(problem.floatingbody.lid_mesh)
        mesh_including_lid = problem.floatingbody.mesh
        hull_mask = trues(mesh_including_lid.nfaces)
    else
        mesh_including_lid = problem.floatingbody.mesh + problem.floatingbody.lid_mesh 
        hull_mask = eachrow(mesh_including_lid.faces) .∈ Ref(eachrow(problem.floatingbody.mesh.faces))
    end
    bc = zeros(eltype(bc_on_hull), mesh_including_lid.nfaces)
    bc[hull_mask] = bc_on_hull

    # Create influneced matrices (include lid)
    S, D = assemble_matrices([Rankine(), RankineReflected(), selected_GF], mesh_including_lid, wavenumber; direct=direct)

    # Solve linear system (include lid)
    potential, sources = solve(D, S, bc; direct=direct)

    # Compute pressure (include lid)
    pressure = 1im * SETTINGS.rho * omega * potential # uses encountered_omega

    # Slice pressure (excludes lid)
    pressure_on_hull = pressure[hull_mask]

    if problem.forward_speed!=0
        # change normals to all be unit vector in x direction
        S, D = assemble_matrices([Rankine(), RankineReflected(), selected_GF], mesh_including_lid, wavenumber; direct=direct, all_normals=[1,0,0])
        if direct
            error("MarineHydro.jl has yet to be developed for nonzero forward speeds 
            with the direct method. Try changing direct to false.")
            # partial_phi_partial_x = S \ (D * potential)
            # sources = S \ potential          
            # partial_phi_partial_x = K * sources
        else
            K = D
            partial_phi_partial_x = K * sources
        end
        additional_pressure = SETTINGS.rho * problem.forward_speed * partial_phi_partial_x
        pressure_on_hull .+= additional_pressure[hull_mask]
    end

    forces = integrate_pressure(problem.floatingbody, problem.influenced_dofs, pressure_on_hull) # NamedTuple of complex forces, where each element corresponds to an influenced dof 
    
    result = make_result(problem, forces)
    return result 
end

# Solve multiple problems (multiple frequencies, radiating dofs, and/or wave directions)
# Equivalent to Capytaine's solve_all() function. Eventually add parallelization  settings here.
function solve_all_problems(problems::Vector{LinearPotentialFlowProblem}; direct::Bool=true, gf::String="Wu")
    
    results = [solve_problem(problem; direct=direct, gf=gf) for problem in problems]
    
    return results
end