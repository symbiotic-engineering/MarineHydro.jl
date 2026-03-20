

# Solve single problem (one frequency and one radiating dof or wave direction)
function solve_problem(problem::LinearPotentialFlowProblem)
    k = problem.omega^2 / SETTINGS.g 
    S, D = assemble_matrix_wu(problem.floatingbody.mesh, k)
    bc = compute_bc(problem)
    potential = solve(D, S, bc)
    pressure = 1im * SETTINGS.rho * problem.omega * potential
    forces = integrate_pressure(problem.floatingbody, pressure) # NamedTuple of complex forces, where each element corresponds to a dof 
    result = make_result(problem, forces)
    return result 
end

# Solve multiple problems (multiple frequencies, radiating dofs, and/or wave directions)
# Equivalent to Capytaine's solve_all() function. Eventually add parallelization  settings here.
function solve_all_problems(problems::Vector{LinearPotentialFlowProblem})
    
    results = Vector{LinearPotentialFlowResult}(undef, length(problems))
    
    for i in 1:length(problems)
        results[i] = solve_problem(problems[i])
    end
    
    return results
end