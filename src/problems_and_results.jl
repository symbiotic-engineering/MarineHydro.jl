using DimensionalData

########################## Problems #########################################

# Abstract type named LinearPotentialFlowProblem. The diffraction and radiation problems will be subtypes of this type. 
abstract type LinearPotentialFlowProblem end

# Define DiffractionProblem struct as a subtype of LinearPotentialFlowProblem
struct DiffractionProblem <: LinearPotentialFlowProblem
    floatingbody::FloatingBody
    omega::Real
    beta::Real
    function DiffractionProblem(floatingbody::FloatingBody,
        omega::Real,
        beta::Real)
        return new(floatingbody, omega, beta)
    end
end

# Define RadiationProblem struct as a subtype of LinearPotentialFlowProblem
struct RadiationProblem <: LinearPotentialFlowProblem
    floatingbody::FloatingBody
    omega::Real
    radiating_dof::Symbol
    function RadiationProblem(floatingbody::FloatingBody,
        omega::Real,
        radiating_dof::Symbol)
        @assert (radiating_dof in keys(floatingbody.dofs)) "the dof Symbol must be a key of floatingbody.dof"
        return new(floatingbody, omega, radiating_dof)
    end
end

########################## Results #########################################

abstract type LinearPotentialFlowResult end

struct DiffractionResult <: LinearPotentialFlowResult
    problem::DiffractionProblem
    forces::NamedTuple
    function DiffractionResult(problem::DiffractionProblem,
        forces::NamedTuple)
        return new(problem, forces)
    end
end
struct RadiationResult <: LinearPotentialFlowResult
    problem::RadiationProblem
    forces::NamedTuple
    function RadiationResult(problem::RadiationProblem,
        forces::NamedTuple)
        return new(problem, forces)
    end
end

# Convert problem and forces for that problem into a results struct
function make_result(problem::RadiationProblem, forces::NamedTuple)
    return RadiationResult(problem,forces)
end
function make_result(problem::DiffractionProblem, forces::NamedTuple)
    return DiffractionResult(problem,forces)
end


# Convert parameters and problem into a Vector of problems
function problems_from_data(parameters::NamedTuple, floatingbody::FloatingBody)

    # There is at least one diffraction problem to solve
    if :wave_directions in keys(parameters)
        diffraction_problems = vec([DiffractionProblem(floatingbody, omega, beta) 
            for beta in parameters[:wave_directions], 
                omega in parameters[:wave_frequencies]])
    else
        diffraction_problems = []
    end

    # There is at least one radiation problem to solve
    if :radiating_dofs in keys(parameters)
        radiation_problems = vec([RadiationProblem(floatingbody, omega, dof) 
            for dof in parameters[:radiating_dofs], 
                omega in parameters[:wave_frequencies]])
    else
        radiation_problems = []

    end

    return vcat(diffraction_problems, radiation_problems)
end


# Convert Vector of results into DimStack
# assemble_data automatically determines what outputs to compute based on what parameters are specified. 
function assemble_data(parameters::NamedTuple, floatingbody::FloatingBody, results::Vector{LinearPotentialFlowResult})

    omegas = parameters[:wave_frequencies]
    betas = parameters[:wave_directions]
    inf_dofs = keys(floatingbody.dofs)
    rad_dofs = keys(floatingbody.dofs)

    added_mass_data = zeros(length(inf_dofs),
        length(rad_dofs),
        length(omegas))
    radiation_damping_data = zeros(length(inf_dofs),
        length(rad_dofs),
        length(omegas))
    excitation_force_data = zeros(ComplexF64,
        length(inf_dofs),
        length(omegas),
        length(betas))
    diffraction_force_data = zeros(ComplexF64,
        length(inf_dofs),
        length(omegas),
        length(betas))
    Froude_Krylov_force_data = zeros(ComplexF64,
        length(inf_dofs),
        length(omegas),
        length(betas))

    for result in results
        if result isa DiffractionResult
            omega = result.problem.omega
            beta = result.problem.beta
            omega_index = findfirst(==(omega),omegas)
            beta_index = findfirst(==(beta),betas)
            diffraction_force_data[:,omega_index,beta_index] = collect(result.forces)
            Froude_Krylov_force_data[:,omega_index,beta_index] = collect(FroudeKrylovForce(floatingbody,omega,beta))
            excitation_force_data[:,omega_index,beta_index] = diffraction_force_data[:,omega_index,beta_index] + Froude_Krylov_force_data[:,omega_index,beta_index]

        elseif result isa RadiationResult
            rad_dof = result.problem.radiating_dof
            omega = result.problem.omega
            rad_dof_index = findfirst(==(rad_dof),rad_dofs)
            omega_index = findfirst(==(omega),omegas)
            added_mass_data[:,rad_dof_index,omega_index] = real(collect(result.forces))/omega^2
            radiation_damping_data[:,rad_dof_index,omega_index] = imag(collect(result.forces))/omega
        end
    end  

    radiation_dims = (Dim{:influenced_dofs}(collect(inf_dofs)), 
        Dim{:radiating_dofs}(collect(rad_dofs)),
        Dim{:wave_frequencies}(omegas))

    diffraction_dims = (Dim{:influenced_dofs}(collect(inf_dofs)),
        Dim{:wave_frequencies}(collect(omegas)),
        Dim{:wave_directions}(betas))


    added_mass_array = DimArray(added_mass_data, radiation_dims)
    radiation_damping_array = DimArray(radiation_damping_data, radiation_dims)
    excitation_force_array = DimArray(excitation_force_data, diffraction_dims)
    diffraction_force_array = DimArray(diffraction_force_data, diffraction_dims)
    Froude_Krylov_force_array = DimArray(Froude_Krylov_force_data, diffraction_dims)


    data = DimStack((
        added_mass = added_mass_array,
        radiation_damping = radiation_damping_array,
        excitation_force = excitation_force_array,
        diffraction_force = diffraction_force_array,
        Froude_Krylov_force = Froude_Krylov_force_array))
    return data 
end



# Compute DimStack of reuslts based on NamedTuple of parameters (with keys wave_frequencies, wave_directions, and radiating_dofs) and floatingbody struct
function compute_hydrodynamic_coefficients(parameters::NamedTuple, floatingbody::FloatingBody)
    problems = problems_from_data(parameters, floatingbody)
    results = solve_all_problems(problems)
    data = assemble_data(parameters, floatingbody, results)
    return data
end