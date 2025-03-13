
struct RankineReflected <: GreensFunction end

function greens(::RankineReflected, element_1, element_2, wavenumber=nothing)
    return greens(Rankine(), free_surface_symmetry(element_1), element_2)
    # = greens(Rankine(), element_1, free_surface_symmetry(element_2))
end

function gradient_greens(::RankineReflected, element_1, element_2, wavenumber=nothing; with_respect_to_first_variable=false)
    if with_respect_to_first_variable
        return gradient_greens(Rankine(), element_1, free_surface_symmetry(element_2); with_respect_to_first_variable)
        # = vertical_reflection(gradient_greens(Rankine(), free_surface_symmetry(element_1), element_2; with_respect_to_first_variable))
    else
        return gradient_greens(Rankine(), free_surface_symmetry(element_1), element_2; with_respect_to_first_variable)
        # = vertical_reflection(gradient_greens(Rankine(), element_1, free_surface_symmetry(element_2); with_respect_to_first_variable))
    end
end

function integral(::RankineReflected, element_1, element_2, wavenumber=nothing)
    return integral(Rankine(), free_surface_symmetry(element_1), element_2)
    # = integral(Rankine(), element_1, free_surface_symmetry(element_2))
end

function integral_gradient(::RankineReflected, element_1, element_2, wavenumber=nothing; with_respect_to_first_variable=false)
    if with_respect_to_first_variable
        return integral_gradient(Rankine(), element_1, free_surface_symmetry(element_2); with_respect_to_first_variable)
        # = vertical_reflection(integral_gradient(Rankine(), free_surface_symmetry(element_1), element_2; with_respect_to_first_variable))
    else
        return integral_gradient(Rankine(), free_surface_symmetry(element_1), element_2; with_respect_to_first_variable)
        # = vertical_reflection(integral_gradient(Rankine(), element_1, free_surface_symmetry(element_2); with_respect_to_first_variable))
    end
end
