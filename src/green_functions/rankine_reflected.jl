
struct RankineReflected <: GreensFunction end

function free_surface_symmetry(e)
    if :vertices in keys(e)
        return (
            center=[e.center[1], e.center[2], -e.center[3]],
            vertices=[
                 e.vertices[4, 1] e.vertices[4, 2] -e.vertices[4, 3]
                 e.vertices[3, 1] e.vertices[3, 2] -e.vertices[3, 3];
                 e.vertices[2, 1] e.vertices[2, 2] -e.vertices[2, 3];
                 e.vertices[1, 1] e.vertices[1, 2] -e.vertices[1, 3];
                 ],  # Inverting order such that order is still consistent with normal vector
            normal=[e.normal[1], e.normal[2], -e.normal[3]],
            radius=e.radius,
            area=e.area,
        )
    else
        (
            center=[e.center[1], e.center[2], -e.center[3]],
        )
    end
end


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
