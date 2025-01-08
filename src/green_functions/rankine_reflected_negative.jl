struct RankineReflectedNegative <: GreensFunction end

function greens(::RankineReflectedNegative, element_1, element_2, wavenumber)
    return -greens(RankineReflected(), element_1, element_2, wavenumber)
end

function gradient_greens(::RankineReflectedNegative, element_1, element_2, wavenumber; with_respect_to_first_variable=false)
    return -gradient_greens(RankineReflected(), element_1, element_2, wavenumber; with_respect_to_first_variable)
end

function integral(::RankineReflectedNegative, element_1, element_2, wavenumber)
    return -integral(RankineReflected(), element_1, element_2, wavenumber)
end

function integral_gradient(::RankineReflectedNegative, element_1, element_2, wavenumber; with_respect_to_first_variable=false)
    return -integral_gradient(RankineReflected(), element_1, element_2, wavenumber; with_respect_to_first_variable)
end
