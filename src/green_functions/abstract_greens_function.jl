
abstract type GreensFunction end

function greens end
function gradient_greens end
function integral end
function integral_gradient end


function greens(gf_tuple::NTuple{N, GreensFunction} where N, element_1, element_2, wavenumber)
    return sum(Complex(greens(gf, element_1, element_2, wavenumber)) for gf in gf_tuple)
end

function gradient_greens(gf_tuple::NTuple{N, GreensFunction} where N, element_1, element_2, wavenumber; with_respect_to_first_variable=false)
    return sum(Vector{ComplexF64}(gradient_greens(gf, element_1, element_2, wavenumber; with_respect_to_first_variable)) for gf in gf_tuple)
end

function integral(gf_tuple::NTuple{N, GreensFunction} where N, element_1, element_2, wavenumber)
    return sum(Complex(integral(gf, element_1, element_2, wavenumber)) for gf in gf_tuple)
end

function integral_gradient(gf_tuple::NTuple{N, GreensFunction} where N, element_1, element_2, wavenumber; with_respect_to_first_variable=false)
    return sum(Vector{ComplexF64}(integral_gradient(gf, element_1, element_2, wavenumber; with_respect_to_first_variable)) for gf in gf_tuple)
end