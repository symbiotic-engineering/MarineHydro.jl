# Slow, but accurate expression of the Green function

struct ExactGuevelDelhommeau <: GreensFunction end

using SpecialFunctions: expint
using Integrals

function integrate(f,a,b,p)
    domain = (a, b)
    prob = IntegralProblem(f, domain,p) 
    sol = Integrals.solve(prob, HCubatureJL())
    return sol.u 
end

function _dimless_wave_term(r, z)
    function integrand(θ,p)
        r,z = p[1],p[2]
        zeta = z + im * r * cos(θ)
        exp_zeta = exp(zeta)
        j_zeta = exp_zeta * (expint(zeta) / π + im)
        Gp = real(j_zeta) + im * real(exp_zeta)
        return Gp
    end
    p = [r,z]
    return 4 * integrate(integrand, 0, π / 2,p)
end

function greens(::ExactGuevelDelhommeau, element_1, element_2, wavenumber)
    with_reduced_coordinates(element_1, element_2, wavenumber) do r, z
        _dimless_wave_term(r, z)
    end
end

"""Derivative of dimless_wave_term with respect to r"""
function _d_dimless_wave_term_dr(r, z)
    function integrand(θ,p)
        r,z = p[1],p[2]
        zeta = z + im*r*cos(θ)
        exp_zeta = exp(zeta)
        j_zeta = exp_zeta*(expint(zeta)/π + im)
        dGpdr = real(im*cos(θ)*(j_zeta - 1/(π*zeta))) + im*real(im*cos(θ)*exp_zeta)
        return dGpdr
    end
    p= [r,z]
    return 4*integrate(integrand, 0, π/2,p)
end

function gradient_greens(gf::ExactGuevelDelhommeau, element_1, element_2, wavenumber; with_respect_to_first_variable=false)
    both_greens_and_gradient_greens(gf, element_1, element_2, wavenumber; with_respect_to_first_variable)[2]
end

function both_greens_and_gradient_greens(::ExactGuevelDelhommeau, element_1, element_2, wavenumber; with_respect_to_first_variable=false)
    with_reduced_coordinates_derivative(element_1, element_2, wavenumber; with_respect_to_first_variable) do r, z
        dGF_dr = _d_dimless_wave_term_dr(r, z)
        GF = _dimless_wave_term(r, z)
        dGF_dz = GF + 2/hypot(r, z)
        return GF, dGF_dr, dGF_dz 
    end
end

function integral(g::ExactGuevelDelhommeau, element_1, element_2, wavenumber)
    # One-point approximation of the integral
    return greens(g, element_1, element_2, wavenumber) * area(element_2)
end

function integral_gradient(g::ExactGuevelDelhommeau, element_1, element_2, wavenumber; with_respect_to_first_variable=false)
    # One-point approximation of the integral
    return gradient_greens(g, element_1, element_2, wavenumber; with_respect_to_first_variable) * area(element_2)
end

function both_integral_and_integral_gradient(g::ExactGuevelDelhommeau, element_1, element_2, wavenumber; with_respect_to_first_variable=false)
    # One-point approximation of the integral
    GF, dGF = both_greens_and_gradient_greens(g::ExactGuevelDelhommeau, element_1, element_2, wavenumber; with_respect_to_first_variable)
    return GF * area(element_2), dGF * area(element_2)
end
