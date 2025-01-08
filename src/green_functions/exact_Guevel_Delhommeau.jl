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
    x  = element_1.center
    xi = element_2.center
    r  = wavenumber * hypot(x[1] - xi[1], x[2] - xi[2])
    z  = wavenumber * (x[3] + xi[3])
    return wavenumber * _dimless_wave_term(r, z)
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

"""Derivative of dimless_wave_term with respect to z"""
function _d_dimless_wave_term_dz(r, z)
    return _dimless_wave_term(r, z) + 2/hypot(r, z)
end

function gradient_greens(::ExactGuevelDelhommeau, element_1, element_2, wavenumber; with_respect_to_first_variable=false)
    x = element_1.center
    xi = element_2.center
    r  = wavenumber * hypot(x[1] - xi[1], x[2] - xi[2])
    z  = wavenumber * (x[3] + xi[3])
    if with_respect_to_first_variable
        if abs(r) > 1e-6
            dr_dx = wavenumber^2 / r * [x[1] - xi[1], x[2] - xi[2]]
        else
            dr_dx = [0.0, 0.0]
        end
        dz_dx = wavenumber
        return wavenumber * [dr_dx * _d_dimless_wave_term_dr(r, z); dz_dx * _d_dimless_wave_term_dz(r, z)]
    else
        if abs(r) > 1e-6
            dr_dxi = wavenumber^2 / r * [xi[1] - x[1], xi[2] - x[2]]
        else
            dr_dxi = [0.0, 0.0]
        end
        dz_dxi = wavenumber
        return wavenumber * [dr_dxi * _d_dimless_wave_term_dr(r, z); dz_dxi * _d_dimless_wave_term_dz(r, z)]
    end
end

function integral(g::ExactGuevelDelhommeau, element_1, element_2, wavenumber)
    # One-point approximation of the integral
    return greens(g, element_1, element_2, wavenumber) * element_2.area
end

function integral_gradient(g::ExactGuevelDelhommeau, element_1, element_2, wavenumber; with_respect_to_first_variable=false)
    # One-point approximation of the integral
    return gradient_greens(g, element_1, element_2, wavenumber; with_respect_to_first_variable) * element_2.area
end
