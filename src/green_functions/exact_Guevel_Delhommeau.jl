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
    x  = center(element_1)
    xi = center(element_2)
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
    x = center(element_1)
    xi = center(element_2)
    r  = wavenumber * hypot(x[1] - xi[1], x[2] - xi[2])
    z  = wavenumber * (x[3] + xi[3])
    dGF_dr = _d_dimless_wave_term_dr(r, z)
    dGF_dz = _d_dimless_wave_term_dz(r, z)
    if with_respect_to_first_variable
        if abs(r) > 1e-6
            dr_dx1 = wavenumber^2 / r * (x[1] - xi[1])
            dr_dx2 = wavenumber^2 / r * (x[2] - xi[2])
        else
            dr_dx1 = zero(x[1])
            dr_dx2 = zero(x[2])
        end
        dz_dx = wavenumber
        return wavenumber * (zero(x) .+ (dr_dx1 * dGF_dr, dr_dx2 * dGF_dr, dz_dx * dGF_dz))
        # The zero(x) is a workaround to convert the following tuple to the same type as `x` (either Vector or SVector).
    else
        if abs(r) > 1e-6
            dr_dxi1 = wavenumber^2 / r * (xi[1] - x[1])
            dr_dxi2 = wavenumber^2 / r * (xi[2] - x[2])
        else
            dr_dxi1 = zero(x[1])
            dr_dxi2 = zero(x[2])
        end
        dz_dxi = wavenumber
        return wavenumber * (zero(x) .+ (dr_dxi1 * dGF_dr, dr_dxi2 * dGF_dr, dz_dxi * dGF_dz))
        # The zero(x) is a workaround to convert the following tuple to the same type as `x` (either Vector or SVector).
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
