abstract type GreensFunction end

"""
    greens(::GreensFunction, element_1, element_2, wavenumber=nothing)

Calculates the Green's function between two points.

# Arguments
- `element_1`: first point, as a `NamedTuple` with the field `center`.
- `element_2`: second point, as a `NamedTuple` with the field `center`.
- `wavenumber`: Defines the wavenumber for waves ~~~~~. Optional for some Greens functions such as `Rankine()`.
"""
function greens end

"""
    gradient_greens(::GreensFunction, element_1, element_2, wavenumber; with_respect_to_first_variable=false)

Calculates the gradient of the Rankine Green's function between two points.

# Arguments
- `element_1`: first point, as a `NamedTuple` with the field `center`.
- `element_2`: second point, as a `NamedTuple` with the field `center`.
- `wavenumber`: Defines the wavenumber for waves ~~~~~. Optional for some Greens functions such as `Rankine()`.
- `with_respect_to_first_variable`: A boolean flag (default is `false`). If `true`, computes the gradient with respect to `element_1`; otherwise, computes the gradient with respect to `element_2`.
"""
function gradient_greens end

"""
    integral(::GreensFunction, element_1, element_2, wavenumber)

Calculates the integral of Green's function over a panel.

# Arguments
- `element_1`: first point, as a `NamedTuple` with the field `center`.
- `element_2`: integration panel, as a `NamedTuple` with the fields `center`, `normal`, `vertices`, `radius`, `normal` and `area`.
- `wavenumber`: Defines the wavenumber for waves ~~~~~. Optional for some Greens functions such as `Rankine()`.
"""
function integral end

"""
    integral_gradient(::GreensFunction, element_1, element_2, wavenumber; with_respect_to_first_variable=false)

Calculates the integral of the gradient of Green's function over a panel.

# Arguments
- `element_1`: first point, as a `NamedTuple` with the field `center`.
- `element_2`: integration panel, as a `NamedTuple` with the fields `center`, `normal`, `vertices`, `radius`, `normal` and `area`.
- `wavenumber`: Defines the wavenumber for waves ~~~~~. Optional for some Greens functions such as `Rankine()`.
- `with_respect_to_first_variable`: A boolean flag (default is `false`). If `true`, computes the gradient with respect to `element_1`; otherwise, computes the gradient with respect to `element_2`.
"""
function integral_gradient end

