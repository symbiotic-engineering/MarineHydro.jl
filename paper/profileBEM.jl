using MarineHydro
using BenchmarkTools
include("BemProgram.jl")
include("MeshGradients_singlebody.jl")
## For other branches of the code
#
# element_1 = (center=[0.0, 0.0, -1.0],)
# element_2 = (
#     center=[1.0, 1.0, -2.0],
#     vertices=[-0.5 -0.5 0.0; 0.5 -0.5 0.0; 0.5 0.5 0.0; -0.5 0.5 0.0] .+ [1.0, 1.0, -2.0]',
#     normal=[0.0, 0.0, 1.0],
#     radius=sqrt(2)/2,
#     area=1.0,
# )

element_1 = (center=[0.0, 0.0, -1.0],)
element_2 = (
    center=[1.0, 1.0, -2.0],
    vertices= [-0.5 -0.5 0.0; 0.5 -0.5 0.0; 0.5 0.5 0.0; -0.5 0.5 0.0] .+ [1.0, 1.0, -2.0]',
    normal=[0.0, 0.0, 1.0],
    radius=sqrt(2)/2,
    area=1.0,
)

function integral_and_gradient(GF, element_1, element_2, wavenumber=nothing)
    return (
        MarineHydro.integral(GF, element_1, element_2, wavenumber),
        MarineHydro.integral_gradient(GF, element_1, element_2, wavenumber)
    )
end
omega = 1.03

println("Rankine")
@btime MarineHydro.integral($(Rankine()), $element_1, $element_2)
@btime MarineHydro.integral_gradient($(Rankine()), $element_1, $element_2)
@btime integral_and_gradient($(Rankine()), $element_1, $element_2)

wavenumber = 1.0
println("GFWu")
@btime MarineHydro.integral($(GFWu()), $element_1, $element_2, $wavenumber)
@btime MarineHydro.integral_gradient($(GFWu()), $element_1, $element_2, $wavenumber)
@btime integral_and_gradient($(GFWu()), $element_1, $element_2, $wavenumber)


println("Gradients of coefficients time")
dof = [0,0,1]
A(radius) = added_mass_program(radius,omega,dof)  
B(radius) = damping_program(radius,omega,dof)



@show A(1.0)
@btime Zygote.gradient(A, 1.0)

nothing
