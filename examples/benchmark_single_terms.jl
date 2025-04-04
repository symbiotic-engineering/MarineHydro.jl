# The goal of this script is to evaluate the time to compute a single term of
# the interaction matrices, for each term of the Green function and for each
# kind of integral.

using MarineHydro
using BenchmarkTools
using StaticArrays

wavenumber = 1.0
green_functions = ["Rankine" => Rankine(), "Wu" => GFWu(), "ExactDelhommeau" => ExactGuevelDelhommeau()]
integrals = ["S   " => MarineHydro.integral, "D   " => MarineHydro.integral_gradient, "both" => MarineHydro.both_integral_and_integral_gradient]

function benchmark(element_1, element_2)
    for (gf_name, gf) in green_functions
        println(" " * gf_name)
        for (term_name, term) in integrals
            print("  " * term_name)
            @btime ($term)($(gf), $element_1, $element_2, $wavenumber)
        end
    end
end

println("WITH STATIC VECTORS")
element_1 = MarineHydro.StaticElement(
    SVector(0.0, 0.0, -1.0),
    @SMatrix([-0.5 -0.5 0.0; 0.5 -0.5 0.0; 0.5 0.5 0.0; -0.5 0.5 0.0]) .+ SVector(0.0, 0.0, -1.0)',
    SVector(0.0, 0.0, 1.0),
    1.0,
    sqrt(2)/2,
)
element_2 = MarineHydro.StaticElement(
    SVector(1.0, 1.0, -2.0),
    @SMatrix([-0.5 -0.5 0.0; 0.5 -0.5 0.0; 0.5 0.5 0.0; -0.5 0.5 0.0]) .+ SVector(1.0, 1.0, -2.0)',
    SVector(0.0, 0.0, 1.0),
    1.0,
    sqrt(2)/2,
    )

benchmark(element_1, element_2)

println()
println("WITH NAMED TUPLES OF VECTORS")
element_1 = (center=[0.0, 0.0, -1.0],)
element_2 = (
    center=[1.0, 1.0, -2.0],
    vertices= [-0.5 -0.5 0.0; 0.5 -0.5 0.0; 0.5 0.5 0.0; -0.5 0.5 0.0] .+ [1.0, 1.0, -2.0]',
    normal=[0.0, 0.0, 1.0],
    radius=sqrt(2)/2,
    area=1.0,
)

benchmark(element_1, element_2)

nothing
