# The goal of this script is to evaluate the time to compute a single term of
# the interaction matrices, for each term of the Green function and for each
# kind of integral.

using BenchmarkTools
using MarineHydro
using StaticArrays

wavenumber = 1.0

green_functions = ["Rankine" => Rankine(), "Wu" => GFWu(), "ExactDelhommeau" => ExactGuevelDelhommeau()]
integrals = ["S" => MarineHydro.integral, "D" => MarineHydro.integral_gradient, "both" => MarineHydro.both_integral_and_integral_gradient]

static_element_1 = MarineHydro.StaticElement(
    SVector(0.0, 0.0, -1.0),
    @SMatrix([-0.5 -0.5 0.0; 0.5 -0.5 0.0; 0.5 0.5 0.0; -0.5 0.5 0.0]) .+ SVector(0.0, 0.0, -1.0)',
    SVector(0.0, 0.0, 1.0),
    1.0,
    sqrt(2)/2,
)
static_element_2 = MarineHydro.StaticElement(
    SVector(1.0, 1.0, -2.0),
    @SMatrix([-0.5 -0.5 0.0; 0.5 -0.5 0.0; 0.5 0.5 0.0; -0.5 0.5 0.0]) .+ SVector(1.0, 1.0, -2.0)',
    SVector(0.0, 0.0, 1.0),
    1.0,
    sqrt(2)/2,
    )

element_1 = (center=[0.0, 0.0, -1.0],)
element_2 = (
    center=[1.0, 1.0, -2.0],
    vertices= [-0.5 -0.5 0.0; 0.5 -0.5 0.0; 0.5 0.5 0.0; -0.5 0.5 0.0] .+ [1.0, 1.0, -2.0]',
    normal=[0.0, 0.0, 1.0],
    radius=sqrt(2)/2,
    area=1.0,
)

suite = BenchmarkGroup()
for (gf_name, gf) in green_functions
    for (term_name, term) in integrals
        suite["StaticElement"][gf_name][term_name] = @benchmarkable ($term)($(gf), $static_element_1, $static_element_2, $wavenumber)
        suite["NamedTuple"][gf_name][term_name] = @benchmarkable ($term)($(gf), $element_1, $element_2, $wavenumber)
    end
end

tune!(suite)
res = run(suite)
BenchmarkTools.save("latest_results.json", res)
