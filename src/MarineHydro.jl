
module MarineHydro

using StaticArrays
using LinearAlgebra
using LinearAlgebra: cross, dot, norm
using ImplicitAD: implicit_linear

const τ̅ = 2π

include("green_functions/abstract_greens_function.jl")
export greens, gradient_greens, integral, integral_gradient
include("green_functions/rankine.jl")
export Rankine
include("green_functions/rankine_reflected.jl")
export RankineReflected
include("green_functions/rankine_reflected_negative.jl")
export RankineReflectedNegative
include("green_functions/wu.jl")
export GFWu
include("green_functions/exact_Guevel_Delhommeau.jl")
export ExactGuevelDelhommeau

include("meshes.jl")
export Mesh, element

include("matrix_assembly.jl")
export assemble_matrices, assemble_matrix_wu, solve

include("waves.jl")
export calculate_radiation_forces, integrate_pressure, DiffractionForce
export FroudeKrylovForce, AiryBC, airy_waves_pressure, airy_waves_velocity,airy_waves_potential
export radiation_bc, solve_hydro_coefficients 

end
