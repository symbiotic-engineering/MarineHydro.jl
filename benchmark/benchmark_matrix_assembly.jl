using BenchmarkTools
using MarineHydro

greens_functions = [
    "Rankine" => Rankine(),
    "Full" => (Rankine(), RankineReflected(), GFWu())
    ]

meshes = [
    "Default" => MarineHydro.Mesh(MarineHydro.example_mesh_from_capytaine()),
    "StaticArrays" => MarineHydro.StaticArraysMesh(MarineHydro.example_mesh_from_capytaine()),
]

methods = [
    "Default" => MarineHydro.assemble_matrices,
    "Broadcasting" => MarineHydro.assemble_matrices_,
]

wavenumber = 1.0

suite = BenchmarkGroup()
for (gf_name, gf) in greens_functions
    for (mesh_name, mesh) in meshes
        for (method_name, method) in methods
            suite[gf_name][mesh_name][method_name] = @benchmarkable ($method)($(gf), $mesh, $wavenumber)
        end
    end
end

tune!(suite)
res = run(suite)
