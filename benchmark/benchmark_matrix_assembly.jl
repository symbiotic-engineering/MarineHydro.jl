using BenchmarkTools
using MarineHydro

resolution = 20

greens_functions = [
    "Rankine" => Rankine(),
    "Full" => (Rankine(), RankineReflected(), GFWu())
    ]

meshes = [
    "Mesh" => MarineHydro.Mesh(MarineHydro.example_mesh_from_capytaine(resolution)),
    "StaticArraysMesh" => MarineHydro.StaticArraysMesh(MarineHydro.example_mesh_from_capytaine(resolution)),
]

methods = [
    "Comprehension" => MarineHydro.assemble_matrices_comprehension,
    "Broadcasting" => MarineHydro.assemble_matrices_broadcasting,
    "Explicit Both" => MarineHydro.assemble_matrices_explicit_both,
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

using PyCall
cpt = pyimport("capytaine")
cpt_gf = cpt.Delhommeau()
cpt_mesh = MarineHydro.example_mesh_from_capytaine(resolution)
suite["Full"]["Capytaine"] = @benchmarkable begin
    ($cpt_gf.evaluate)($cpt_mesh, $cpt_mesh, $0.0, $Inf, $wavenumber)
end

tune!(suite)
res = run(suite, verbose=true)
