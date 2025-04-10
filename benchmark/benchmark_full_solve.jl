using BenchmarkTools
using MarineHydro
using CUDA
using DataFrames
using Plots

using LinearAlgebra.BLAS
BLAS.set_num_threads(1)

greens_functions = (Rankine(), RankineReflected(), GFWu())

assembly_methods = [
    "GPU" => mesh -> MarineHydro.assemble_matrices_broadcasting(greens_functions, mesh, 1.0; arrtype=cu),
    "CPU" => mesh -> MarineHydro.assemble_matrices_broadcasting(greens_functions, mesh, 1.0),
]

resolutions = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90]

suite = BenchmarkGroup()
for (assembly_method_name, assembly_method) in assembly_methods
    for res in resolutions
        mesh = MarineHydro.StaticArraysMesh(MarineHydro.example_mesh_from_capytaine(res))
        bc = ones(mesh.nfaces)
        if assembly_method_name == "GPU"
            bc = cu(bc)
        end
        suite[assembly_method_name][mesh.nfaces] = @benchmarkable begin
            @CUDA.sync begin
                S, D = ($assembly_method)($mesh)
                # solve(D, S, $bc)
                D \ $bc
            end
        end
    end
end

# tune!(suite)
results = run(suite, verbose=true)


df = DataFrame(
               (method=met, nfaces=nfaces, time=mean(val).time)
                 for (met, met_results) in results for (nfaces, val) in met_results)

df = sort(df, [:method, :nfaces])

fig = plot(xscale=:log10, xlabel="Number of faces", yscale=:log10, ylabel="Full resolution time (s)")
for (group_name, group) in pairs(groupby(df, :method))
    plot!(group[!, :nfaces], group[!, :time]/1e9, label=group_name[:method])
end
fig
