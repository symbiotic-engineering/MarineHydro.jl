using MarineHydro
using PyCall
using CUDA
using Test

cpt = pyimport("capytaine")
radius = 1.0
resolution = (5, 5)
cptmesh = cpt.mesh_sphere(radius=radius, center=(0, 0, 0), resolution=resolution, name="floating sphere").immersed_part()
mesh = Mesh(cptmesh)
smesh = MarineHydro.StaticArraysMesh(cptmesh)

greens_functions = (Rankine(), RankineReflected(), GFWu())

S, D = assemble_matrices(greens_functions, mesh, 1.0)
S_, D_ = MarineHydro.assemble_matrices_(greens_functions, mesh, 1.0)
S__, D__ = MarineHydro.assemble_matrices_(greens_functions, smesh, 1.0)
S_gpu, D_gpu = MarineHydro.assemble_matrices_(greens_functions, smesh, 1.0; arrtype=CuArray)


@test S ≈ S_ ≈ S__ ≈ Array(S_gpu)
@test D ≈ D_ ≈ D__ ≈ Array(D_gpu)
