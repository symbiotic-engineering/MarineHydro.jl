using MarineHydro
using PyCall
using LinearAlgebra

cpt = pyimport("capytaine")
cptmesh = cpt.mesh_sphere(radius=1.0, center=(0, 0, 0), resolution=(10, 10)).immersed_part()

mesh = Mesh(cptmesh + cptmesh.translated_x(5.0))
total_nfaces = mesh.nfaces

face_is_in_sphere_1(j) = j <= Int(total_nfaces/2)
face_is_in_sphere_2(j) = j > Int(total_nfaces/2)

# Just a check that the functions above are correct for this problem
for j in 1:total_nfaces
    if face_is_in_sphere_1(j)
        @assert mesh.centers[j, 1] < 1.0
    elseif face_is_in_sphere_2(j)
        @assert mesh.centers[j, 1] > 4.0
    else
        @error "That should not happen"
    end
end

omega = 2.03
k = omega^2 / 9.8
S, D = assemble_matrices((Rankine(), RankineReflected(), GFWu()), mesh, k)

sphere_1_heave_normal = [face_is_in_sphere_1(j) ? element(mesh, j).normal' * [0, 0, 1] : 0.0 for j in 1:total_nfaces]
sphere_2_heave_normal = [face_is_in_sphere_2(j) ? element(mesh, j).normal' * [0, 0, 1] : 0.0 for j in 1:total_nfaces]

Am = zeros((2, 2))
Bm = zeros((2, 2))

# radiation of first sphere
potential = MarineHydro.solve(D, S, -1im * omega * sphere_1_heave_normal)
pressure = 1im * 1000 * omega * potential 
force_on_sphere_1 = -sum(pressure .* sphere_1_heave_normal .* mesh.areas)
A11 = real(force_on_sphere_1)/omega^2
B11 = imag(force_on_sphere_1)/omega
force_on_sphere_2 = -sum(pressure .* sphere_2_heave_normal .* mesh.areas)
A12 = real(force_on_sphere_2)/omega^2
B12 = imag(force_on_sphere_2)/omega
Am = [A11 A12; A12 A11]
Bm = [B11 B12; B12 B11]
# # radiation of second sphere (very symmetric here, but that is just for demonstration)
# potential = MarineHydro.solve(D, S, -1im * omega * sphere_2_heave_normal)
# pressure = 1im * 1000 * omega * potential 
# force_on_sphere_1 = -sum(pressure .* sphere_1_heave_normal .* mesh.areas)
# Am[2, 1] = real(force_on_sphere_1)/omega^2
# Bm[2, 1] = imag(force_on_sphere_1)/omega
# force_on_sphere_2 = -sum(pressure .* sphere_2_heave_normal .* mesh.areas)
# Am[2, 2] = real(force_on_sphere_2)/omega^2
# Bm[2, 2] = imag(force_on_sphere_1)/omega

@show Bm

## Just for reference
cptbody = cpt.FloatingBody(cptmesh, name="sphere")
cptbody.add_translation_dof(name="Heave")
cpt_twobodies = cptbody + cptbody.translated_x(5.0, name="other_sphere")
xr = pyimport("xarray")
test_matrix = xr.Dataset(coords=Dict("omega" => [omega], "radiating_dof" => collect(keys(cpt_twobodies.dofs))))
ds = cpt.BEMSolver().fill_dataset(test_matrix, cpt_twobodies, method="direct")
@show ds.added_mass.values
@show ds.radiation_damping.values
