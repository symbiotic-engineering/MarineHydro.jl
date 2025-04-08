############
#  Meshes  #
############

struct Mesh
    vertices::AbstractMatrix
    faces::AbstractMatrix
    centers::AbstractMatrix
    normals::AbstractMatrix
    areas::AbstractVector
    radii::AbstractVector
    nvertices::Number
    nfaces::Number
    function Mesh(vertices::AbstractMatrix, faces::AbstractMatrix,
                  centers::AbstractMatrix, normals::AbstractMatrix, areas::AbstractVector,
                  radii::AbstractVector,nvertices::Number,
                  nfaces::Number
                  )
        # TODO: centers, areas, and radii(?) could be calculated from list of vertices and faces
        nvertices = size(vertices)[1]
        @assert (size(vertices) == (nvertices, 3)) "each vertex needs 3 coordinates"
        nfaces = size(faces)[1]
        @assert (size(faces) == (nfaces, 4)) "only quadrilateral panels are allowed"
        @assert (size(centers) == (nfaces, 3)) "centers must be [nfaces x 3]"
        @assert (size(centers) == (nfaces, 3)) "centers must be [nfaces x 3]"
        @assert (length(areas) == nfaces) "areas vector must have length nfaces"
        @assert (length(radii) == nfaces) "radii vector must have length nfaces"
        return new(vertices, faces, centers, normals, areas, radii, nvertices, nfaces)
    end
end

struct StaticArraysMesh
    vertices::Vector{SVector{3, Float64}}
    faces::Vector{SVector{4, Int64}}
    centers::Vector{SVector{3, Float64}}
    normals::Vector{SVector{3, Float64}}
    areas::Vector{Float64}
    radii::Vector{Float64}
    nvertices::Int64
    nfaces::Int64
end

###################################
#  Loading meshes from Capytaine  #
###################################


using PyCall

function example_mesh_from_capytaine()
    cpt = pyimport("capytaine")
    radius = 1.0
    resolution = (4, 4)
    return cpt.mesh_sphere(radius=radius, center=(0, 0, 0), resolution=resolution, name="floating sphere").immersed_part()
end

function Mesh(mesh::PyObject)
    mesh = Mesh(
        mesh.vertices, mesh.faces, 
        mesh.faces_centers, mesh.faces_normals, mesh.faces_areas, mesh.faces_radiuses,mesh.nb_vertices,mesh.nb_faces
    )
    return mesh
end

function StaticArraysMesh(mesh::PyObject)
    mesh = StaticArraysMesh(
        [SVector(v[1], v[2], v[3]) for v in eachrow(mesh.vertices)],
        [SVector(f[1], f[2], f[3], f[4]) .+ 1 for f in eachrow(mesh.faces)],
        [SVector(c[1], c[2], c[3]) for c in eachrow(mesh.faces_centers)],
        [SVector(n[1], n[2], n[3]) for n in eachrow(mesh.faces_normals)],
        mesh.faces_areas,
        mesh.faces_radiuses,
        mesh.nb_vertices,
        mesh.nb_faces
    )
    return mesh
end


##############
#  Elements  #
##############

# Mesh elements can be provided individually as NamedTuples, for instance to test the Green's function
center(nt::NamedTuple) = nt[:center]
vertices(nt::NamedTuple) = nt[:vertices]
radius(nt::NamedTuple) = nt[:radius]
normal(nt::NamedTuple) = nt[:normal]
area(nt::NamedTuple) = nt[:area]
faces(nt::NamedTuple) = nt[:faces]


# When using a Mesh, the element is stored as a reference to the Mesh and a index
struct LazyElement
    mesh::Mesh
    index::Integer
end

function element(mesh::Mesh, J::Int)
    return LazyElement(mesh, J)
end

center(e::LazyElement) = e.mesh.centers[e.index, :] #, e.mesh.centers[e.index, 1], e.mesh.centers[e.index, 2]
normal(e::LazyElement) = e.mesh.normals[e.index, :] #, e.mesh.normals[e.index, 1], e.mesh.normals[e.index, 2]
area(e::LazyElement) = e.mesh.areas[e.index]
radius(e::LazyElement) = e.mesh.radii[e.index]
vertices(e::LazyElement) = e.mesh.vertices[e.mesh.faces[e.index, :] .+ 1, :]


# When using a StaticArraysMesh, the element is stored as a set of static arrays (although LazyElement might also work)
struct StaticElement
    center::SVector{3, Float64}
    vertices::SMatrix{4, 3, Float64, 12}
    normal::SVector{3, Float64}
    area::Float64
    radius::Float64
end

function element(mesh::StaticArraysMesh, J::Int)
    return StaticElement(
        mesh.centers[J],
        hcat(mesh.vertices[mesh.faces[J]]...)',  # makes a 4x3 SMatrix
        mesh.normals[J],
        mesh.areas[J],
        mesh.radii[J],
    )
end

center(e::StaticElement) = e.center
normal(e::StaticElement) = e.normal
area(e::StaticElement) = e.area
radius(e::StaticElement) = e.radius
vertices(e::StaticElement) = e.vertices


########################
#  Reflected elements  #
########################

struct ReflectedElement{T}
    element::T
end
free_surface_symmetry(e) = ReflectedElement(e)

center(e::ReflectedElement{T}) where T <: Union{NamedTuple, LazyElement} = (c = center(e.element); [c[1], c[2], -c[3]])
normal(e::ReflectedElement{T}) where T <: Union{NamedTuple, LazyElement} = (n = normal(e.element); [n[1], n[2], -n[3]])
area(e::ReflectedElement) = area(e.element)
radius(e::ReflectedElement) = radius(e.element)
function vertices(e::ReflectedElement{T} where T <: Union{NamedTuple, LazyElement})
    v = vertices(e.element)
    ET = eltype(v)
    return ET[
        v[4, 1]  v[4, 2]  -v[4, 3];
        v[3, 1]  v[3, 2]  -v[3, 3];
        v[2, 1]  v[2, 2]  -v[2, 3];
        v[1, 1]  v[1, 2]  -v[1, 3];
        ]   # Inverting order such that order is still consistent with normal vector
end

center(e::ReflectedElement{StaticElement}) = (c = center(e.element); SVector(c[1], c[2], -c[3]))
normal(e::ReflectedElement{StaticElement}) = (n = normal(e.element); SVector(n[1], n[2], -n[3]))
area(e::ReflectedElement{StaticElement}) = area(e.element)
radius(e::ReflectedElement{StaticElement}) = radius(e.element)
function vertices(e::ReflectedElement{StaticElement})
    v = vertices(e.element)
    return @SMatrix [
        v[4, 1] v[4, 2] -v[4, 3]
        v[3, 1] v[3, 2] -v[3, 3];
        v[2, 1] v[2, 2] -v[2, 3];
        v[1, 1] v[1, 2] -v[1, 3];
    ]   # Inverting order such that order is still consistent with normal vector
end

