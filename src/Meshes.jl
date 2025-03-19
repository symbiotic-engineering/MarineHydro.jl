using PyCall
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



function Mesh(mesh::PyObject)
    mesh = Mesh(
        mesh.vertices, mesh.faces, 
        mesh.faces_centers, mesh.faces_normals, mesh.faces_areas, mesh.faces_radiuses,mesh.nb_vertices,mesh.nb_faces
    )
    return mesh
end


function element(mesh::Mesh, J::Int)
    return (
        center = mesh.centers[J, :],
        vertices = mesh.vertices[mesh.faces[J, :] .+ 1, :],
        radius = mesh.radii[J],
        normal = mesh.normals[J, :],
        area = mesh.areas[J],
        faces = mesh.faces[J, :] .+ 1
    )
end

center(nt::NamedTuple) = nt[:center]
vertices(nt::NamedTuple) = nt[:vertices]
radius(nt::NamedTuple) = nt[:radius]
normal(nt::NamedTuple) = nt[:normal]
area(nt::NamedTuple) = nt[:area]
faces(nt::NamedTuple) = nt[:faces]


struct LazyElement
    mesh::Mesh
    index::Integer
end

center(e::LazyElement) = e.mesh.centers[e.index, :] #, e.mesh.centers[e.index, 1], e.mesh.centers[e.index, 2]
normal(e::LazyElement) = e.mesh.normals[e.index, :] #, e.mesh.normals[e.index, 1], e.mesh.normals[e.index, 2]
area(e::LazyElement) = e.mesh.areas[e.index]
radius(e::LazyElement) = e.mesh.radii[e.index]
vertices(e::LazyElement) = e.mesh.vertices[e.mesh.faces[e.index, :] .+ 1, :]


struct ReflectedElement{T}
    element::T
end
free_surface_symmetry(e) = ReflectedElement(e)

center(e::ReflectedElement) = (c = center(e.element); [c[1], c[2], -c[3]])
normal(e::ReflectedElement) = (n = normal(e.element); [n[1], n[2], -n[3]])
area(e::ReflectedElement) = area(e.element)
radius(e::ReflectedElement) = radius(e.element)
function vertices(e::ReflectedElement)
    v = vertices(e.element)
    return [
        v[4, 1] v[4, 2] -v[4, 3]
        v[3, 1] v[3, 2] -v[3, 3];
        v[2, 1] v[2, 2] -v[2, 3];
        v[1, 1] v[1, 2] -v[1, 3];
    ]   # Inverting order such that order is still consistent with normal vector
end
