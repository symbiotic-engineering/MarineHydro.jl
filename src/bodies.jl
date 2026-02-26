

struct FloatingBody
    mesh::Mesh
    dofs::Dict # {name::Str : dofs::AbstractMatrix (nfaces,3)}

    function FloatingBody(mesh::Mesh, dofs::Dict)
        #add assert statements

        return new(mesh, dofs)
    end
end



function FloatingBody(mesh::Mesh, rigid_dof_list::Vector{String}, rotation_center::AbstractVector)
    dofs = Dict()

    for dof_name in rigid_dof_list
        if dof_name in ["Surge", "Sway", "Heave"]
            dofs[dof_name] = translational_dofs(mesh, dof_name)
        elseif dof_name in ["Roll", "Pitch", "Yaw"]
            dofs[dof_name] = rotational_dofs(mesh, dof_name, rotation_center)
        end
    end
    return FloatingBody(mesh, dofs)
end

function translational_dofs(mesh::Mesh, dof_name::String)
    num_panels = mesh.nfaces
    dof = zeros(num_panels, 3)
    if dof_name=="Surge"
        dof[:,1] .= 1
    elseif dof_name=="Sway"
        dof[:,2] .= 1
    elseif dof_name=="Heave"
        dof[:,3] .= 1
    end
    return dof
end

function rotational_dofs(mesh::Mesh, dof_name::String, rotation_center::AbstractVector)
    face_centers = mesh.centers
    if dof_name=="Roll"
        axis_of_rot = [1, 0, 0]
    elseif dof_name=="Pitch"
        axis_of_rot = [0, 1, 0]
    elseif dof_name=="Yaw"
        axis_of_rot = [0, 0, 1]
    end
    pos_vec = face_centers .- rotation_center'
    dof_vecs = cross.(Ref(axis_of_rot), eachrow(pos_vec))
    dof = stack(dof_vecs)' # make vector of vectors into a matrix
    return dof
end
