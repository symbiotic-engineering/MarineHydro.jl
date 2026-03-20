
using LinearAlgebra: cross, dot, norm

struct FloatingBody
    mesh::Mesh
    dofs::NamedTuple 
    body_name::String

    function FloatingBody(mesh::Mesh, dofs::NamedTuple, body_name::String)
        # add assert statements

        return new(mesh, dofs, body_name)
    end
end



function FloatingBody(mesh::Mesh, rigid_dof_list::Vector{String}, rotation_center::AbstractVector, body_name::String)
    
    # generator
    dof_pairs = (Symbol(name) => if name in ["Surge", "Sway", "Heave"]
            translational_dofs(mesh, name)
        else
            rotational_dofs(mesh, name, rotation_center)
        end for name in rigid_dof_list)
            
    # convert Pair to NamedTuple using ; and splat
    dofs = (; dof_pairs...)

    return FloatingBody(mesh, dofs, body_name)
end




# rigid_dof_list can contain symbols or strings 
function FloatingBody(mesh::Mesh, rigid_dof_list::Vector{Symbol}, rotation_center::AbstractVector, body_name::String)
    return FloatingBody(mesh, string.(rigid_dof_list), rotation_center, body_name)
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
    dof = copy(stack(dof_vecs)') # make vector of vectors into a matrix
    return dof
end

# If rotation center not specified, assume it is at origin.
function FloatingBody(mesh::Mesh, rigid_dof_list::Vector{String}, body_name::String)
    rotation_center = [0.0,0.0,0.0]
    for dof in rigid_dof_list
        if dof in ["Roll","Pitch","Yaw"]
            display("Setting origin as rotation center.")
        end
    end
    return FloatingBody(mesh, rigid_dof_list, rotation_center, body_name)
end

