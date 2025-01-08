
using ChainRulesCore
using FiniteDifferences
using PyCall
using Zygote

cpt = pyimport("capytaine")
resolution = (6,6)


function cptMesh(radius)
        cptmesh = cpt.mesh_sphere(name="sphere", radius = radius,center = (0,0,0),resolution = resolution)
        cptmesh.keep_immersed_part(inplace = true)
    return cptmesh
end


function get_center_x(J,mesh)
    centerx = mesh.centers[J,1]
    return centerx
end


function mesh_sizes(radius) 
    cptmesh = cptMesh(radius)
    return cptmesh.nb_faces, cptmesh.nb_vertices
end


h = 1e-3

∂J_r_fd(f,r) = (f(r+h) .- f(r-h)) ./ (2*h)



function sphere_vertices(radius) #To Do : a non differetiable argument for resolution
    cptmesh = cptMesh(radius)
    return reduce(vcat,cptmesh.vertices)
end
function sphere_centers(radius)
    cptmesh = cptMesh(radius)

    return reduce(vcat,cptmesh.faces_centers)
end

function sphere_normals(radius)
    cptmesh = cptMesh(radius)
    return reduce(vcat,cptmesh.faces_normals)
end
function sphere_radii(radius)
    cptmesh = cptMesh(radius)
    return reduce(vcat,cptmesh.faces_radiuses)
end
function sphere_areas(radius)
    cptmesh = cptMesh(radius)
    return reduce(vcat,cptmesh.faces_areas)
end

function sphere_faces(radius)
    cptmesh = cptMesh(radius)
    return reduce(vcat,cptmesh.faces)
end

# # #for one body
# define_rrule_with_finite_differences(sphere_areas, [1,2])
# define_rrule_with_finite_differences(sphere_radii, [1,2])
# define_rrule_with_finite_differences(sphere_normals, [1,2])
# define_rrule_with_finite_differences(sphere_centers, [1,2])
# define_rrule_with_finite_differences(sphere_vertices, [1,2])
# define_rrule_with_finite_differences(sphere_faces, [1,2])




∂J_r_fd(f,r; h=1e-5) = (f(r+h) .- f(r-h)) ./ (2*h)
#central_fdm.(2, 1)(f, r) #2nd order Central FD

#sensitivities with respect to dx1...fast way to do mesh convergence?




"""Rules for Zygote finite differencing of the mesh"""

function ChainRulesCore.rrule(::typeof(sphere_vertices), r)
    y = sphere_vertices(r)
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function
        dv = ∂J_r_fd(sphere_vertices,r)  # Finite difference derivative with respect to `r`
        dx = dv' * dy    #pullback
        return (df, dx)  
    end
    
    return y, f_pullback
end


#   function ChainRulesCore.frule(::typeof(sphere_vertices), r) #zygote by default does not use frule
#     y = NoTangent()
#     dx = ∂J_r_fd(sphere_vertices,r,dx1)' * dy

#     return y, dx
#   end



function ChainRulesCore.rrule(::typeof(sphere_centers), r)
    y = sphere_centers(r)
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function
        dv = ∂J_r_fd(sphere_centers, r)  # Finite difference derivative with respect to `r`
        dx = dv' * dy    #pullback
        return (df, dx)  
    end
    return y, f_pullback
  end


function ChainRulesCore.rrule(::typeof(sphere_normals), r)
    y = sphere_normals(r)
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function
        dv = ∂J_r_fd(sphere_normals, r)  # Finite difference derivative with respect to `r`
        dx = dv' * dy    #pullback
        return (df, dx)  
    end
    return y, f_pullback
  end


function ChainRulesCore.rrule(::typeof(sphere_radii), r)
    y = sphere_radii(r)
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function
        dv = ∂J_r_fd(sphere_radii, r)  # Finite difference derivative with respect to `r`
        dx = dv' * dy    #pullback
        return (df, dx)  
    end
    return y, f_pullback
  end




function ChainRulesCore.rrule(::typeof(sphere_areas), r)
    y = sphere_areas(r)
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function
        dv = ∂J_r_fd(sphere_areas, r)  # Finite difference derivative with respect to `r`
        dx = dv' * dy    #pullback
        return (df, dx)  
    end
    return y, f_pullback
  end



function ChainRulesCore.rrule(::typeof(sphere_faces), r)
    y = sphere_faces(r)
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function
        dv = ∂J_r_fd(sphere_faces, r)  # Finite difference derivative with respect to `r`
        dx = dv' * dy    #pullback
        return (df, dx)  
    end
    return y, f_pullback
  end



  function differentiableMesh(r1)
    nf,nv =   Zygote.@ignore mesh_sizes(r1) #hack..need to know how many first to reshape
    vertices = reshape(sphere_vertices(r1),(nv,3))
    faces =  reshape(sphere_faces(r1),(nf,4))
    centers = reshape(sphere_centers(r1),(nf,3))
    normals = reshape(sphere_normals(r1),(nf,3))
    areas = sphere_areas(r1)
    radii = sphere_radii(r1)
    nvertices = size(vertices,1)
    nfaces = size(centers,1)
    return Mesh(vertices,faces,centers,normals,areas,radii,nvertices,nfaces) #
end


  # test differentiability
# function func(r1,dx1)
#     m = differentiableMesh(r1,dx1)
#     area = m.areas
#     normal = m.normals
#     return sum(area .* normal)
# end

# Zygote.jacobian(x->func(x[1],x[2]),[2.0,3.0])