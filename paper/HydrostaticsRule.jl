using Zygote
using ChainRulesCore
using MarineHydro 
cpt = pyimport("capytaine")

∂J_r_fd_hs(f,r; h=1e-5) = (f(r+h) .- f(r-h)) ./ (2*h)

function hydrostatics(r)
    cptmesh1 = cpt.mesh_sphere(name="sphere", radius = r,resolution=(6,6),center = (0,0,0))
    cptbody = cpt.FloatingBody(cptmesh1,center_of_mass = (0,0,-0.1))
    cptbody.rotation_center=(0, 0, -0.1)
    cptbody.add_translation_dof(name="Heave")
    cptbody.keep_immersed_part()
    C1 = cptbody.compute_hydrostatics()["hydrostatic_stiffness"].sel(radiating_dof = "Heave").item()
    return C1
end

function ChainRulesCore.rrule(::typeof(hydrostatics), r)
    y = hydrostatics(r)
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function
        dv = ∂J_r_fd_hs(hydrostatics, r)  # Finite difference derivative with respect to `r`
        dx = dv' * dy    #pullback
        return (df, dx)  
    end
    return y, f_pullback
end

function inertia(r)
    cptmesh1 = cpt.mesh_sphere(name="sphere",  radius = r,resolution=(6,6),center = (0,0,0))
    cptbody = cpt.FloatingBody(cptmesh1,center_of_mass = (0,0,-0.3))
    cptbody.rotation_center=(0, 0, -0.3)
    cptbody.add_translation_dof(name="Heave")
    cptbody.keep_immersed_part()
    I1 = cptbody.compute_hydrostatics()["inertia_matrix"]
    return I1
end

function ChainRulesCore.rrule(::typeof(inertia), r)
    y = inertia(r)
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function
        dv = ∂J_r_fd_hs(inertia, r)  # Finite difference derivative with respect to `r`
        dx = dv' * dy    #pullback
        return (df, dx)  
    end
    
    return y, f_pullback
end





