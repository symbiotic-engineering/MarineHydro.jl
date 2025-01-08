using BEM
using Zygote
using FiniteDifferences
using Test
using ColorTypes
using Plots
include("MeshGradients_singlebody.jl")
include("BemProgram.jl")

orange = RGB(230/255,159/255,0/255) 
vermillion = RGB(213/255, 94/255, 0/255) 
bluishgreen = RGB(0/255, 158/255, 115/255) 

green_functions = (
        Rankine(),
        RankineReflected(),
        GFWu(),
    )
omega = 1.03
r = 1.0
dof = [0,0,1]
A(radius) = added_mass_program(radius,omega,dof)  
B(radius) = damping_program(radius,omega,dof)
@show A(1.0)
A_r_grad1, = Zygote.gradient(A,r)
fd_grad1 = FiniteDifferences.central_fdm(2, 1)(A, r)
@test A_r_grad1 !== nothing
@test typeof(A_r_grad1) == Float64
@test isapprox(A_r_grad1, fd_grad1, atol=1e-6,rtol=1e-6)


B_r_grad1, = Zygote.gradient(B,r)
fd_grad1 = FiniteDifferences.central_fdm(2, 1)(B, r)
@test A_r_grad1 !== nothing
@test typeof(B_r_grad1) == Float64
@test isapprox(B_r_grad1, fd_grad1, atol=1e-6,rtol=1e-6)

radius_range = [1,2,3,4,5]
r = 1.0
w = 1.03
heave = [0,0,1] #heave
surge = [1,0,0]

added_mass_program_fd(radius) = added_mass_program(radius,w,heave)
fd_julia = [central_fdm(3, 1)(added_mass_program_fd, r) for r in radius_range] #check w omega
_gradients = Array{Any}(undef, length(radius_range))

for (i, r) in enumerate(radius_range)
        grads = Zygote.gradient(X -> added_mass_program(X[1], X[2], X[3]), [r, w, heave])[1][1]
        _gradients[i] = grads
    end
println(fd_julia)
println(_gradients )
print(_gradients .- fd_julia )
plot(collect(radius_range), _gradients , xlabel="r(radius)", ylabel="∂A_∂r", label="AD", marker = "*",color = vermillion,linestyle = :dash)
plot!(collect(radius_range), fd_julia, xlabel="r(radius)", ylabel="∂A_∂r", label="FD",marker = "*",color = bluishgreen,linestyle = :solid)

savefig("/home/cornell/BEMJulia/BEM.jl/paper/Plots/fd_ad_added_mass_heave_with_radius.pdf")

#### damping ###
damping_fd(radius) = damping_program(radius,w,heave)
fd_julia = [central_fdm(3, 1)(damping_fd, r) for r in radius_range] #check w omega
_gradients = Array{Any}(undef, length(radius_range))

for (i, r) in enumerate(radius_range)
    grads = Zygote.gradient(X -> damping_program(X[1], X[2], X[3]), [r, w, heave])[1][1]
    _gradients[i] = grads
end
println(fd_julia)
println(_gradients )
print(_gradients .- fd_julia )
plot(collect(radius_range), _gradients , xlabel="r(radius)", ylabel="∂B_∂r", label="AD", marker = "*",color = vermillion,linestyle = :dash)
plot!(collect(radius_range), fd_julia, xlabel="r(radius)", ylabel="∂B_∂r", label="FD",marker = "*",color = bluishgreen,linestyle = :solid)

savefig("/home/cornell/BEMJulia/BEM.jl/paper/Plots/fd_ad_damping_heave_with_radius.pdf")



#                                 ####### with respect to omega ###################
radius = 1
omega_range  = collect(range(0.01, stop=4, step=0.4))
fd_julia = [central_fdm(3, 1)(omega_added_mass_bem_program, omega) for omega in omega_range] #check w omega
_gradients = Array{Any}(undef, length(omega_range))
for (i, omega) in enumerate(omega_range)
        grads = Zygote.gradient(X -> omega_added_mass_bem_program(X[1], X[2], X[3]), [omega,radius, heave])[1][1]
        _gradients[i] = grads
    end



plot(omega_range, fd_julia, xlabel="ω", ylabel="∂A_∂ω", label="FD",marker = "*",color = bluishgreen,linestyle = :solid)
plot!(omega_range, _gradients , xlabel="ω", ylabel="∂A_∂ω", label="AD", marker = "*",color = vermillion,linestyle = :dash)
savefig("/home/cornell/BEMJulia/BEM.jl/paper/Plots/fd_ad_added_mass_omega_heave.pdf")