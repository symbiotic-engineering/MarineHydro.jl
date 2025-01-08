using Optim
using Zygote
using Plots
include("/home/cornell/BEMJulia/MarineHydro.jl/paper/Power.jl")
# Define the objective function
# Wrapper to combine inputs for optimization
function objective(x)
    p = -1*power(x[1], x[2])
    return p
end


function g!(G, x)
    grad, = Zygote.jacobian(objective, x)
    f = objective(x)
   # print("current:$f")
    for i in eachindex(G)
        G[i] = grad[i]
    end
end

# Initial guesses for r and dx
r0 = 2.0 
dx0 = 3.0   
x0 = vcat(r0, dx0)  
lower = [1.0, 1.0]
upper = [4.0, 4.0]
@show objective(x0)

# Storage for convergence data
obj_values = []  # To store objective values
x_history = []   # To store x at each iteration

# # Callback function to store intermediate results
# callback = function (state)
#     push!(obj_values, state.value)
#     push!(x_history, copy(state.x))
#     return false 
# end

# Perform the optimization
res = optimize(
    objective, 
    g!,
    lower, 
    upper, 
    x0, 
    Fminbox(LBFGS()), 
    Optim.Options(iterations=100, g_tol=1e-3, show_trace=true, store_trace=true)
)

# Display the results
println("Optimal x: ", Optim.minimizer(res))
println("Minimum objective value: ", Optim.minimum(res))
print(res)
