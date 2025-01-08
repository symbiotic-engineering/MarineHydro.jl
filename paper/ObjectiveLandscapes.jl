using Plots
include("power.jl")

# Define the bounds
lower = [1.0, 1.0]
upper = [5.0, 5.0]


r0_range = collect(lower[1]:0.1:upper[1])  
dx0_range = collect(lower[2]:0.1:upper[2])  

z_values = [power(r0, dx0) for r0 in r0_range, dx0 in dx0_range]

contour(r0_range, dx0_range, log.(z_values), title="Power per volume Objective space",
     xlabel="r0", ylabel="dx0", fill=true)

    #objective for each iteration 
arrays = [100, 200, 300, 600, 800]

r0_each_iter = [3.1, 3.2, 3.3, 3.5, 3.6]
dx0_each_iter = [2.1, 2.3, 2.5, 2.7, 3.0]

# Plot the points and connect them with a line (trace)
scatter!(r0_trace, dx0_trace, color=:red, markersize=5)
plot!(r0_trace, dx0_trace, color=:blue, linewidth=2)

savefig("objective_contour.pdf")