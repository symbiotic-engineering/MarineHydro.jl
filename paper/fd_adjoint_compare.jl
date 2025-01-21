using Zygote
using BenchmarkTools
using FiniteDiff
using Plots
using ColorTypes

# Define colorblind-friendly colors using hex values
orange = RGB(230/255,159/255,0/255)  
vermillion = RGB(213/255, 94/255, 0/255) 
bluishgreen = RGB(0/255, 158/255, 115/255)  


include("./Power.jl")

function objective(x)
    p = -1*power(x[1], x[2])
    return p
end


# Dummy design variables
function augmented_objective(x)
    return objective(x[1:2])  # Only the first two variables affect the objective
end

function finite_diff_gradient(f, x)
    return FiniteDiff.finite_difference_gradient(f, x)
end

function zygote_gradient(f, x)
    grad, = Zygote.jacobian(f, x)
    return grad
end

@show finite_diff_gradient(objective, [1.0,5.0])
@show zygote_gradient(objective, [1.0,5.0])


function compare_gradients(x_dummy)
    zyg_grad = zygote_gradient(augmented_objective, x_dummy)
    fd_grad = finite_diff_gradient(augmented_objective, x_dummy)
    fd_benchmark = @benchmark finite_diff_gradient(augmented_objective, $x_dummy)
    zyg_benchmark = @benchmark zygote_gradient(augmented_objective, $x_dummy)
    return fd_benchmark, zyg_benchmark
end


# List of sizes for design variables
sizes = [10,40,100] 

# Store benchmark times
fd_times = Float64[]
zyg_times = Float64[]

# Run benchmarks for different sizes
for size in sizes
    x = [2.0, 3.0]  
    x_dummy = vcat(x, rand(size .- 2))
    fd_benchmark, zyg_benchmark = compare_gradients(x_dummy)
    #  median times for each method
    push!(fd_times, median(fd_benchmark).time / 1e9)  # Convert time to seconds
    push!(zyg_times, median(zyg_benchmark).time / 1e9)
end

# Plotting the results
plot(sizes, fd_times, label="Finite Difference", xlabel="Number of Design Variables", ylabel="Time (seconds)", lw=2, color = bluishgreen)
plot!(sizes, zyg_times, label="AD", xlabel="Number of Design Variables", ylabel="Time (seconds)", lw= 2, color = vermillion)
savefig("Time_AD_FD_.pdf")
