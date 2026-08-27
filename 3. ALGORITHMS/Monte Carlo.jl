using Combinatorics
using Graphs
using LaTeXStrings
using Plots
using Base.Threads
using Random
using Statistics
using DelimitedFiles
using CSV
using DataFrames


function MonteCarlo_epsilon_mean_std(A::AbstractMatrix{<:Real}, num_samples::Int) # This matrix may be squared and symmetric
    n = size(A, 1) # Number of vertices
    A_dense = Matrix{Float64}(A) # Convert the adjacency matrix to a dense format for faster access
    mean_epsilon = 0.0 # Initialize the mean
    mean_sqr_epsilon = 0.0 # Initialize the mean of the squares for variance calculation
    for _ in 1:num_samples
        p = randperm(n) # Generate a random permutation of the vertices
        s = 0.0

        @inbounds for i in 1:n-1
            @inbounds for j in i+1:n # Only consider the upper triangle of the matrix to avoid double counting
                    d = A_dense[i, j] - A_dense[p[i], p[j]] # Compute the difference for the (i, j) entry
                    s += 2*d^2 # Accumulate the squared differences and multiply by 2 to account for the lower triangle
            end
        end

        epsilon = s / (n * (n - 1)) # Compute the normalized frobenius distance for the current permutation
        mean_epsilon += epsilon # Accumulate the mean
        mean_sqr_epsilon += epsilon^2 # Accumulate the variance
     end
    
    mean_epsilon /= num_samples # Compute the mean
    std_epsilon = sqrt(mean_sqr_epsilon / num_samples - mean_epsilon^2) # Compute the standard deviation

    error_mean = std_epsilon /  sqrt(num_samples) # Standard error of the mean
    error_std = std_epsilon / sqrt(2 * num_samples) # Standard error of the variance (using the fact that the variance of the sample variance is approximately 2*std^4/n for large n)


    return mean_epsilon, std_epsilon, error_mean, error_std
end

function theoretical_mean_std(A::AbstractMatrix{<:Real})
    density = Graphs.density(SimpleGraph(A)) # Compute the density of the graph
    theoretical_mean = 2 * density * (1 - density) # Theoretical mean
    degrees = vec(sum(A, dims=2)) # Compute the degree of each vertex
    m = sum(degrees) / 2 # Total number of links in the graph
    n = size(A, 1) # Number of vertices
    
    theoretical_variance =  (16 / (n^2 * (n-1)^2))*(m * density * (1 - m*density) + (4 * (m * (m-1) - sum(degrees.*(degrees .- 1)))^2)/(n * (n-1) * (n-2) * (n-3)) + ((sum(degrees.*(degrees .- 1)))^2)/(n * (n-1) * (n-2))) # Theoretical variance
    
    
    

    return theoretical_mean, sqrt(theoretical_variance)
end


p = 0.0
num_runs = 20
num_samples = 100000
k=2

sizes = [100, 150, 200, 250, 300, 350, 400]
n_plot = collect(100:100:2000)
theoretical_means = Float64[]
theoretical_stds = Float64[]
for size in n_plot
    #k = 2 * round(Int, 0.101*(size-1)/2)
    #g = erdos_renyi(size, p, seed=1) # Erdos-Renyi graph 
    #g = watts_strogatz(size, k, p; seed=10) # Watts-Strogatz graph 
    g = static_scale_free(size, round(Int, 0.3*size*(size-1)/2), 2, seed=1)
    #g = barabasi_albert(size, k, seed=1) # Barabasi-Albert graph
    A = adjacency_matrix(g)
    theoretical_mean, theoretical_std = theoretical_mean_std(A)
    push!(theoretical_means, theoretical_mean)
    push!(theoretical_stds, theoretical_std)
end
data = hcat(n_plot, theoretical_means, theoretical_stds)
writedlm("Divergence/t_config_2.csv", data, ',')

"""
MC_means = zeros(Float64, length(sizes))
MC_stds = zeros(Float64, length(sizes))
MC_errors_mean = zeros(Float64, length(sizes))
MC_errors_std = zeros(Float64, length(sizes))

@threads for idx in eachindex(sizes)
    time_taken = @elapsed begin
        size = sizes[idx]
        #g = erdos_renyi(size, p, seed=1)
        #g = watts_strogatz(size, k, p; seed=10)
        #g = barabasi_albert(size, k, seed=1)
        g = static_scale_free(size, round(Int, 0.3*size*(size-1)/2), 2, seed=1)
        A = adjacency_matrix(g)
        mean_epsilon, std_epsilon, error_mean, error_std = MonteCarlo_epsilon_mean_std(A, num_samples)
        
        MC_means[idx] = mean_epsilon
        MC_stds[idx] = std_epsilon
        MC_errors_mean[idx] = error_mean
        MC_errors_std[idx] = error_std
    end
    println("Time taken to size $(size): ", time_taken, " seconds")
end
data = hcat(sizes, MC_means, MC_stds, MC_errors_mean, MC_errors_std)
writedlm("Divergence/mc_config_2.csv", data, ',')
"""
