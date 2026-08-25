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

function compute_orbits(first_bin::Matrix{Int32}, epsilon_bins::Vector{})
    bin_max = length(epsilon_bins)
    orbits = Vector{Vector{Vector{Int}}}(undef, bin_max)
    n = size(first_bin, 1)
    g = SimpleGraph(n)
    for k in 1:bin_max
        for j in 1:n
            for i in 1:j-1
                if first_bin[i,j] == k
                    add_edge!(g, i, j)
                end
            end
        end
        orbits[k] = connected_components(g)
    end
    return orbits
end

function epsilon_aut_size(A::AbstractMatrix{<:Real}) # This matrix may be whatever category of matrix, but it should be squared and symmetric.
    n = size(A, 1) # Number of vertices
    first_bin =  fill(typemax(Int32), n, n)
    A_dense = Matrix{Float64}(A) # Convert the adjacency matrix to a dense format for faster access
    epsilon_bins = [4 * m / (n * (n - 1)) for m in 0:(n * (n - 1) ÷ 4)] # Define bins for the histogram based on the possible values of epsilon
    histogram = zeros(Int, length(epsilon_bins)) # Initialize the histogram vector
    time_taken = @elapsed begin
        for p in permutations(1:n) # Generate all permutations of the vertices
            s = 0.0

            @inbounds for i in 1:n-1
               @inbounds for j in i+1:n # Only consider the upper triangle of the matrix to avoid double counting
                    d = A_dense[i, j] - A_dense[p[i], p[j]] # Compute the difference for the (i, j) entry
                    s += 2*d^2 # Accumulate the squared differences and multiply by 2 to account for the lower triangle
                end
            end

            epsilon = s / (n * (n - 1)) # Compute the normalized frobenius distance for this permutation
            histogram[searchsortedlast(epsilon_bins, epsilon)] += 1 # Increment the appropriate bin in the histogram
        end
    end    
    automorphism_size = cumsum(histogram) # Compute the cumulative sum to get the size of the automorphism group for each bin
    println("Time taken to compute epsilon vector: ", time_taken, " seconds")
    return epsilon_bins, automorphism_size, histogram  # Return the epsilon bins and the normalized automorphism sizes
end


num_samples = 1000

deg = CSV.read("zs/deg.csv", DataFrame, header=1)[:,1] # Read the degree distribution




@threads for i in 1:num_samples
    g = random_configuration_model(11, deg, seed=i) # Generate the configuration model
    A =adjacency_matrix(g)

    epsilon_bins, automorphism_size, _ = epsilon_aut_size(A)

    data = hcat(epsilon_bins, automorphism_size)
    writedlm("zs/rand$(i).csv", data, ',') 
end
