using Combinatorics
using Graphs
using LaTeXStrings
using Plots
using Base.Threads
using Random
Random.seed!(80)
using Statistics
using DelimitedFiles
using CSV
using DataFrames

const δ = 0.525
const N_TRANS = 3000
const N_MEASURE = 1000
const R = 50
const M = 1000
# const β = 0.72pi


# =====================================================================
# Función para cargar órbitas desde un archivo CSV.
# =====================================================================

function load_orb(path::String)
    orbits_by_eps = Dict{Float64, Vector{Vector{Int}}}()
    current_eps = nothing
    current_orbits = Vector{Vector{Int}}()

    for line in eachline(path)
        line = strip(line)
        isempty(line) && continue

        if startswith(line, "epsilon")
            if current_eps !== nothing
                orbits_by_eps[current_eps] = current_orbits
            end
            current_eps = parse(Float64, split(line, "=")[2])
            current_orbits = Vector{Vector{Int}}()
        else
            orbit = parse.(Int, split(line, ","))
            push!(current_orbits, orbit)
        end

        if current_eps !== nothing
            orbits_by_eps[current_eps] = current_orbits
        end
    end
    return orbits_by_eps
end

# =====================================================================
# Función para cargar la matriz de adyacencia un- o weighted desde un CSV con los edges, sin saber el nombre de las columnas.
# =====================================================================

function load_A(path::String, weighted::Bool, zero_indexed::Bool, directed::Bool)
    df = CSV.read(path, DataFrame)
    rename!(df, Symbol.(strip.(string.(names(df)))))
    
    sources = collect(df[!, 1])
    targets = collect(df[!, 2])
    weights = weighted ? collect(Float64.(df[!, 3])) : ones(Float64, length(sources))

    if zero_indexed
        sources .+= 1
        targets .+= 1
    end

    n = maximum(vcat(sources, targets))
    A = zeros(Float64, n, n)

    for k in eachindex(sources)
        i, j, w = sources[k], targets[k], weights[k]
        i == j && continue
        A[i, j] = w
        if !directed
            A[j, i] = w
        end
    end

    return A
end


# =====================================================================
# Map: x_i^{t+1} = [ beta*I(x_i^t) + sigma * sum_j A_ij I(x_j^t) + delta ] mod 2pi
#   I(x) = (1 - cos x) / 2
# =====================================================================


function map_step!(x::Vector{Float64}, A::Matrix{Float64}, β::Float64, σ::Float64, δ::Float64)
    n = length(x)
    I_x = (1 .- cos.(x)) ./ 2
    coupling = A * I_x
    for i in 1:n
        x[i] = mod(β * I_x[i] + σ * coupling[i] + δ, 2pi)
    end
    return x
end

# =====================================================================
# wrap: envuelve una diferencia angular a (-pi, pi]
# =====================================================================

wrap(x::Float64) = mod(x + π, 2π) - π

# =====================================================================
# Δx_{RMS} por órbita
# =====================================================================

function orbit_rms(hist::Matrix{Float64}, orbit::Vector{Int})
    T = size(hist, 2); nk = length(orbit)
    acc = 0.0
    @inbounds for t in 1:T
        sr = 0.0; si = 0.0
        for i in orbit
            sr += cos(hist[i, t]); si += sin(hist[i, t])
        end
        xb = atan(si, sr)
        s = 0.0
        for i in orbit
            d = wrap(hist[i, t] - xb); s += d * d
        end
        acc += s / nk
    end
    return sqrt(acc / T)
end

# =====================================================================
# Δx_{TSE}: media de Δx_RMS sobre las órbitas no triviales
# =====================================================================

function tse(hist::Matrix{Float64}, parts::Vector{Vector{Int}})
    tot = 0.0; nt = 0
    for orb in parts
        length(orb) == 1 && continue
        nt += 1; tot += orbit_rms(hist, orb)
    end
    return nt == 0 ? NaN : tot / nt
end

# =====================================================================
# Particiones de órbitas en el null model
# =====================================================================

function null_partition(orbits::Vector{Vector{Int}}, rng=Random.default_rng())
    nodes = vcat(orbits...)
    sizes =length.(orbits)
    shuffled = shuffle(rng, nodes)
    null_orbits = Vector{Vector{Int}}()
    idx = 1
    for s in sizes
        push!(null_orbits, shuffled[idx:idx+s-1])
        idx += s
    end
    return null_orbits
end

# =====================================================================
# safe_mean / safe_std: como mean/std pero devuelven NaN en vez de
# reventar si el vector (tras quitar NaNs) queda vacío o con 1 elemento.
# =====================================================================
 
safe_mean(v) = (f = filter(!isnan, v); isempty(f) ? NaN : mean(f))
safe_std(v)  = (f = filter(!isnan, v); length(f) < 2 ? NaN : std(f))





# =====================================================================
# Uso del código
# =====================================================================

const path_orbits = "Karate/orb_0.csv"
const path_edges = "Karate/edges.csv"

const β_values = [0.1π, 0.2π, 0.3π, 0.4π, 0.5π, 0.6π, 0.7π, 0.8π, 0.9π, 1.0π, 1.1π, 1.2π, 1.3π, 1.4π, 1.5π]
const σ_vals = [0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.2, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24]
function main()
    orbits_by_eps = load_orb(path_orbits)
    A = load_A(path_edges, true, false, false)
    n = size(A, 1)
    epsilons = sort(collect(keys(orbits_by_eps)))
    n_eps = length(epsilons)
    outdir = "sync/Karate"
    mkpath(outdir)
 
    @threads for σ_idx in eachindex(σ_vals)
            σ = σ_vals[σ_idx]
            rng = Xoshiro(1000 * σ_idx)
            println("Processing σ = $σ")
        for β in β_values
    
            real_tse_mat = fill(NaN, n_eps, R)   # TSE real, por epsilon y realización
            null_tse_mat = fill(NaN, n_eps, R)   # media del null, por epsilon y realización
            z_mat        = fill(NaN, n_eps, R)   # z-score, por epsilon y realización
    
            for r in 1:R
                x = 2π .* rand(rng, n)
                x_hist = zeros(Float64, n, N_MEASURE)
    
                for t in 1:(N_TRANS + N_MEASURE)
                    map_step!(x, A, β, σ, δ)
                    if t > N_TRANS
                        x_hist[:, t - N_TRANS] = x
                    end
                end
    
                for (i, ϵ) in enumerate(epsilons)
                    orbits = orbits_by_eps[ϵ]
                    real_val = tse(x_hist, orbits)
                    real_tse_mat[i, r] = real_val
                    isnan(real_val) && continue
    
                    null_vals = Vector{Float64}(undef, M)
                    for j in 1:M
                        null_orbits = null_partition(orbits, rng)
                        null_vals[j] = tse(x_hist, null_orbits)
                    end
                    filter!(!isnan, null_vals)
                    isempty(null_vals) && continue
    
                    μ_null = mean(null_vals)
                    σ_null = std(null_vals)
                    null_tse_mat[i, r] = μ_null
                    σ_null > 0 || continue
    
                    z_mat[i, r] = (real_val - μ_null) / σ_null
                end
            end
    
           
            mean_real = [safe_mean(real_tse_mat[i, :]) for i in 1:n_eps]
            std_real  = [safe_std(real_tse_mat[i, :])  for i in 1:n_eps]
            mean_null = [safe_mean(null_tse_mat[i, :]) for i in 1:n_eps]
            std_null = [safe_std(null_tse_mat[i, :]) for i in 1:n_eps]
            mean_z    = [safe_mean(z_mat[i, :])        for i in 1:n_eps]
            std_z     = [safe_std(z_mat[i, :])         for i in 1:n_eps]
    
            df = DataFrame(
                epsilon       = epsilons,
                n_orbits      = [length(orbits_by_eps[ϵ]) for ϵ in epsilons],
                mean_real_tse = mean_real,
                std_real_tse  = std_real,
                mean_null_tse = mean_null,
                std_null_tse  = std_null,
                mean_z        = mean_z,
                std_z         = std_z
            )
            outfile = joinpath(outdir, "pecora_map_sigma_$(σ)_beta_$(β).csv")
            CSV.write(outfile, df)
            println("σ = $σ, β = $β listo -> $outfile")
        end
    end
end
 
main()

            

