using Combinatorics
using Graphs
using Random
using Statistics
using DelimitedFiles
using CSV
using DataFrames
using Base.Threads
 

# ---------------------------------------
# ÓRBITAS PARA CADA BIN
# ---------------------------------------
 
function compute_orbits(first_bin_0::Matrix{Int32}, n::Int, bin_max::Int)
    orbits = Vector{Vector{Vector{Int}}}(undef, bin_max)
 
    g = SimpleGraph(n)
    for k in 1:bin_max
        for j in 1:n
            for i in 1:j-1
                if first_bin_0[i,j] == k
                    add_edge!(g, i, j)
                end
            end
        end
        orbits[k] = connected_components(g)
    end
    return orbits
end
 
# ---------------------------------------
# COTA INFERIOR RESTANTE POR GRADO
# Misma logica que la version para simple graph,
# pero ahora se utiliza el grado de salida.
# ---------------------------------------
 
@inline function cota_inferior_restante_grado(grado_ord::Vector{Int}, usado::BitVector, cand_global::Vector{Int},
    k::Int, n::Int)::Int
 
    cota_k = 0
    j_idx = 1
 
    for i in k:n
        while j_idx <= n && usado[cand_global[j_idx]]
            j_idx += 1
        end
        j_idx > n && break
 
        cota_k += abs(grado_ord[i] - grado_ord[cand_global[j_idx]])
        j_idx += 1
    end
 
    return cota_k
end
 
# ---------------------------------------
# PODA POR ADYACENCIA (DIRIGIDA) Y POR GRADOS
#
# Ahora, en cada asignacion pi(k) = v, comparamos DOS entradas de la
# matriz por cada nodo i ya asignado, no una sola:
#   - Aw[k, i]  vs  Aw[v, pi(i)]   (arista k -> i)
#   - Aw[i, k]  vs  Aw[pi(i), v]   (arista i -> k)
# Cada discrepancia suma +1 (no +2 como en la version simetrica).
# ---------------------------------------
 
function poda_recursiva!(Aw::Matrix{Int8}, grado_ord::Vector{Int}, cand_nodo::Vector{Vector{Int}},
    cost_max::Int, n::Int, pi_vec::Vector{Int}, usado::BitVector, k::Int, hist::Vector{Int},
    perm_complete::Base.RefValue{Int}, perm_prunned::Base.RefValue{Int}, first_bin::Matrix{Int32},
    coste_adyacencia_parcial::Int, coste_grados_parcial::Int, cand_global::Vector{Int})
 
 
     if coste_grados_parcial + cota_inferior_restante_grado(grado_ord, usado, cand_global, k, n) > cost_max
        perm_prunned[] += 1
        return
     end
 
        if k > n
            perm_complete[] += 1
            bin_idx = (coste_adyacencia_parcial ÷ 2) + 1   # <- arreglo
            hist[bin_idx] += 1
            @inbounds for i in 1:n
                j = pi_vec[i]
                if bin_idx < first_bin[i,j]
                    first_bin[i,j] = bin_idx
                end
            end
            return
        end
 
     @inbounds for v in cand_nodo[k]
        usado[v] && continue
        coste_k = abs(grado_ord[k] - grado_ord[v])
        nuevo_coste_grados_parcial = coste_grados_parcial + coste_k
        if nuevo_coste_grados_parcial > cost_max
            perm_prunned[] += 1
            break
        end
 
        coste_A = 0
        broke = false
        @inbounds for i in 1:k-1
            j = pi_vec[i]
 
            # Direccion k -> i
            if Aw[k, i] != Aw[v, j]
                coste_A += 1
                if coste_adyacencia_parcial + coste_A > cost_max
                    broke = true
                    break
                end
            end
 
            # Direccion i -> k 
            if Aw[i, k] != Aw[j, v]
                coste_A += 1
                if coste_adyacencia_parcial + coste_A > cost_max
                    broke = true
                    break
                end
            end
        end
 
        nuevo_coste_adyacencia_parcial = coste_adyacencia_parcial + coste_A
        if broke
            perm_prunned[] += 1
            continue
        end
 
        pi_vec[k] = v
        usado[v] = true
 
        poda_recursiva!(Aw, grado_ord, cand_nodo, cost_max, n, pi_vec,
            usado, k+1, hist, perm_complete, perm_prunned, first_bin,
            nuevo_coste_adyacencia_parcial, nuevo_coste_grados_parcial, cand_global)
 
        usado[v] = false
 
    end
    return
end
 

 
function epsilon_aut(A, epsilon_max)
    n = size(A, 1)
 
    
    grado_out = vec(sum(A, dims=2))
    grado_in  = vec(sum(A, dims=1))
    grado = grado_out   # En lugar de ordenar por grado total, ordenamos por grado de salida
 
    nodos_ordenados = sortperm(grado, rev = true)
    Aw = Matrix{Int8}(A[nodos_ordenados, nodos_ordenados])
    grado_ord = grado[nodos_ordenados]
 

    nbins = n * (n - 1) ÷ 2 + 1
    epsilon_bins = [2m / (n * (n - 1)) for m in 0:(n*(n-1)÷2)]
    cost_max::Int = floor(Int, epsilon_max * n * (n - 1) + 1e-10)
    bin_max::Int = (cost_max ÷ 2) + 1
 
    cand_global = sortperm(grado_ord, rev = true)
 
    cand_nodo = Vector{Vector{Int}}(undef, n)
    for k in 1:n
        cand_nodo[k] = sortperm(abs.(grado_ord .- grado_ord[k]))
    end
 
   
    pre_prunned = 0
    spawn_pairs = Tuple{Int, Int, Int, Int}[]
    for v1 in 1:n
        r1 = abs(grado_ord[v1] - grado_ord[1])
        if r1 > cost_max
            pre_prunned += (n-1)
            continue
        end
        for v2 in 1:n
            v2 == v1 && continue
            r2 = abs(grado_ord[v2] - grado_ord[2])
            if r1 + r2 > cost_max
                pre_prunned += 1
                continue
            end
            coste_A = 0
            if Aw[1, 2] != Aw[v1, v2]
                coste_A += 1
            end
            if Aw[2, 1] != Aw[v2, v1]
                coste_A += 1
            end
            if coste_A > cost_max
                pre_prunned += 1
                continue
            end
            push!(spawn_pairs, (v1, v2, coste_A, r1 + r2))
        end
    end
    println("$pre_prunned Pre-prunned branches of $(n*(n-1))")
    n_tasks = length(spawn_pairs)
    task_order = sortperm(spawn_pairs, by = p -> (p[3], p[4]))
 
    T = Threads.nthreads()
    println("=" ^ 50)
    println("\n Spawning $n_tasks tasks across $T threads")
    println("=" ^ 50)
 
    hist_t = [zeros(Int, nbins) for _ in 1:T]
    perm_prunned_t = [Ref(0) for _ in 1:T]
    perm_complete_t = [Ref(0) for _ in 1:T]
    first_bin_t = [fill(typemax(Int32), n, n) for _ in 1:T]
 
    free_buf = Channel{Int}(T)
    for t in 1:T
        put!(free_buf, t)
    end
 
    time_taken = @elapsed begin
        Threads.@threads :dynamic for ii in 1:n_tasks
            idx = task_order[ii]
            v1, v2, coste_A, coste_grados = spawn_pairs[idx]
            bid = take!(free_buf)
            try
                pi_vec_l = Vector{Int}(undef, n)
                usado_l = falses(n)
                pi_vec_l[1] = v1; usado_l[v1] = true
                pi_vec_l[2] = v2; usado_l[v2] = true
                poda_recursiva!(Aw, grado_ord, cand_nodo, cost_max, n,
                    pi_vec_l, usado_l, 3, hist_t[bid], perm_complete_t[bid], perm_prunned_t[bid], first_bin_t[bid],
                    coste_A, coste_grados, cand_global)
            finally
                put!(free_buf, bid)
            end
        end
    end
 
    hist = zeros(Int, nbins)
    perm_prunned = 0
    perm_complete = 0
    first_bin = fill(typemax(Int32), n, n)
 
    @inbounds for t in 1:T
        hist .+= hist_t[t]
        perm_prunned += perm_prunned_t[t][]
        perm_complete += perm_complete_t[t][]
        for j in 1:n
            for i in 1:j-1
                if first_bin_t[t][i,j] < first_bin[i,j]
                    first_bin[i,j] = first_bin_t[t][i,j]
                end
            end
        end
    end
    aut_eps = cumsum(hist)
    perm_prunned += pre_prunned
 
    first_bin_0 = fill(typemax(Int32), n, n)
    @inbounds for j in 2:n
        for i in 1:j-1
            b = first_bin[i,j]
            if b < typemax(Int32)
                i_0, j_0 = nodos_ordenados[i] < nodos_ordenados[j] ? (nodos_ordenados[i], nodos_ordenados[j]) : (nodos_ordenados[j], nodos_ordenados[i])
                if b < first_bin_0[i_0, j_0]
                    first_bin_0[i_0, j_0] = b
                end
            end
        end
    end
 
    orbitas = compute_orbits(first_bin_0, n, bin_max)
 
    println("\nTotal complete permutations: $perm_complete")
    println("Total prunned permutations: $perm_prunned")
    println("Elapsed time: $(round(time_taken, digits=2)) seconds")
 
    return epsilon_bins, aut_eps, orbitas
end
