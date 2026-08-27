using Combinatorics
using Graphs
using Random
using Statistics
using DelimitedFiles
using CSV
using DataFrames
using Base.Threads
using LinearAlgebra

# ---------------------------------------
# ÓRBITAS POR ϵ
# ---------------------------------------

function compute_orbits_w(first_energy_0::Matrix{Float64}, n::Int, energias::Vector{Float64})
    m = length(energias)
    orbits = Vector{Vector{Vector{Int}}}(undef, m)

    g = SimpleGraph(n)
    current = connected_components(g)
    for (idx, e) in enumerate(energias)
        changed = false
        for j in 1:n
            for i in 1:j-1
                if first_energy_0[i, j] == e
                    changed |= add_edge!(g, i, j)
                end
            end
        end
        changed && (current = connected_components(g))
        orbits[idx] = current
    end
    return orbits
end

# ---------------------------------------
# COTA POR CANDIDATO
#
# Para dos filas W[k,:] y W[v,:], de entre todas las formas de
# emparejar sus n entradas, la que minimiza la suma de diferencias al
# cuadrado es emparejarlas ambas ordenadas en el mismo sentido:
#
#   sum_j (W[k,j]-W[v,tau(j)])^2  >=  sum_k (a_(k) - b_(k))^2
#
# donde a_(k), b_(k) son las entradas de las filas k y v ordenadas de
# mayor a menor. Esta cota es MAS AJUSTADA que la basada solo en la
# fuerza total (s_k-s_v)^2/n, porque usa la distribucion completa de
# pesos de cada fila, no solo su suma. Para grafos binarios, esta
# misma cota colapsa exactamente a |k_v-k_v'|,recuperando
# el caso no ponderado como caso particular.
# ---------------------------------------

@inline function coste_fila_ordenada(sorted_rows::Vector{Vector{Float64}}, k::Int, v::Int)::Float64
    a = sorted_rows[k]
    b = sorted_rows[v]
    s = 0.0
    @inbounds for idx in eachindex(a)
        d = a[idx] - b[idx]
        s += d * d
    end
    return s
end

# ---------------------------------------
# COTA INFERIOR RESTANTE POR FUERZA
#
# En lugar de utilizar la cota basada en la rearrengement
# inequality, que requiere O(n^3) para cada nodo del árbol 
# de búsqueda, se utiliza una cota más barata basada en la
# fuerza total de los nodos: 
#
# (s_k-s_v)^2/n,
#
# nacida de la desigualdad de Cauchy-Schwarz. Esta cota es más
# débil que la de rearrengement, pero es mucho más barata de
# evaluar (O(n) en lugar de O(n^3)).
# ---------------------------------------

@inline function cota_inferior_restante_fuerza(fuerza_ord::Vector{Float64}, usado::BitVector, cand_global::Vector{Int},
    k::Int, n::Int)::Float64

    cota_k = 0.0
    j_idx = 1

    for i in k:n
        while j_idx <= n && usado[cand_global[j_idx]]
            j_idx += 1
        end
        j_idx > n && break

        cota_k += (fuerza_ord[i] - fuerza_ord[cand_global[j_idx]])^2 / n
        j_idx += 1
    end

    return cota_k
end

# ---------------------------------------
# PODA POR ADYACENCIA Y POR FUERZA/FILA-ORDENADA 
# ---------------------------------------

function poda_recursiva_ponderada!(Ww::Matrix{Float64}, fuerza_ord::Vector{Float64},
    sorted_rows::Vector{Vector{Float64}}, cand_nodo::Vector{Vector{Int}},
    cost_max::Float64, n::Int, pi_vec::Vector{Int}, usado::BitVector, k::Int, energias::Vector{Float64},
    perm_complete::Base.RefValue{Int}, perm_prunned::Base.RefValue{Int}, first_energy::Matrix{Float64},
    coste_adyacencia_parcial::Float64, coste_fuerza_parcial::Float64, cand_global::Vector{Int},
    max_complete::Union{Int, Nothing} = nothing)


    if coste_fuerza_parcial + cota_inferior_restante_fuerza(fuerza_ord, usado, cand_global, k, n) > cost_max
        perm_prunned[] += 1
        return
    end

    if k > n
        perm_complete[] += 1
        if max_complete !== nothing && perm_complete[] > max_complete
            throw(ErrorException("max_complete=$max_complete permutaciones completas superado en un hilo; reduce ε_max"))
        end
        energia = coste_adyacencia_parcial / (n * (n - 1))
        push!(energias, energia)
        @inbounds for i in 1:n
            j = pi_vec[i]
            if energia < first_energy[i, j]
                first_energy[i, j] = energia
            end
        end
        return
    end

    @inbounds for v in cand_nodo[k]
        usado[v] && continue

        coste_k = coste_fila_ordenada(sorted_rows, k, v)
        nuevo_coste_fuerza_parcial = coste_fuerza_parcial + coste_k
        if nuevo_coste_fuerza_parcial > cost_max
            perm_prunned[] += 1
            break
        end

        coste_A = 0.0
        broke = false
        @inbounds for i in 1:k-1
            d = Ww[k, i] - Ww[v, pi_vec[i]]
            if d != 0.0
                coste_A += 2 * d^2
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

        poda_recursiva_ponderada!(Ww, fuerza_ord, sorted_rows, cand_nodo, cost_max, n, pi_vec,
            usado, k+1, energias, perm_complete, perm_prunned, first_energy,
            nuevo_coste_adyacencia_parcial, nuevo_coste_fuerza_parcial, cand_global,
            max_complete)

        usado[v] = false

    end
    return
end



function epsilon_aut_ponderado(W::AbstractMatrix{<:Real}, epsilon_max::Float64;
    max_complete::Union{Int, Nothing} = nothing)
    n = size(W, 1)
    @assert size(W, 2) == n "W debe ser cuadrada"
    @assert issymmetric(W) "W debe ser simétrica (grafo no dirigido)"

    fuerza = vec(sum(W, dims=2))

    nodos_ordenados = sortperm(fuerza, rev = true)
    Ww = Matrix{Float64}(W[nodos_ordenados, nodos_ordenados])
    fuerza_ord = fuerza[nodos_ordenados]


    sorted_rows = [sort(Ww[i, :], rev=true) for i in 1:n]

    cost_max::Float64 = epsilon_max * n * (n - 1)

    # cand_global ordena por fuerza: es el que usa
    # cota_inferior_restante_fuerza (nivel "rama", agregado sobre TODAS las
    # posiciones aun no decididas), donde solo la cota barata por fuerza es
    # factible de evaluar en cada nodo del arbol de busqueda (la version
    # ajustada por fila supondria resolver un problema de asignacion optima
    # -- algoritmo hungaro, O(n^3) -- en cada nodo, inviable).
    cand_global = sortperm(fuerza_ord, rev = true)

    # cand_nodo[k], en cambio, ordena los candidatos para UNA posicion k ya
    # fijada frente a UN candidato v ya fijado (nivel "candidato"): aqui SI
    # es barato evaluar la cota ajustada coste_fila_ordenada(k,v) para los n
    # candidatos (cada evaluacion es O(n) porque sorted_rows ya esta
    # ordenado), asi que ordenamos directamente por esa metrica -- la MISMA
    # que se compara contra cost_max en el bucle de poda_recursiva_ponderada!.
    # Esto es imprescindible: como la metrica de corte (coste_fila_ordenada)
    # NO es una funcion monotona de |Δfuerza| (dos filas pueden tener la
    # misma fuerza total y formas muy distintas), ordenar por |Δfuerza| como
    # antes invalidaba el "break" del bucle -- podia cortar la iteracion
    # antes de llegar a un candidato con coste_fila_ordenada bajo, perdiendo
    # permutaciones validas. Ordenando por la propia metrica de corte, el
    # "break" vuelve a ser valido (monotonia por construccion).
    cand_nodo = Vector{Vector{Int}}(undef, n)
    for k in 1:n
        costes_k = [coste_fila_ordenada(sorted_rows, k, v) for v in 1:n]
        cand_nodo[k] = sortperm(costes_k)
    end

    pre_prunned = 0
    spawn_pairs = Tuple{Int, Int, Float64, Float64}[]
    for v1 in 1:n
        r1 = coste_fila_ordenada(sorted_rows, 1, v1)
        if r1 > cost_max
            pre_prunned += (n-1)
            continue
        end
        for v2 in 1:n
            v2 == v1 && continue
            r2 = coste_fila_ordenada(sorted_rows, 2, v2)
            if r1 + r2 > cost_max
                pre_prunned += 1
                continue
            end
            d = Ww[1, 2] - Ww[v1, v2]
            coste_A = 2 * d^2
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

    energias_t = [Float64[] for _ in 1:T]
    perm_prunned_t = [Ref(0) for _ in 1:T]
    perm_complete_t = [Ref(0) for _ in 1:T]
    first_energy_t = [fill(Inf, n, n) for _ in 1:T]

    free_buf = Channel{Int}(T)
    for t in 1:T
        put!(free_buf, t)
    end

    time_taken = @elapsed begin
        Threads.@threads :dynamic for ii in 1:n_tasks
            idx = task_order[ii]
            v1, v2, coste_A, coste_fuerza = spawn_pairs[idx]
            bid = take!(free_buf)
            try
                pi_vec_l = Vector{Int}(undef, n)
                usado_l = falses(n)
                pi_vec_l[1] = v1; usado_l[v1] = true
                pi_vec_l[2] = v2; usado_l[v2] = true
                poda_recursiva_ponderada!(Ww, fuerza_ord, sorted_rows, cand_nodo, cost_max, n,
                    pi_vec_l, usado_l, 3, energias_t[bid], perm_complete_t[bid], perm_prunned_t[bid], first_energy_t[bid],
                    coste_A, coste_fuerza, cand_global, max_complete)
            finally
                put!(free_buf, bid)
            end
        end
    end

    all_energies = Float64[]
    perm_prunned = 0
    perm_complete = 0
    first_energy = fill(Inf, n, n)

    @inbounds for t in 1:T
        append!(all_energies, energias_t[t])
        perm_prunned += perm_prunned_t[t][]
        perm_complete += perm_complete_t[t][]
        for j in 1:n
            for i in 1:j-1
                if first_energy_t[t][i,j] < first_energy[i,j]
                    first_energy[i,j] = first_energy_t[t][i,j]
                end
            end
        end
    end
    perm_prunned += pre_prunned

    sort!(all_energies)

    epsilon_bins = Float64[]
    aut_eps = Int[]
    running = 0
    N = length(all_energies)
    i = 1
    while i <= N
        e = all_energies[i]
        j = i
        while j <= N && all_energies[j] == e
            j += 1
        end
        running += (j - i)
        push!(epsilon_bins, e)
        push!(aut_eps, running)
        i = j
    end

    first_energy_0 = fill(Inf, n, n)
    @inbounds for j in 2:n
        for i in 1:j-1
            e = first_energy[i,j]
            if e < Inf
                i_0, j_0 = nodos_ordenados[i] < nodos_ordenados[j] ? (nodos_ordenados[i], nodos_ordenados[j]) : (nodos_ordenados[j], nodos_ordenados[i])
                if e < first_energy_0[i_0, j_0]
                    first_energy_0[i_0, j_0] = e
                end
            end
        end
    end

    orbitas = compute_orbits_w(first_energy_0, n, epsilon_bins)

    println("\nTotal complete permutations: $perm_complete")
    println("Total prunned permutations: $perm_prunned")
    println("Elapsed time: $(round(time_taken, digits=2)) seconds")

    return epsilon_bins, aut_eps, orbitas
end
