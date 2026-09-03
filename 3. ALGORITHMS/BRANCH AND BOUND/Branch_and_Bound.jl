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
# Al tener la matriz first_bin (primer bin al que i, j son intercambiables)
# podemos obtener las órbitas como las componentes conexas del grafo con edges 
# (parecido al proceso de percolación) entre los nodos que son intercambiables
# para el bin asociado al ϵ que queremos.
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
#
# Habiendo llegado a asignar k-1 nodos, podemos obtener una cota
# inferior de lo que costará asignar el resto de los nodos. 
# Ordenamos los nodos que quedan por asignar y los que todavía
# pueden ser asignados por grado. El coste mínimo de lo que queda por 
# asignar será la suma de diferencias de grado: 
#  \Sum_{i=k+1}^n |grado_nodo_i - grado_nodo_asignable_i|
# ---------------------------------------

@inline function cota_inferior_restante_grado(grado_ord::Vector{Int}, usado::BitVector, cand_global::Vector{Int},
    k::Int, n::Int)::Int

    cota_k = 0
    j_idx = 1

    # Recorremos los nodos que quedan por asignar 
    for i in k:n
        # Avanzamos hasta el siguiente nodo que pueda ser asignado.
        # Como cand_global está ordenado por grado, nos aseguramos
        # de que el nodo que asignemos tiene el mayor grado posible,
        # lo que minimiza la diferencia de grado.
        while j_idx <= n && usado[cand_global[j_idx]]
            j_idx += 1
        end
        j_idx > n && break  # Si no quedan nodos disponibles, salimos del bucle

        # Calculamos la diferencia de grado y la sumamos a la cota
        cota_k += abs(grado_ord[i] - grado_ord[cand_global[j_idx]])
        j_idx += 1
    end

    return cota_k
end

# ---------------------------------------
# PODA POR ADYACENCIA Y POR GRADOS
#
# En cada asignación π(k) = v aplicamos dos criterios de poda:
#
# (1) Coste incremental por aristas: contamos cuántas aristas entre
#     la columna k y las columnas ya asignadas 1...k-1 no coinciden
#     bajo π(k)=v. Cada discrepancia suma 2 al coste exacto acumulado.
#     Si coste_adyacencia_parcial + cota_A > cost_max, descartamos v.
#     Nótese que el coste de la permutación completa es exactamente
#     coste_adyacencia_parcia
# (2) Cota inferior por grados: calculamos el mínimo posible coste
#     adicional para completar las posiciones k...n usando solo
#     diferencias de grado. Si coste_grados_parcial + cota_inferior_restante_grado > cost_max,
#     descartamos toda la rama sin probar ningún candidato.
#
# Ambos criterios se aplican independientemente: basta con que uno
# de los dos supere cost_max para descartar la rama o el candidato.
# ---------------------------------------

function poda_recursiva!(Aw::Matrix{Int8}, grado_ord::Vector{Int}, cand_nodo::Vector{Vector{Int}}, 
    cost_max::Int, n::Int, pi_vec::Vector{Int}, usado::BitVector, k::Int, hist::Vector{Int},
    perm_complete::Base.RefValue{Int}, perm_prunned::Base.RefValue{Int}, first_bin::Matrix{Int32},
    coste_adyacencia_parcial::Int, coste_grados_parcial::Int, cand_global::Vector{Int})
     

     # Si el coste acumulado por grado más el coste mínimo de lo que queda por asignar es mayor 
     # que el coste máximo, descartamos toda la rama sin probar ningún candidato.
     if coste_grados_parcial + cota_inferior_restante_grado(grado_ord, usado, cand_global, k, n) > cost_max
        perm_prunned[] += 1
        return # Volvemos al anterior nodo
     end

     # Si ya hemos asignado todos los nodos, guardamos la permutación
     if k > n
        perm_complete[] += 1
        # Guardamos la permutación en su coste asociado
        # El coste de la permutación siempre es múltiplo de 4 
        # por la construcción de la métrica.
        bin_idx = (coste_adyacencia_parcial ÷ 4) + 1 
        hist[bin_idx] += 1
        # Guardamos el primer bin en el que pueden intercambiarse i, j
        @inbounds for i in 1:n
            j = pi_vec[i]
            if bin_idx < first_bin[i,j]
                first_bin[i,j] = bin_idx
            end
        end
        return # Volvemos al anterior nodo
     end

     # Recorremos los candidatos ordenados por grado
     @inbounds for v in cand_nodo[k]
        usado[v] && continue  # Si el nodo ya está asignado, lo saltamos
        # Calculamos el coste incremental por grado para asignar π(k) = v
        coste_k = abs(grado_ord[k] - grado_ord[v])
        nuevo_coste_grados_parcial = coste_grados_parcial + coste_k
        # Si el coste acumulado por grado más el coste mínimo de la permutación
        # que acabamos de asignar es mayor que el coste máximo, descartamos v sin probarlo.
        if nuevo_coste_grados_parcial > cost_max
            perm_prunned[] += 1

            break # No tiene sentido seguir probando candidatos con diferencia de grados aún mayores
        end

        # Calculamos el coste incremental por adyacencia al asignar π(k) = v
        coste_A = 0
        broke = false
        # Calculamos el coste por adyacencia al haber asignado π(k) = v, comparando la fila k con las filas ya asignadas 1...k-1
        # Si las entradas en la mtriz son distinas, representa un coste de +2, por simetría.
        @inbounds for i in 1:k-1
            if Aw[k, i] != Aw[v, pi_vec[i]]
                coste_A += 2
                if coste_adyacencia_parcial + coste_A > cost_max
                    broke = true
                    break # No tiene sentido seguir comparando si ya superamos el coste máximo
                end
            end
        end

        nuevo_coste_adyacencia_parcial = coste_adyacencia_parcial + coste_A
        if broke
            perm_prunned[] += 1
            continue # Probamos otro candidato
        end

        # Actualizamos el vector permutación y el vector de usados
        pi_vec[k] = v
        usado[v] = true

        # Llamada recursiva para asignar el siguiente nodo
        poda_recursiva!(Aw, grado_ord, cand_nodo, cost_max, n, pi_vec,
            usado, k+1, hist, perm_complete, perm_prunned, first_bin,
            nuevo_coste_adyacencia_parcial, nuevo_coste_grados_parcial, cand_global)

        # Si la rama fue descartada, ponemos v de nuevo disponible
        usado[v] = false

    end
    return 
end

"""
Para hacer la poda por grado y adyacencia, dado un coste máximo (asociado al número de missmatches
que queremos permitir), se llama a la función:

epsilon_aut(A::Matrix{Int8}, epsilon_max::Float64),

que devuelve la curva aut_eps y las permutaciones asociadas a cada coste.
"""

function epsilon_aut(A, epsilon_max)
    n = size(A, 1)

    # Calculamos el grado de cada nodo 
    grado = vec(sum(A, dims=2))

    # Ordenamos los nodos por grado y obtenemos el vector de grados ordenados
    nodos_ordenados = sortperm(grado, rev = true)
    # Permutamos la matriz para ordernarla por grado, así podemos acceder a las filas y columnas en el orden de grado.
    # Luego reconvertimos con el vector de nodos ordenados
    Aw = Matrix{Int8}(A[nodos_ordenados, nodos_ordenados])
    grado_ord = grado[nodos_ordenados]

    # Construimos el histograma de costes
    nbins = n * (n - 1) ÷ 4 + 1 # El coste máximo posible es n*(n-1) si todas las adyacencias son distintas, pero como cada discrepancia suma 2, el número de bins es la mitad más uno para incluir el caso de coste máximo.
    epsilon_bins = [4m / (n * (n - 1)) for m in 0:(n*(n-1)÷4)]
    cost_max::Int = floor(Int, epsilon_max * n * (n - 1) + 1e-10) # Redondeamos hacia abajo al múltiplo de 4 más cercano
    bin_max:: Int = (cost_max ÷ 4) + 1

    # Ordemaos los candidatos por grado descendente
    cand_global = sortperm(grado_ord, rev = true)

    # Para cada posición k, precomputamos el orden de los candidatos diferencia de grado ascendente.
    cand_nodo = Vector{Vector{Int}}(undef, n)
    for k in 1:n
        cand_nodo[k] = sortperm(abs.(grado_ord .- grado_ord[k])) 
    end

    # Inicializamos el prunning mirando todas las combinaciones de intercambio
    # entre los nodos de mayor grado, que son los que más contribuyen al coste por grado.
    # Esto nos da una cota inicial de coste máximo que nos permite podar muchas ramas desde el principio.
    pre_prunned = 0
    spawn_pairs = Tuple{Int, Int, Int, Int}[] #(v1, v2, coste_adyacencia_parcial, coste_grados_parcial)
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
            # Calculamos el coste por adyacencia de asignar π(1)=v1 y π(2)=v2
            coste_A = 0
            if Aw[1, 2] != Aw[v1, v2] # Si la adyacencia es distinta, el coste es 2
                coste_A += 2
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

    # Paralelizamos la búsqueda por branches
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

    # Sumamos los resultados de cada hilo
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

    # Mapeamos los índices de los nodos al original
    first_bin_0 = fill(typemax(Int32), n, n)
    @inbounds for j in 2:n
        for i in 1:j-1
            b = first_bin[i,j]
            if b < typemax(Int32)
                # Sólo nos quedamos con la parte triangular superior, que es la necesaria para construir las órbitas
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
