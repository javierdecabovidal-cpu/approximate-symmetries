using Graphs
using Random

include("Prunning_d.jl")   # version dirigida de epsilon_aut



function edge_swap_randomization(g2::DiGraph, nswaps; seed=801) # Función para hacer un shuffle de aristas que preserve los grados de salida y entrada.
    rng = MersenneTwister(seed)
    j = 0
    for _ in 1:nswaps
        edges_vec = collect(edges(g2))

        e1 = rand(rng, edges_vec)
        e2 = rand(rng, edges_vec)

        u, v = src(e1), dst(e1)
        x, y = src(e2), dst(e2)

        length(Set([u, v, x, y])) < 4 && continue

        has_edge(g2, u, y) && continue
        has_edge(g2, x, v) && continue

        rem_edge!(g2, u, v)
        rem_edge!(g2, x, y)

        add_edge!(g2, u, y)
        add_edge!(g2, x, v)

        j += 1
    end
    println("  Se hicieron $j swaps de aristas")

    return g2
end

# Verificación de que el swap preserva los grados de salida y entrada

function check_degree_preserved(g_orig::DiGraph, g_null::DiGraph)
    ok_out = outdegree(g_orig) == outdegree(g_null)
    ok_in  = indegree(g_orig) == indegree(g_null)
    if !ok_out || !ok_in
        println("  El swap NO preservo el grado correctamente ",
                "(outdeg ok=$ok_out, indeg ok=$ok_in)")
    end
    return ok_out && ok_in
end


df = CSV.read("validation/demp/edges.csv", DataFrame)
rename!(df, Symbol.(strip.(string.(names(df)))))
edges_df = select(df, [:source, :target])

# Si el CSV es 0-indexed descomenta:
#edges_df.source .+= 1
#edges_df.target .+= 1

n = maximum(vcat(edges_df.source, edges_df.target))
g = DiGraph(n)  

for row in eachrow(edges_df)
    row.source == row.target && continue   # Quitar autolazos
    add_edge!(g, row.source, row.target)
end

println("Red cargada: n=$(nv(g)), m=$(ne(g))")
println("outdegree: min=$(minimum(outdegree(g))), max=$(maximum(outdegree(g)))")
println("indegree:  min=$(minimum(indegree(g))), max=$(maximum(indegree(g)))")



num_samples = 60
c = 4
epsilon_max = 2*c/n/(n-1) + 1e-9

for i in 1:num_samples
    g_null = edge_swap_randomization(g, 150000; seed=i)

    check_degree_preserved(g, g_null)

    A = adjacency_matrix(g_null)   
    println("Sample $i: n=$(nv(g_null)), m=$(ne(g_null))")

    eb, aut, orb = epsilon_aut(A, epsilon_max)  
    mask = eb .<= epsilon_max
    CSV.write("validation/demp/aut_$(i).csv",
        DataFrame(
            epsilon = eb[mask],
            aut = aut[mask]
        ))
end
