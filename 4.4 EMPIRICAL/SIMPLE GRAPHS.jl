using Graphs
using Random

include("Prunning.jl")

function edge_swap_randomization(g, nswaps; seed=801) # Si hay problemas para generar los config model, se usa esta función como alternativa para hacer el shuffle degree-preserving.
    g2 = deepcopy(g)
    rng = MersenneTwister(seed)

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
    end

    return g2
end

df = CSV.read("validation/jazz/edges.csv", DataFrame)
rename!(df, Symbol.(strip.(string.(names(df)))))
edges_df = select(df, [:source, :target])


 edges_df.source .+= 1
 edges_df.target .+= 1

n = maximum(vcat(edges_df.source, edges_df.target))
g_orig = SimpleGraph(n)


for row in eachrow(edges_df)
    add_edge!(g_orig, row.source, row.target)
end

num_samples = 60
deg = CSV.read("validation/jazz/deg.csv", DataFrame, header=1)[:,2]
n = length(deg)
c = 6
epsilon_max = 4*c/n/(n-1)+ 1e-9

@threads for i in 1:num_samples
    #g_rand = edge_swap_randomization(g_orig, 150000; seed=i)
    g = random_configuration_model(n, deg, seed=i)
    A =adjacency_matrix(g)
    println("Sample $i: n=$(nv(g)), m=$(ne(g))")

    eb, aut, orb = epsilon_aut(A, epsilon_max)
    mask = eb .<= epsilon_max
    CSV.write("validation/jazz/aut_$(i).csv",
    DataFrame(
        epsilon = eb[mask],
        aut = aut[mask]
    ))
end

