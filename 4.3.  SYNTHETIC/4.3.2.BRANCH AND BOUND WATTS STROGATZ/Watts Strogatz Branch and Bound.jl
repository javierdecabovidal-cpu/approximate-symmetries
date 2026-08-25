using Graphs
using Random

include("Prunning.jl")

p = 0.2
samples = 100
c = 5
n = 22
k = 6
for i in 1:samples
    g_0 = watts_strogatz(n, k, 0.0)
    g = watts_strogatz(n, k, 0.0)
    j = 0
    while Graphs.Experimental.has_isomorph(g_0, g)
        g = watts_strogatz(n, k, p, seed=samples*j+i)
        j += 1
    end
    edges_list = [(src(e), dst(e)) for e in edges(g)]
    CSV.write("ws_try/edges_02_ws$(i).csv",
          DataFrame(source = first.(edges_list),
                    target = last.(edges_list)))

    A = adjacency_matrix(g)
    epsilon_max = 4*c/n/(n-1)+ 1e-9
    eb, aut, orb = epsilon_aut(A, epsilon_max)
    mask = eb .<= epsilon_max
    CSV.write("ws_try/aut_02_$(i).csv",
        DataFrame(
            epsilon = eb[mask],
            aut = aut[mask]
        )
    )
end