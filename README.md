# Approximate Symmetries in Networks

Repository for the Master's thesis "Approximate Symmetries in Networks: Theory, Validation and Applications" (2026), Javier de Cabo Vidal and Lucas Lacasa, Institute for Cross-Disciplinary Physics and Complex Systems (IFISC), University of the Balearic Islands.

This repository provides the data-generation pipeline used in the thesis. It includes the three approaches described there — brute force enumeration, Monte Carlo sampling of the permutation space, and a branch and bound algorithm — together with the scripts that generate all the raw data files behind the results reported in the thesis. Post-processing steps (e.g. hypothesis testing against null models, statistical analysis, and figure generation) are not included in this repository. Other networks, network models and tolerance ranges reported in the thesis follow the same structure and analysis procedure, differing only in the input graph and in the maximum number of allowed mismatches. The empirical networks analyzed in the thesis are not redistributed here; they are publicly available from Netzschleuder, the Network Repository and EcoBase. For any questions regarding the implementation or methodology, please feel free to contact me [javier.de-cabo2@estudiant.uib.cat](mailto:javier.de-cabo2@estudiant.uib.cat)

## Structure

Each folder corresponds to a section of the thesis.

| Path | Section | Contents |
|---|---|---|
| `3. ALGORITHMS/` | Sec. 3 | The three counting algorithms |
| ├ `Brute Force.jl` | 3.1 | Exhaustive enumeration of the $n!$ permutations |
| ├ `Monte Carlo.jl` | 3.2 | Uniform sampling of the permutation space |
| └ `BRANCH AND BOUND/` | 3.3 | Branch and bound, in its simple, directed and weighted variants |
| `4.2. CONCENTRATION OF THE MEASURE/` | Sec. 4.2 | Sampled vs. theoretical mean and variance of the energy |
| `4.3. SYNTHETIC/` | Sec. 4.3 | Perturbed rings (4.3.1) and Watts-Strogatz networks (4.3.2) |
| `4.4 EMPIRICAL/` | Sec. 4.4 | ε-automorphism curves for simple and directed empirical networks |
| `4.5. CLUSTER SYNCHRONIZATION/` | Sec. 4.5 | Coupled-map simulations on the ε-orbit partitions |

## 3. ALGORITHMS TO COUNT APPROXIMATE SYMMETRIES

In the folder `3. ALGORITHMS` you will find the three approaches developed to compute the ε-automorphism curve.
### 3.1. Brute force

This code contains the function that computes the curve for simple, directed and weighted networks, as well as the function that computes the orbits.

The function is called as:

```julia
Brute_Force(A, Directed, Weighted)
```

where `A` is the adjacency matrix of the network, and `Directed` and `Weighted` are booleans.

It returns:

- **Unweighted case:** the vector of epsilons, the histogram of energies, the CDF of the energies (i.e. the ε-automorphisms) and the `first_bin` matrix, which is used to compute the orbits.
- **Weighted case:** the vector of epsilons with at least one permutation, the CDF of the energies (the ε-automorphisms) and the `first_bin` matrix.

The orbits are obtained by calling `compute_orbits(first_bin, epsilon_bins)`, or `compute_orbits_weighted(first_bin, epsilons)` if the network is weighted. Both return the orbits for each epsilon.

### 3.2. Monte Carlo

The Monte Carlo function returns the mean value and the standard deviation of the energy. It takes the adjacency matrix of a simple graph and the number of samples:

```julia
MonteCarlo_epsilon_mean_std(A, num_samples)
```

It returns the mean energy, the standard deviation, and the Monte Carlo error of both.

### 3.3. Branch and bound

Instead of going through all the permutations, this method builds them node by node and drops a branch as soon as its cost is already above the epsilon we are looking for. It gives the same curve as the brute force, but only up to a maximum epsilon that we choose in advance, `epsilon_max`. The lower the `epsilon_max`, the faster it runs.

There is one file per type of network:

| File | Network | Function |
|---|---|---|
| `Branch_and_Bound.jl` | simple | `epsilon_aut` |
| `Branch_and_Bound_directed.jl` | directed | `epsilon_aut` |
| `Branch_and_Bound_weighted.jl` | weighted | `epsilon_aut_ponderado` |

The search is split between threads, so Julia has to be started with more than one:

```bash
julia -t auto            # or: julia -t 8
```

The packages used are `Combinatorics`, `Graphs`, `Random`, `Statistics`, `DelimitedFiles`, `CSV`, `DataFrames`, and `LinearAlgebra` for the weighted case.

The simple and the directed files both define a function called `epsilon_aut`, so they should not be loaded in the same session: the second one overwrites the first.

For simple and directed networks the function is called as:

```julia
epsilon_aut(A, epsilon_max)
```

where `A` is the adjacency matrix and `epsilon_max` is the largest epsilon we want to reach. In the simple case each mismatch counts 2, so the bins are spaced by `4/(n(n-1))`; in the directed case the two directions of each pair are checked separately, each mismatch counts 1, and the bins are spaced by `2/(n(n-1))`.

For weighted networks:

```julia
epsilon_aut_ponderado(W, epsilon_max)
```

where `W` is the weight matrix, which has to be square and symmetric. Here the energy is the sum of the squared differences of the weights, divided by `n(n-1)`, so it no longer falls on a fixed grid. The nodes are sorted by strength, and the degree bound of the unweighted case is replaced by a bound on the weights. There is also an optional argument, `epsilon_aut_ponderado(W, epsilon_max; max_complete = 10^7)`, that stops the run with an error if a thread finds more complete permutations than that, which is useful when `epsilon_max` turns out to be too large.

The three functions return the same three things:

- `epsilon_bins`: the values of epsilon. For simple and directed networks it is the full grid; for weighted ones it only contains the energies that actually appear, in increasing order.
- `aut_eps`: the curve, i.e. the number of permutations with energy up to each epsilon.
- `orbits`: the orbits at each epsilon, computed as the connected components of the graph that joins the nodes that are already interchangeable. `orbits[k]` are the orbits at `epsilon_bins[k]`, with the nodes labelled as in the original matrix.

Here the orbits are computed inside the function, so `compute_orbits` does not have to be called separately.

For simple and directed networks, only the values with `epsilon_bins[k] <= epsilon_max` mean anything: past that point everything has been pruned and the curve stays flat. `orbits` is already cut there, so the other two can be cut as well:

```julia
k_max = length(orbits)
epsilon_bins, aut_eps = epsilon_bins[1:k_max], aut_eps[1:k_max]
```

Each run prints how many branches were pruned before the search starts, how many tasks and threads are used, how many permutations were completed and pruned, and the time it took.
