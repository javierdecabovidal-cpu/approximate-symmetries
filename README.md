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

This code contains the function that computes the curve for simple, directed and weighted
networks, as well as the function that computes the orbits.

The function is called as:

```julia
Brute_Force(A, Directed, Weighted)
```

where `A` is the adjacency matrix of the network, and `Directed` and `Weighted` are booleans.

It returns:

- **Unweighted case:** the vector of epsilons, the histogram of energies, the CDF of the
  energies (i.e. the ε-automorphisms) and the `first_bin` matrix, which is used to compute
  the orbits.
- **Weighted case:** the vector of epsilons with at least one permutation, the CDF of the
  energies (the ε-automorphisms) and the `first_bin` matrix.

The orbits are obtained by calling `compute_orbits(first_bin, epsilon_bins)`, or
`compute_orbits_weighted(first_bin, epsilons)` if the network is weighted. Both return the
orbits for each epsilon.

### 3.2. Monte Carlo

The Monte Carlo function returns the mean value and the standard deviation of the energy.
It takes the adjacency matrix of a simple graph and the number of samples:

```julia
MonteCarlo_epsilon_mean_std(A, num_samples)
```

It returns the mean energy, the standard deviation, and the Monte Carlo error of both.

### 3.3. Branch and bound

The branch-and-bound approach explores the tree of partial permutations and discards a
branch as soon as its accumulated cost is guaranteed to exceed the target ε. It therefore
returns the same ε-automorphism curve as the brute-force method, but only up to a maximum
value `epsilon_max` that the user fixes in advance. Two bounds are combined: an incremental
adjacency cost (computed against the nodes already assigned) and a lower bound on the cost
of the nodes still to be assigned, obtained from the sorted degree (or strength) sequence.

Three implementations are provided, one per type of network:

| File | Network | Entry point |
|---|---|---|
| `Branch_and_Bound.jl` | simple (undirected, unweighted) | `epsilon_aut` |
| `Branch_and_Bound_directed.jl` | directed | `epsilon_aut` |
| `Branch_and_Bound_weighted.jl` | weighted (undirected) | `epsilon_aut_ponderado` |

All three are parallelised over the branches with `Threads.@threads`, so Julia must be
started with several threads, otherwise the code runs on a single core:

```bash
julia -t auto            # or: julia -t 8
```

The dependencies are `Combinatorics`, `Graphs`, `Random`, `Statistics`, `DelimitedFiles`,
`CSV`, `DataFrames` and, for the weighted version, `LinearAlgebra`.

The undirected and the directed files both define a function called `epsilon_aut`, so they
must not be included in the same Julia session: the second `include` silently overwrites the
first one.

#### Undirected and directed networks

```julia
include("Branch_and_Bound.jl")            # or Branch_and_Bound_directed.jl

epsilon_bins, aut_eps, orbits = epsilon_aut(A, epsilon_max)
```

where `A` is the adjacency matrix of the network (symmetric in the undirected case) and
`epsilon_max` is the largest ε to be explored, given as a fraction of the total number of
off-diagonal entries. The nodes are internally relabelled by decreasing degree (by
out-degree in the directed case), and all the results are mapped back to the original
labelling before being returned.

In the undirected case each mismatched pair contributes 2 to the cost, so the energy of a
permutation is always a multiple of 4 and the bins are spaced by `4/(n(n-1))`. In the
directed case the two directions of every pair are compared separately, each mismatch
contributes 1, and the bins are spaced by `2/(n(n-1))`.

#### Weighted networks

```julia
include("Branch_and_Bound_weighted.jl")

epsilon_bins, aut_eps, orbits = epsilon_aut_ponderado(W, epsilon_max)
```

Here `W` is the weight matrix, which must be square and symmetric; both conditions are
checked with an assertion. The energy of a permutation is the sum of the squared weight
differences, normalised by `n(n-1)`. Since the energies now form a continuous spectrum, the
nodes are ordered by strength and the degree bound is replaced by two bounds on the weights:
a cheap one derived from the Cauchy-Schwarz inequality, used at branch level, and a tighter
one based on the sorted rows (the rearrangement inequality), used to rank and cut the
candidates at each position. For a binary matrix the second bound reduces exactly to the
degree difference, so the unweighted case is recovered as a particular case.

An optional keyword argument caps the number of complete permutations stored per thread and
raises an error when it is exceeded, which is useful to abort a run whose `epsilon_max` turns
out to be too large:

```julia
epsilon_bins, aut_eps, orbits = epsilon_aut_ponderado(W, epsilon_max; max_complete = 10^7)
```

#### Output

The three functions return the same three objects:

- `epsilon_bins`: the values of ε. For the unweighted versions this is a regular grid
  covering the whole range `[0, 1]`; for the weighted version it contains only the distinct
  energies that were actually found, in increasing order.
- `aut_eps`: the ε-automorphism curve, i.e. the cumulative number of permutations whose
  energy is at most the corresponding ε.
- `orbits`: the orbits at each ε, obtained as the connected components of the graph whose
  edges join the nodes that have already become interchangeable. `orbits[k]` is the list of
  orbits at `epsilon_bins[k]`, and each orbit is a list of nodes in the original labelling.

Unlike the brute-force version, the orbits are computed inside the function, so
`compute_orbits` does not need to be called separately.

For the unweighted versions, only the entries with `epsilon_bins[k] <= epsilon_max` are
meaningful: beyond that point every branch has been pruned and the curve stays flat
artificially. The `orbits` vector is already truncated at that point, so a convenient way to
restrict the other two is

```julia
k_max = length(orbits)
epsilon_bins, aut_eps = epsilon_bins[1:k_max], aut_eps[1:k_max]
```

Each run also prints the number of branches pruned before the search starts, the number of
tasks and threads, the total number of complete and pruned permutations, and the elapsed
time.



