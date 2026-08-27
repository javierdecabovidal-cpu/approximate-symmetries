# Approximate Symmetries in Networks

Repository for the Master's thesis "Approximate Symmetries in Networks: Theory, Validation and Applications" (2026), Javier de Cabo Vidal and Lucas Lacasa, Institute for Cross-Disciplinary Physics and Complex Systems (IFISC), University of the Balearic Islands.

This repository provides the data-generation pipeline used in the thesis. It includes the three approaches described there — brute force enumeration, Monte Carlo sampling of the permutation space, and a branch and bound algorithm — together with the scripts that generate all the raw data files behind the results reported in the thesis. Post-processing steps (e.g. hypothesis testing against null models, statistical analysis, and figure generation) are not included in this repository. Other networks, network models and tolerance ranges reported in the thesis follow the same structure and analysis procedure, differing only in the input graph and in the maximum number of allowed mismatches. The empirical networks analyzed in the thesis are not redistributed here; they are publicly available from Netzschleuder, the Network Repository and EcoBase. For any questions regarding the implementation or methodology, please feel free to contact me [javier.de-cabo2@estudiant.uib.cat](mailto:javier.de-cabo2@estudiant.uib.cat)

## Structure

Each folder corresponds to a section of the thesis.

| Path | Section | Contents |
|---|---|---|
| `3. ALGORITHMS/` | Sec. 3 | The three counting algorithms |
| ├ `Brute Force.jl` | 3.1 | Exhaustive enumeration of the $n!$ permutations |
| ├ `Monte Carlo.jl` | 3.1 | Uniform sampling of the permutation space |
| └ `BRANCH AND BOUND/` | 3.1 | Branch and bound, in its simple, directed and weighted variants |
| `4.2. CONCENTRATION OF THE MEASURE/` | Sec. 4.2 | Sampled vs. theoretical mean and variance of the energy |
| `4.3. SYNTHETIC/` | Sec. 4.3 | Perturbed rings (4.3.1) and Watts-Strogatz networks (4.3.2) |
| `4.4 EMPIRICAL/` | Sec. 4.4 | ε-automorphism curves for simple and directed empirical networks |
| `4.5. CLUSTER SYNCHRONIZATION/` | Sec. 4.5 | Coupled-map simulations on the ε-orbit partitions |

## ALGORITHMS TO COUNT APPROXIMATE SYMMETRIES

In the folder 3. ALGORITHMS, you will find the three approaches developed to compute the epsilon-automorphism curve: 

# 3.1. Brute Force

It includes the function to compute the curve in simple, directed and weighted networks, as well as the function to compute the orbits. T
To call the function you must include: 

Brute_Force(A::(Adjacency matrix of the network), Directed(true or false), Undirected(true or false))

It returns:
- For the unweighted case: the vector of epsilons, the histogram energies, the CDF of the energies (the epsilon-automorphisms) and the first_bin matrix (used to compute the orbits).
- For the weighted case: the vector of epsilons with at least one permutation, the CDF of the energies (the epsilon-automorphisms) and the first_bin matrix (used to compute the orbits).

To compute the orbits we call to the function compute_orbits(first_bin, epsilon_bins) (if the network is weighted, we call compute_orbits_weighted(first_bin, epsilons)). It returns the orbits for each epsilon. 


# 3.1 Monte Carlo


The Monte Carlo function returns the mean value and standard deviation of the energy. We must entry with an adjacency matrix of a simple graph by introducing the matrix and the number of samples:

MonteCarlo_epsilon_mean_std(A(Adjacency matrix), num_samples(number of samples))

It returns the mean energy, the standard deviation and the Monte Carlo Error of both of them.



