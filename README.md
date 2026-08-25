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
