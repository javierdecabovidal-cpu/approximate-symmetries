# Approximate Symmetries in Networks

Repository for the Master's thesis "Approximate Symmetries in Networks: Theory, Validation and Applications" (2026), Javier de Cabo Vidal and Lucas Lacasa, Institute for Cross-Disciplinary Physics and Complex Systems (IFISC), University of the Balearic Islands.

This repository provides the data-generation pipeline used in the thesis. It includes the three approaches described there — brute force enumeration, Monte Carlo sampling of the permutation space, and a branch and bound algorithm — together with the scripts that generate all the raw data files behind the results reported in the thesis. Post-processing steps (e.g. hypothesis testing against null models, statistical analysis, and figure generation) are not included in this repository. Other networks, network models and tolerance ranges reported in the thesis follow the same structure and analysis procedure, differing only in the input graph and in the maximum number of allowed mismatches. The empirical networks analyzed in the thesis are not redistributed here; they are publicly available from Netzschleuder, the Network Repository and EcoBase. For any questions regarding the implementation or methodology, please feel free to contact me [javier.de-cabo2@estudiant.uib.cat](mailto:javier.de-cabo2@estudiant.uib.cat)

## Setup

```
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```
