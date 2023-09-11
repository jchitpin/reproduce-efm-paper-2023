# Uniquely identifying EFM weights by cycle-history Markov chain

This repository contains the scripts necessary to regenerate all figures and
results from the manuscript by J. G. Chitpin and T. J. Perkins. Certain values,
such as mean reconstruction error for the Markovian weights, are not exported
as a text file but can be viewed by interactively running the scripts. The
beginning of each script will list what is computed/exported.

## Requirements

The following is required:

* Julia (minimum version 1.6).
* Gurobi (must be version 9.12; free academic license available)
* A shell to run the scripts/workflows.
* TeX distribution (like TeX Live or MiKTeX) to compile figures.


## Notes on reproducibility

* Depending on your version of package dependency `DifferentialEquations.jl`,
  the numerically calculated flux vectors for Figure 4 may vary slightly
  (by roughly 10^-16). Thus, the EFM weights may vary slightly from
  MarkovWeightedEFMs.jl and the other optimization-based approaches. This code
  was last run with `DifferentialEquations v7.4.0`.
* All experiments conducted on a Ryzen 5950X with 16 cores allocated to Julia
  and MATLAB with 32 GB of memory (16 GB recommended to run all scripts).

## Installation

Download the repository and install all necessary Julia packages by running the
following commands in your desired installation directory.

1. `$ cd /home/<username>/<directory>/`  
2. `$ git clone jchitpin/reproduce-efm-paper-2023`  
3. `$ cd reproduce-efm-paper-2023/src/`  
4. `$ julia install-julia-packages.jl # or run line by line in Julia REPL`


## Workflow to reproduce results

The scripts in the following subsections should be run in order to regenerate
the intermediate data files. All scripts should be run in their current working
directory (`reproduce-efm-paper-2023/src/`).

### Figure 1 and 2

1. `$ julia main-efm-weights-example-markov.jl`
2. `$ julia main-efm-weights-example-optimization.jl`

Figures are generated via:

1. `$ sh figure-01.sh`
2. `$ sh figure-02.sh`

### Figure 3

Figure is generated via:

1. `sh figure-03.sh`

### Figure 4

1. `$ julia main-sphingo-network-validation.jl`
2. `$ julia main-efm-weights-sphingo-markov.jl`
3. `$ julia main-efm-weights-sphingo-optimization.jl`

Figures are generated via:

1. `sh figure-04.sh`

### Supplementary Figure 1

1. `$ julia main-sphingo-network-validation.jl`
2. `$ julia main-efm-weights-sphingo-markov.jl`
3. `$ julia main-efm-weights-sphingo-optimization.jl`

Figures are generated via:

1. `sh supplementary-01.sh`

### Supplementary Figure 2

1. `$ julia main-sphingo-network-validation.jl`
2. `$ julia main-efm-weights-sphingo-markov.jl`
3. `$ julia main-efm-weights-sphingo-optimization.jl`

Figures are generated via:

1. `sh supplementary-02.sh`

### Supplementary Figure 3/Table 1

1. `$ julia main-sphingo-network-validation.jl`
2. `$ julia main-efm-weights-sphingo-markov.jl`

Table is generated via:

1. `sh supplementary-03.sh`

## Reference

Justin G. Chitpin and Theodore J. Perkins, *A Markov constraint to uniquely identify elementary flux mode weights in unimolecular metabolic networks*, biorXiv preprint **biorXiv:2022.07.25.501464**, doi: https://doi.org/10.1101/2022.07.25.501464.


## Acknowledgements

We acknowledge the support of the Natural Sciences and Engineering
Research Council of Canada (NSERC), Discovery grant RGPIN-2019-0660 to
T.J.P. J.G.C. was supported by an NSERC CREATE Matrix Metabolomics Scholarship
and an NSERC Alexander Graham Bell Canada Graduate Scholarship.




