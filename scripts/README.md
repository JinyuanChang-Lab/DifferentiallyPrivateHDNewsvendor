# R Scripts

This directory contains the R scripts used to reproduce the simulation studies,
real-data analysis, tables, and figures reported in the paper.

All scripts should be run from the root directory of the replication package.
For example:

```bash
Rscript scripts/simulations/main_sec2_iteration.R
```
 
 

## Directory Structure

```text
scripts/
├── figures/
│   ├── plot_grid.R
│   ├── plot_LDP.R
│   ├── plot_NN.R
│   └── plot_sec5_examples.R
├── real_data/
│   ├── main_sec6_real_data.R
│   ├── prepare_real_data.R
│   └── supp_NN_real_data.R
├── simulations/
│   ├── main_sec2_iteration.R
│   ├── main_sec5_example1.R
│   ├── main_sec5_example2.R
│   ├── main_sec5_example3.R
│   ├── main_sec5_example4.R
│   ├── supp_grid.R
│   ├── supp_LDP.R
│   ├── supp_NN_np.R
│   ├── supp_NN_priv.R
│   └── supp_runtime_comparison.R
└── README.md
```

## Software Environment

All computations are implemented in R. The required R packages are listed at
the beginning of each script.

### Windows Environment

The scripts run on the Windows computer were tested using R version 4.4.3 and
the following package versions:

| Package | Tested version |
|---|---:|
| `cowplot` | 1.1.3 |
| `doParallel` | 1.0.17 |
| `doRNG` | 1.8.6.2 |
| `doSNOW` | 1.0.20 |
| `dplyr` | 1.1.4 |
| `foreach` | 1.5.2 |
| `ggplot2` | 3.5.1 |
| `ggtext` | 0.1.2 |
| `here` | 1.0.2 |
| `lbfgs` | 1.2.1.2 |
| `lubridate` | 1.9.4 |
| `patchwork` | 1.3.0 |
| `purrr` | 1.0.4 |
| `qrnn` | 2.1.1 |
| `quantregForest` | 1.3.7.1 |
| `RSpectra` | 0.16.2 |
| `tidyr` | 1.3.1 |

The Windows scripts also use the `grid` package included with R.

### Linux Server Environment

The computationally intensive simulation scripts run on the Linux server were
tested using R version 4.5.2 and the following package versions:

| Package | Tested version |
|---|---:|
| `here` | 1.0.2 |
| `lbfgs` | 1.2.1.2 |
| `qrnn` | 2.1.1 |
| `quantreg` | 6.1 |
| `quantregForest` | 1.3.7.1 |

The Linux scripts also use the `parallel` package included with R.


## Computational Details

A complete mapping between the results reported in the paper, their output
locations, and the corresponding R scripts is provided in the root
`README.md`.

The following table reports the computational environment and approximate
elapsed runtime for the data-processing, simulation, and real-data analysis
scripts. Plotting scripts in `scripts/figures/` are omitted because their runtimes are negligible.

| Script                                          | Computational environment        |                         Approximate runtime |
| ----------------------------------------------- | -------------------------------- | ------------------------------------------: |
| **Results in the main paper**                                  |                                  |                                             |
| `scripts/simulations/main_sec2_iteration.R`     | Windows PC, Intel i7-14700F      |                        less than 10 seconds |
| `scripts/real_data/main_sec6_real_data.R`       | Windows PC, Intel i7-14700F      |                                  40 minutes |
| `scripts/real_data/prepare_real_data.R`         | Windows PC, Intel i7-14700F      |                                   5 minutes |
| `scripts/simulations/main_sec5_example1.R`      | Linux server using 167 CPU cores |                                  10 minutes |
| `scripts/simulations/main_sec5_example2.R`      | Linux server using 167 CPU cores |                                   2 minutes |
| `scripts/simulations/main_sec5_example3.R`      | Linux server using 167 CPU cores |                                   2 minutes |
| `scripts/simulations/main_sec5_example4.R`      | Linux server using 167 CPU cores |                                   5 minutes |
| **Results in the Electronic Companion**                   |                                  |                                             |
| `scripts/real_data/supp_NN_real_data.R`         | Windows PC, Intel i7-14700F      |                                  20 minutes |
| `scripts/simulations/supp_runtime_comparison.R` | Windows PC, Intel i7-14700F      |                                  10 minutes |
| `scripts/simulations/supp_grid.R`               | Linux server using 167 CPU cores |                                  10 minutes |
| `scripts/simulations/supp_LDP.R`                | Linux server using 167 CPU cores |                                  30 minutes |
| `scripts/simulations/supp_NN_np.R`              | Linux server using 167 CPU cores |                                   5 minutes |
| `scripts/simulations/supp_NN_priv.R`            | Linux server using 167 CPU cores |                                  26 hours |
 
### Important Note for Server Execution

When using a computing server, copy the entire replication package, except for
the `data/` directory, to the server rather than uploading individual scripts.
The `data/` directory is not required because the computationally intensive
scripts run on the server do not use the empirical data. This preserves the
remaining project directory structure and ensures that all relative paths can
be resolved correctly. Run these scripts from the root directory of the
replication package on the server.

After the computations are completed, copy the generated result files back to
the same relative locations in the local replication package. In particular,
files generated in `results/figure_data/` on the server should be copied to
`results/figure_data/` in the local package. The plotting scripts in
`scripts/figures/` can then be run locally to generate the final figures.