# Repository Layout

```
./
├── data/ - generated datasets and fit results (read-only inputs to scripts).
│   ├── multi-res-10000draws.RData - 10k simulated 100×100 Gaussian field draws (by accident, took 10 mins).
│   ├── multi-res-500draws_no_nug.RData - 500 simulated field draws without nugget.
│   ├── multi-res-100draws_no_nug.RData - 100 simulated field draws without nugget.
│   ├── simulation_res.rds - MLE simulation fit results (alpha₀, φ, σ², etc.).
│   ├── simulation_res_reml.rds - REML fit results counterpart.
│   ├── simulat_res_fake_nll.rds - alternative fit results using a modified negative log-likelihood.
│   ├── x_simulation.rds - sampled coarse-grid observations/coords extracted from the simulated fields.
│   └── x_bundle.rds - packaged observation bundle used for parameter fitting in src/fit_x.R.
├── src/ - your active R code.
│   ├── simulate_grid.R - builds large Gaussian random fields on a 100×100 grid (with/without nugget) and saves draw stacks to data/.
│   ├── sample_from_grid.R - samples coarse-resolution boxes from the simulated fields, producing x_obs/coords bundles (x_simulation.rds).
│   ├── likelihood.R - negative log-likelihood definitions and MLE fitting (handles resolution aggregation and optional nugget).
│   ├── restricted_likelihood.R - REML-specific concentrated likelihood and fitting wrapper mirroring likelihood.R.
│   ├── fit_x.R - fits covariance parameters to a bundled observation set (x_bundle.rds) via MLE. (Individual version of `simulate_x.R`)
│   ├── simulate_x.R - sequential MLE/REML Monte Carlo loop over slices in x_simulation.rds; writes simulation_res*.rds.
│   ├── simulate_x_parallel.R - parallelized version of simulate_x.R using parallel package.
│   ├── analyse_simulation.R - post-analysis of simulation results (bias/variance/RMSE, quick plots/histograms).
│   ├── kriging.R - exponential covariance helper for kriging and small utilities.
│   └── run_kriging.R - perform kriging for any point in the grid.
├── lib/ - supervisor-supplied legacy helpers (avoid editing).
│   ├── BFuns.R - core MLE helper (profile nll, optimization wrapper, covariance construction).
│   ├── SimuData.R - Not used
│   ├── Test_Template.R - Not used
│   └── 3_CreateLocs.R - Not used
├── lib_data/ - supervisor’s reference data.
│   └── loc800.Rdata - Not used
├── README.md - this file.
├── .gitignore - ignores temporary outputs, R workspace files, etc.
└── .gitattributes - Used to enable Git LFS to store data file (under data/) in git.
```

## Why do I have parallel simulation?

Parallel processing of the simulation / optimisation procedure. However, the speeds up is not obvious and sometime can hurt performance for small number of simulations (<100). However, for the 500 number of simulations, I can see a clear speed-up for both REML and MLE.

```sh
# For MLE Simulation

Rscript src/simulate_x.R --num-core=4 --num-sim=500   1503.64s user 507.47s system 392% cpu 8:32.74 total
Rscript src/simulate_x_parallel.R --num-core=1 --num-sim=500 --data-file="data/x_simulation.rds"  1365.00s user 352.80s system 99% cpu 28:48.44 total

# For REML Simulation
Rscript src/simulate_x_parallel.R --use-reml=TRUE --num-core=4 --num-sim=500   1797.23s user 626.58s system 385% cpu 10:28.30 total
```
