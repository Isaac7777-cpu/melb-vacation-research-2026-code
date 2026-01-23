# Repository Layout

```
src/            # my R simulation code
data/           # generated data, not committed due to size.
lib/            # supervisor’s code (leave unchanged up to formatting)
lib_data/       # supervisor's data (leave unchanged)
results/        # derived outputs (figures, tables, logs)
README.md       # this summary
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
