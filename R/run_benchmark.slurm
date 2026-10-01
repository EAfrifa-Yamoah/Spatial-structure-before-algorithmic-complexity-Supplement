#!/bin/bash
#SBATCH --job-name=sdm_bench
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=sdm_bench_%j.out
#SBATCH --error=sdm_bench_%j.err
# Adjust partition, account and module names to the campus HPC.
# #SBATCH --partition=work
# #SBATCH --account=<project>

module load r/4.3   # or whatever provides R >= 4.3 with sf, TMB and compilers

# Package set (install once in an interactive session on the login node):
#   CRAN: dplyr tidyr sf ranger mgcv TMB RSQLite dbplyr DBI future.apply ggplot2
#         gstat Matrix terra lme4 units s2 RcppEigen Rcpp rnaturalearthdata
#         rnaturalearth blockCV sdmTMB fmesher lubridate checkmate rappdirs
#   GitHub: remotes::install_github("DFO-NOAA-Pacific/surveyjoin")
# Data: surveyjoin::cache_data(); surveyjoin::load_sql_data()  (login node, once;
#   needs outbound HTTPS to github.com; the cache lives in ~/.cache/surveyjoin)

export BENCH_RUN_FULL=TRUE
export BENCH_WORKERS=${SLURM_CPUS_PER_TASK}
export BENCH_REPS=10
export OMP_NUM_THREADS=1        # one thread per worker; parallelism is across replicates

cd "$SLURM_SUBMIT_DIR"
Rscript R/benchmark_scaffold.R
