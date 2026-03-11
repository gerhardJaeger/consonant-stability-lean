J=$(( $(nproc) / 4 ))
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 STAN_NUM_THREADS=2
find data/simulated/*.csv -type f -maxdepth 1 -print0 \
  | sort -z \
  | xargs -0 -n1 -P "$J" -I{} Rscript R/OU_stan_simulations.r "{}"
