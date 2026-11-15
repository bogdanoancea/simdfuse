#!/usr/bin/env bash
#
# simdfuse on the HPC cluster: install on the control node, measure on a
# compute node.
#
#   usage:  ./simdfuse-hpc.sh            install, smoke-test, submit the job
#           ./simdfuse-hpc.sh install    install and smoke-test only
#           ./simdfuse-hpc.sh submit     submit the measurement job only
#
# Run it from the control node. It deliberately does NOT benchmark there:
# login nodes are shared and throttled, so numbers taken on one are not worth
# having. The measurement goes to a compute node through sbatch.
#
set -euo pipefail

# ---------------------------------------------------------------- settings
REPO="${REPO:-$HOME/r-projects/simdfuse}"
RMODULE="${RMODULE:-r/4.6.2}"
OUTDIR="${OUTDIR:-$HOME/simdfuse-results}"
PARTITION="${PARTITION:-cpu-standard}"
NODE="${NODE:-}"                  # e.g. NODE=hpc-node10 to pin one node
WALLTIME="${WALLTIME:-04:00:00}"
N_THREAD_1SOCK="${N_THREAD_1SOCK:-1.5e8}"   # ~2290 tasks: ~35 per thread at 64
N_THREAD_2SOCK="${N_THREAD_2SOCK:-5e8}"     # ~7629 tasks: ~30 per thread at 256

say() { printf '\n\033[1;34m== %s\033[0m\n' "$*"; }
die() { printf '\n\033[1;31m!! %s\033[0m\n' "$*" >&2; exit 1; }

load_r() {
  # `module` is a shell function, so a non-interactive shell may not have it
  if ! command -v module >/dev/null 2>&1; then
    for f in /etc/profile.d/modules.sh /usr/share/Modules/init/bash \
             /etc/profile.d/lmod.sh /usr/share/lmod/lmod/init/bash; do
      [ -r "$f" ] && . "$f" && break
    done
  fi
  command -v module >/dev/null 2>&1 || die "no 'module' command; load R yourself and re-run"
  module load "$RMODULE"
  command -v Rscript >/dev/null 2>&1 || die "Rscript not on PATH after 'module load $RMODULE'"
}

# ---------------------------------------------------------------- install
do_install() {
  load_r
  say "R and library paths"
  R --version | head -1
  # The module's own library is read-only, so make sure the personal one exists.
  # Once the directory is there R picks it up automatically and no -l is needed.
  local userlib; userlib="$(Rscript -e 'cat(Sys.getenv("R_LIBS_USER"))')"
  mkdir -p "$userlib"
  echo "personal library: $userlib"
  Rscript -e 'cat(.libPaths(), sep="\n")'

  say "optional packages the benchmark compares against"
  # requireNamespace-guarded in the benchmark, so missing ones only drop rows
  Rscript -e 'for (p in c("collapse","data.table")) if (!requireNamespace(p, quietly=TRUE))
                install.packages(p, lib=Sys.getenv("R_LIBS_USER"),
                                 repos="https://cloud.r-project.org")
              cat("collapse:", requireNamespace("collapse", quietly=TRUE),
                  " data.table:", requireNamespace("data.table", quietly=TRUE), "\n")'

  say "source"
  [ -d "$REPO" ] || die "no repository at $REPO (git clone it first)"
  git -C "$REPO" pull --ff-only 2>/dev/null \
    || echo "note: could not fast-forward; using the working tree as it stands"
  git -C "$REPO" --no-pager log --oneline -1

  say "build and install"
  # Installing from the source directory rather than a tarball skips vignette
  # building, so rmarkdown is not needed. --preclean is NOT optional: the five
  # sf_kernels_<isa>.c files #include sf_kernels.h and make does not infer that
  # prerequisite, so stale objects would be relinked silently.
  R CMD INSTALL --preclean "$REPO"

  say "which kernel the dispatcher chose"
  Rscript -e 'library(simdfuse); str(sf_info())'
  echo "on EPYC 7713 (Zen 3) expect isa=avx2+fma and lanes=4; Zen 3 has no AVX-512"

  say "correctness, without testthat"
  Rscript -e '
    library(simdfuse)
    set.seed(1)
    x <- runif(1e5, .5, 3); y <- runif(1e5, .5, 3); mu <- mean(x); s <- sd(x)
    i <- sample(c(1:9L, NA_integer_), 1e5, TRUE)
    ok <- c(
      fused   = identical(sf_eval((x + y*0.5 - x*0.25)*(x - y) + x/y - 1),
                          (x + y*0.5 - x*0.25)*(x - y) + x/y - 1),
      math    = identical(sf_eval(exp(x) + log(y)), exp(x) + log(y)),
      power   = identical(sf_eval((x - mu)^2), (x - mu)^2),
      modulo  = identical(sf_eval(x %% y), x %% y),
      integer = identical(sf_eval(i * 2), as.double(i) * 2),
      na_kept = identical(is.na(sf_eval(i + 0)), is.na(as.double(i))),
      loglik  = isTRUE(all.equal(
                  sf_sum(~ -0.5*log(2*pi) - log(s) - 0.5*((x-mu)/s)^2),
                  sum(-0.5*log(2*pi) - log(s) - 0.5*((x-mu)/s)^2))),
      det_thr = identical(sf_sum(quote(x), threads = 1),
                          sf_sum(quote(x), threads = 0)))
    print(ok)
    if (!all(ok)) { cat("FAILURES ABOVE\n"); quit(status = 1) }
    cat("all correctness checks passed\n")'
  say "install OK"
}

# ---------------------------------------------------------------- measure
write_rscript() {
  mkdir -p "$OUTDIR"
  # quoted heredoc: R's $ must not be touched by the shell
  cat > "$OUTDIR/run-bench.R" <<'RSCRIPT'
library(simdfuse)
source(system.file("benchmarks", "benchmarks.R", package = "simdfuse"))

args  <- commandArgs(TRUE)
which <- args[1]                       # "1socket" or "2socket"
nthr  <- as.numeric(args[2])
out   <- args[3]

cat("\n== sf_info ==\n"); print(unlist(sf_info()))
cat("\n== sf_bench_sysinfo ==\n"); print(sf_bench_sysinfo())

p <- sf_bench_profile("full")
p$n_thread <- nthr

if (identical(which, "1socket")) {
  # One socket, data bound to the same NUMA node: the clean scaling curve,
  # with no remote-memory traffic in it.
  cat("\n== tuning (block sizes and per-function fast math) ==\n")
  tn <- sf_tune(n = 1e7)
  saveRDS(tn, file.path(dirname(out), "tune.rds"))
  p$ns <- c(p$ns, 1e8)
  res <- sf_bench_all(p, verbose = TRUE)
} else {
  # Both sockets. R allocates single-threaded, so first touch puts every page
  # on one node; threads on the far socket then read across the interconnect.
  # Comparing this curve with the single-socket one is what separates the
  # package's scaling from the machine's memory topology.
  res <- sf_bench_threads(p)
}
write.csv(res, out, row.names = FALSE)
cat("\nwrote", out, "with", nrow(res), "rows\n")

# Flag anything too noisy to interpret -- a lesson this project learned the
# hard way, after a 1.27x swing in a figure nothing had changed.
if ("spread" %in% names(res)) {
  bad <- res[!is.na(res$spread) & res$spread > 0.15, ]
  if (nrow(bad)) {
    cat("\n", nrow(bad), " of ", nrow(res),
        " timings varied by more than 15% between repetitions:\n", sep = "")
    print(bad[order(-bad$spread), c("bench","case","engine","seconds","spread")],
          row.names = FALSE)
  } else cat("\nall timings within 15% between repetitions\n")
}
RSCRIPT
}

do_submit() {
  write_rscript
  local nodeopt=""
  [ -n "$NODE" ] && nodeopt="#SBATCH --nodelist=$NODE"
  cat > "$OUTDIR/job.sbatch" <<SBATCH
#!/usr/bin/env bash
#SBATCH --job-name=simdfuse
#SBATCH --partition=$PARTITION
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --exclude=hpc-ctrl01,hpc-db01
#SBATCH --time=$WALLTIME
#SBATCH --output=$OUTDIR/slurm-%j.log
$nodeopt
set -euo pipefail
for f in /etc/profile.d/modules.sh /usr/share/Modules/init/bash \\
         /etc/profile.d/lmod.sh /usr/share/lmod/lmod/init/bash; do
  [ -r "\$f" ] && . "\$f" && break
done
module load $RMODULE

# Pin to physical cores. On this box CPUs 0-63 and 128-191 are node 0, where
# 128-191 are SMT siblings of 0-63, so without this threads land on sibling
# hyperthreads while physical cores sit idle.
export OMP_PROC_BIND=close OMP_PLACES=cores
# R may be linked against a threaded BLAS that would spawn its own pool
# during data setup and compete with the measurement.
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=\${SLURM_CPUS_ON_NODE:-256}

echo "node=\$(hostname)  cpus=\$(nproc)  job=\$SLURM_JOB_ID"
lscpu | grep -iE 'model name|socket|core\\(s\\) per|thread\\(s\\) per|^NUMA node' || true

cd "$REPO"

echo; echo "################ single socket, memory bound to the same node"
if command -v numactl >/dev/null 2>&1; then
  numactl --cpunodebind=0 --membind=0 \\
    Rscript "$OUTDIR/run-bench.R" 1socket $N_THREAD_1SOCK "$OUTDIR/epyc-1socket.csv"
else
  echo "numactl not found: running unbound, so this curve carries NUMA effects"
  Rscript "$OUTDIR/run-bench.R" 1socket $N_THREAD_1SOCK "$OUTDIR/epyc-1socket.csv"
fi

echo; echo "################ both sockets, thread sweep only"
Rscript "$OUTDIR/run-bench.R" 2socket $N_THREAD_2SOCK "$OUTDIR/epyc-2socket.csv"

echo; echo "done. results in $OUTDIR"
ls -la "$OUTDIR"
SBATCH

  say "submitting"
  local jid
  jid="$(sbatch --parsable "$OUTDIR/job.sbatch")"
  echo "job $jid submitted to $PARTITION${NODE:+ on $NODE}"
  cat <<EOF

monitor:
  squeue -j $jid
  tail -f $OUTDIR/slurm-$jid.log
cancel:
  scancel $jid

when it finishes, these are the files to send back:
  $OUTDIR/epyc-1socket.csv
  $OUTDIR/epyc-2socket.csv
  $OUTDIR/tune.rds
  $OUTDIR/slurm-$jid.log
EOF
}

case "${1:-all}" in
  install) do_install ;;
  submit)  do_submit ;;
  all)     do_install; do_submit ;;
  *)       die "usage: $0 [install|submit|all]" ;;
esac
