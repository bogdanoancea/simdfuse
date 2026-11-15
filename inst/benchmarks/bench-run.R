#!/usr/bin/env Rscript
# Full-profile benchmark run, for the author's reference figures and for anyone
# who wants the large-n numbers.
#
# Note that vignette("performance") does NOT read the CSV this writes. The
# vignette sources benchmarks.R and runs the same functions live on whatever
# machine builds it, so a reader sees their own hardware rather than ours.
#
# Usage:  Rscript inst/benchmarks/bench-run.R
#         SIMDFUSE_BENCH_SCALE=quick Rscript inst/benchmarks/bench-run.R

suppressMessages(library(simdfuse))

args <- commandArgs(trailingOnly = FALSE)
here <- grep("^--file=", args, value = TRUE)
here <- if (length(here)) dirname(sub("^--file=", "", here[1])) else "."
source(file.path(here, "benchmarks.R"))

p <- sf_bench_profile(Sys.getenv("SIMDFUSE_BENCH_SCALE", "full"))
cat("profile: ", p$scale, "\n", sep = "")
si <- sf_bench_sysinfo()
cat(sprintf("machine: %s, R %s, isa %s (%d lanes), %d threads\n",
            si$platform, si$r_version, si$isa, si$lanes, si$max_threads))
if (!is.na(si$cache)) cat("cache:   ", si$cache, "\n", sep = "")
cat("\n")

set.seed(20260114)
res <- sf_bench_all(p, verbose = TRUE)
out <- file.path(here, "results.csv")
write.csv(res, out, row.names = FALSE)
cat("\nwrote ", nrow(res), " rows to ", out, "\n", sep = "")
