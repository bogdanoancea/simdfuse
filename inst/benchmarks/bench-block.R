#!/usr/bin/env Rscript
# Block-size sweep on its own, at the full profile. See bench-run.R for why the
# performance vignette does not read the CSV this writes.
#
# Usage:  Rscript inst/benchmarks/bench-block.R

suppressMessages(library(simdfuse))

args <- commandArgs(trailingOnly = FALSE)
here <- grep("^--file=", args, value = TRUE)
here <- if (length(here)) dirname(sub("^--file=", "", here[1])) else "."
source(file.path(here, "benchmarks.R"))

p <- sf_bench_profile(Sys.getenv("SIMDFUSE_BENCH_SCALE", "full"))
set.seed(20260114)
res <- sf_bench_block(p)
res$block <- as.integer(res$engine)
res$ms <- res$seconds * 1000

for (lab in unique(res$case)) {
  s <- res[res$case == lab, ]
  cat(sprintf("%-9s best block = %6d (%7.1f ms);  128 = %7.1f ms (%+.1f%%)\n",
              lab, s$block[which.min(s$ms)], min(s$ms),
              s$ms[s$block == 128L], 100 * (s$ms[s$block == 128L] / min(s$ms) - 1)))
}
out <- file.path(here, "block.csv")
write.csv(res[, c("case", "block", "ms")], out, row.names = FALSE)
cat("\nwrote ", nrow(res), " rows to ", out, "\n", sep = "")
