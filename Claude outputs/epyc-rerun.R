# Re-run of the EPYC measurement with the tuning actually applied.
#
# The first run was taken at the defaults (block 128, block_red NA), which the
# run's own block sweep shows costs 2.18x on reductions and 1.43x on an
# eight-operator expression. Options live in the R session, so sf_tune() has to
# run in the SAME Rscript as the benchmark -- a separate invocation resets them.
library(simdfuse)
source(system.file("benchmarks", "benchmarks.R", package = "simdfuse"))

cat("== before tuning ==\n")
print(unlist(sf_options()[c("block","block_red","fast_exp","fast_log","isa")]))

# Tune at the size the work is actually done at: below the last-level cache the
# block size barely registers and the answer does not transfer upward. This
# machine has 256 MB of L3 per socket, so 1e7 (80 MB) is cache-resident and
# would tune for the wrong regime entirely.
cat("\n== tuning at n = 1e8 ==\n")
tn <- sf_tune(n = 1e8)

cat("\n== after tuning (these are what the benchmark will use) ==\n")
print(unlist(sf_options()[c("block","block_red","fast_exp","fast_log")]))

p <- sf_bench_profile("full")
p$ns       <- c(p$ns, 1e8)
p$n_thread <- 1.5e8

cat("\n== full sweep ==\n")
res <- sf_bench_all(p, verbose = TRUE)

# The reduce panel at 1e7 measures L3, not memory: 80 MB against 256 MB of L3.
# Repeat it past the cache so the figure means something for large data.
cat("\n== reductions out of cache, n = 5e8 ==\n")
p2 <- p; p2$n <- 5e8
res <- rbind(res, sf_bench_all(p2, which = "reduce", verbose = TRUE))

out <- path.expand("~/simdfuse-epyc-tuned.csv")
write.csv(res, out, row.names = FALSE)
cat("\nwrote", out, "with", nrow(res), "rows\n")

cat("\n== environment as measured ==\n")
print(sf_bench_sysinfo())
print(unlist(sf_info()))

# Anything this noisy is not a measurement. 22% of the first run's timings
# exceeded this, and two exceeded 100%.
if ("spread" %in% names(res)) {
  bad <- res[!is.na(res$spread) & res$spread > 0.15, ]
  cat(sprintf("\n%d of %d timings varied more than 15%% between repetitions\n",
              nrow(bad), sum(!is.na(res$spread))))
  if (nrow(bad)) print(bad[order(-bad$spread), c("bench","case","engine","seconds","spread")],
                       row.names = FALSE)
}
