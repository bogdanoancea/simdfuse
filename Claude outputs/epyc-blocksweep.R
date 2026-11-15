# Does the shallow-expression loss on EPYC come from the block size?
# sf_eval() takes its expression by substitution, so the quoted call has to be
# bound to a plain variable first -- sf_eval(exs[[nm]]) tries to compile the
# subscript expression itself.
library(simdfuse)
cat("isa:", sf_info()$isa, " lanes:", sf_info()$lanes,
    " threads:", sf_info()$max_threads,
    " block:", getOption("simdfuse.block"),
    " block_red:", getOption("simdfuse.block_red"), "\n\n")

set.seed(1); n <- 1e8
a <- runif(n, .5, 2); b <- runif(n, .5, 2)
tm <- function(f, r = 5) { f()
  min(replicate(r, { t <- as.numeric(Sys.time()); f(); as.numeric(Sys.time()) - t })) }

exs <- list("1 op"  = quote(a + b),
            "2 ops" = quote(a + b * 0.5),
            "4 ops" = quote((a + b * 0.5) * (a - b)),
            "8 ops" = quote((a + b * 0.5 - a * 0.25) * (a - b) + a / b - 1))
red <- quote((a - 1.25)^2)

base <- sapply(exs, function(e) tm(function() eval(e)))
base_red <- tm(function() sum((a - 1.25)^2))

cat("elementwise, single thread, speedup over base R\n")
cat(sprintf("%8s", "block")); for (nm in names(exs)) cat(sprintf("%10s", nm))
cat(sprintf("%12s\n", "sum reduce"))
for (blk in c(128L, 512L, 1024L, 4096L, 16384L, 65536L)) {
  sf_options(block = blk, block_red = blk)
  cat(sprintf("%8d", blk))
  for (nm in names(exs)) {
    e <- exs[[nm]]                                   # bind, then pass
    cat(sprintf("%9.2fx", base[[nm]] / tm(function() sf_eval(e, threads = 1))))
  }
  cat(sprintf("%11.2fx\n", base_red / tm(function() sf_sum(red, threads = 1))))
}

cat("\nand what the tuner picks on this machine:\n")
r <- sf_tune(n = 1e8)
