#!/usr/bin/env Rscript
# Standalone figures from a full-profile run, for the design note. The vignette
# draws its own figures from its own live measurements; this script exists to
# produce publication-sized PNGs from results.csv.
#
# Usage:  Rscript inst/benchmarks/bench-figures.R [results.csv] [outdir]

args <- commandArgs(trailingOnly = TRUE)
csv <- if (length(args) >= 1) args[1] else "results.csv"
out <- if (length(args) >= 2) args[2] else "fig"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
res <- read.csv(csv, stringsAsFactors = FALSE)

pal <- c("#1b6ca8", "#d1495b", "#66a182", "#edae49", "#2e4057")
open_png <- function(f, w = 1700, h = 1150) {
  png(file.path(out, f), width = w, height = h, res = 190)
  par(mar = c(4.2, 4.4, 2.6, 1.1), las = 1, bty = "l", mgp = c(2.6, 0.6, 0))
}
sec <- function(bench, case, engine) {
  s <- res[res$bench == bench & res$case == case & res$engine == engine, ]
  if (nrow(s)) s$seconds[1] else NA_real_
}

## 1. fusion
f <- subset(res, bench == "fusion")
ops <- c("1 op", "2 ops", "4 ops", "8 ops", "16 ops"); opn <- c(1, 2, 4, 8, 16)
ns <- sort(unique(f$n))
sp <- lapply(ns, function(nn) {
  s <- f[f$n == nn, ]
  vapply(ops, function(o) s$seconds[s$case == o & s$engine == "base R"] /
                          s$seconds[s$case == o & s$engine == "simdfuse"], 0)
})
open_png("fig1-fusion.png")
plot(NA, xlim = c(1, 16), ylim = range(0.5, unlist(sp) * 1.1, na.rm = TRUE),
     log = "x", xaxt = "n", xlab = "arithmetic operations in the expression",
     ylab = "speedup over base R  (x)",
     main = "Fusion: the win grows with expression depth")
axis(1, at = opn, labels = opn); abline(h = 1, col = "grey60", lty = 2)
for (i in seq_along(ns))
  lines(opn, sp[[i]], col = pal[i], lwd = 2.4, type = "b", pch = 19, cex = 0.9)
legend("topleft", sprintf("n = %.0e", ns), col = pal[seq_along(ns)], lwd = 2.4,
       pch = 19, bty = "n", cex = 0.9)
dev.off()

## 2. instruction set sweep
e <- subset(res, bench == "isa_exp")
ord <- c("scalar", "sse2", "avx2", "avx512", "neon")
e <- e[order(match(e$case, ord)), ]
open_png("fig2-isa.png", 1500, 1050)
bp <- barplot(e$seconds * 1e3, names.arg = e$case, col = pal[1], border = NA,
              ylab = sprintf("milliseconds (n = %.0e)", e$n[1]),
              main = "Runtime dispatch: vectorised exp() by instruction set",
              ylim = c(0, max(e$seconds * 1e3) * 1.2))
text(bp, e$seconds * 1e3, sprintf("%.0f ms", e$seconds * 1e3), pos = 3, cex = 0.85)
text(bp, e$seconds * 1e3 / 2, sprintf("%.1fx", e$seconds[1] / e$seconds),
     col = "white", font = 2, cex = 0.95)
dev.off()

## 3. reductions
r <- subset(res, bench == "reduce")
cs <- unique(r$case)
eng <- c("base R", "simdfuse 1 thread", "simdfuse all threads")
m <- sapply(cs, function(c0) vapply(eng, function(en) sec("reduce", c0, en), 0)) * 1e3
open_png("fig3-reductions.png", 1750, 1100)
bp <- barplot(m, beside = TRUE, col = pal[1:3], border = NA, names.arg = cs,
              ylab = sprintf("milliseconds (n = %.0e)", r$n[1]),
              main = "Fused reductions, no output vector allocated",
              ylim = c(0, max(m, na.rm = TRUE) * 1.22), cex.names = 0.78)
legend("topright", eng, fill = pal[1:3], border = NA, bty = "n", cex = 0.82)
for (j in seq_len(ncol(m)))
  text(bp[2:3, j], m[2:3, j], sprintf("%.1fx", m[1, j] / m[2:3, j]), pos = 3, cex = 0.7)
dev.off()

## 4. skeletons
vals <- c(sec("stencil", "5-point mean", "stats::filter") /
            sec("stencil", "5-point mean", "simdfuse"),
          sec("segreduce", "sum by group", "base rowsum") /
            sec("segreduce", "sum by group", "simdfuse"),
          sec("scan", "cumsum(log(x))", "base R") /
            sec("scan", "cumsum(log(x))", "simdfuse"))
alt <- c(sec("stencil", "5-point mean", "data.table::frollmean") /
           sec("stencil", "5-point mean", "simdfuse"),
         sec("segreduce", "sum by group", "collapse::fsum") /
           sec("segreduce", "sum by group", "simdfuse"), NA)
lab <- c("5-point stencil\n(vs stats::filter)", "sum by group\n(vs base rowsum)",
         "cumsum(log(x))\n(vs base)")
open_png("fig4-skeletons.png", 1600, 1080)
bp <- barplot(rbind(vals, alt), beside = TRUE, col = pal[c(1, 3)], border = NA,
              names.arg = lab, ylab = "speedup  (x)", cex.names = 0.8,
              main = "Parallel skeletons, fused with the expression",
              ylim = c(0, max(vals, alt, na.rm = TRUE) * 1.22))
abline(h = 1, col = "grey60", lty = 2)
legend("topright", c("vs base R", "vs fastest package alternative"),
       fill = pal[c(1, 3)], border = NA, bty = "n", cex = 0.82)
text(bp[1, ], vals, sprintf("%.1fx", vals), pos = 3, cex = 0.78)
ok <- !is.na(alt)
text(bp[2, ok], alt[ok], sprintf("%.1fx", alt[ok]), pos = 3, cex = 0.78)
dev.off()

## 5. block size
b <- subset(res, bench == "block")
if (nrow(b)) {
  b$block <- as.integer(b$engine)
  cs <- unique(b$case)
  open_png("fig5-block.png", 1600, 1080)
  rel <- lapply(cs, function(c0) { s <- b[b$case == c0, ]; s$seconds / min(s$seconds) })
  plot(NA, xlim = range(b$block), ylim = c(0.95, max(unlist(rel)) * 1.05), log = "x",
       xlab = "block size (doubles per register)",
       ylab = "time relative to the best for that expression",
       main = "Choosing the block size")
  for (i in seq_along(cs)) {
    s <- b[b$case == cs[i], ]
    lines(s$block, s$seconds / min(s$seconds), col = pal[i], lwd = 2.2,
          type = "b", pch = 19, cex = 0.8)
  }
  abline(v = 128, col = "grey45", lty = 3)
  legend("topleft", cs, col = pal[seq_along(cs)], lwd = 2.2, pch = 19,
         bty = "n", cex = 0.85)
  dev.off()
}
cat("figures written to ", out, ": ", paste(list.files(out), collapse = ", "), "\n", sep = "")
