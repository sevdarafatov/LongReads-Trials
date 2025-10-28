library(data.table)

# ---------- 1) Load ----------
sample235 <- fread("DMR_results_2024A235.tsv")
sample250 <- fread("DMR_results_2024A250.tsv")

# ---------- 2) Significant calls ----------
p_cutoff <- 0.05
db_cutoff <- 0.20

sig_sample235 <- sample235[p < p_cutoff & abs(DeltaBeta) > db_cutoff,
  .(chr, start, end, DeltaBeta_sample235 = DeltaBeta, p_sample235 = p)]
sig_sample250 <- sample250[p < p_cutoff & abs(DeltaBeta) > db_cutoff,
  .(chr, start, end, DeltaBeta_sample250 = DeltaBeta, p_sample250 = p)]

# Ensure integer coordinates for foverlaps
sig_sample235[, `:=`(start = as.integer(start), end = as.integer(end))]
sig_sample250[, `:=`(start = as.integer(start), end = as.integer(end))]

setkey(sig_sample235, chr, start, end)
setkey(sig_sample250, chr, start, end)

# ---------- 3) Shared significant DMRs (>3 bp overlap) ----------
ov <- foverlaps(sig_sample235, sig_sample250, type = "any", nomatch = 0L)

# Compute overlap in bp
ov[, overlap_bp := pmin(end, i.end) - pmax(start, i.start) + 1L]
min_overlap_bp <- 3L
ov <- ov[overlap_bp > min_overlap_bp]

# Pick the best (largest overlap) per sample235 region
ov_best <- ov[, .SD[which.max(overlap_bp)], by = .(i.chr, i.start, i.end)]

# Tidy columns
shared <- ov_best[, .(
  chr                = i.chr,
  start              = i.start,
  end                = i.end,
  DeltaBeta_sample235 = i.DeltaBeta_sample235,
  p_sample235        = i.p_sample235,
  chr_sample250      = chr,
  start_sample250    = start,
  end_sample250      = end,
  DeltaBeta_sample250 = DeltaBeta_sample250,
  p_sample250        = p_sample250,
  overlap_bp
)]

# ---------- 4) Metrics ----------
direction_consistency <- mean(sign(shared$DeltaBeta_sample235) == sign(shared$DeltaBeta_sample250))
correlation <- cor(shared$DeltaBeta_sample235, shared$DeltaBeta_sample250)

# ---------- 5) Unique DMRs (sample235 only) ----------
sample235_flags <- foverlaps(sig_sample235, sig_sample250, type = "any", nomatch = NA)
sample235_flags[, overlap_bp := ifelse(
  is.na(start), 0L,
  pmin(end, i.end) - pmax(start, i.start) + 1L
)]
sample235_unique <- sample235_flags[overlap_bp <= min_overlap_bp,
  .(chr = i.chr, start = i.start, end = i.end,
    DeltaBeta_sample235 = i.DeltaBeta_sample235,
    p_sample235 = i.p_sample235)]
sample235_unique <- unique(sample235_unique)

# ---------- 6) Summary + Save ----------
cat("Significant DMRs (sample235):", nrow(sig_sample235), "\n")
cat("Significant DMRs (sample250):", nrow(sig_sample250), "\n")
cat("Shared significant DMRs (>3 bp):", nrow(shared), "\n")
cat("Unique significant DMRs (sample235 only):", nrow(sample235_unique), "\n")
cat("Direction consistency (shared):", round(100 * direction_consistency, 2), "%\n")
cat("Δβ correlation (shared):", round(correlation, 3), "\n")

fwrite(shared, "Shared_Significant_DMRs_min3bp_samples.tsv", sep = "\t")
fwrite(sample235_unique, "Unique_Significant_DMRs_sample235Only.tsv", sep = "\t")
