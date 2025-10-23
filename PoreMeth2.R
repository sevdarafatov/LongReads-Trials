#############################################
#   PoreMeth2 Analysis Pipeline (Reusable)
#   Description: Runs full PoreMeth2 workflow
#   Input: modkit entropy TSV files (test & control)
#############################################

# ===== 1. Load required libraries =====
if (!require("data.table")) install.packages("data.table")
if (!require("devtools")) install.packages("devtools")

# Install PoreMeth2 if not yet installed
if (!require("PoreMeth2")) devtools::install_github("genepi-freiburg/PoreMeth2")
library(PoreMeth2)
library(data.table)

# ===== 2. Define file paths =====
# (Change these to your filenames for new samples)
file_control <- "2024A250_modkit.output_ds_sorted.entropy.file.tsv"  # Control sample (e.g. WT)
file_test    <- "2024A235_modkit.output_ds_sorted.entropy.file.tsv"  # Test sample (e.g. KI)

# ===== 3. Import modkit entropy files =====
TableControl <- fread(file_control)
TableTest    <- fread(file_test)

# ===== 4. Quality check for single experiment =====
# Control sample
png("QC_Control.png", width = 1800, height = 1200, res = 200)
PoreMeth2SingleExpQualityPlot(TableControl)
dev.off()

# Test sample
png("QC_Test.png", width = 1800, height = 1200, res = 200)
PoreMeth2SingleExpQualityPlot(TableTest)
dev.off()

# Paired QC plot (side-by-side)
png("QC_Paired_Control_vs_Test.png", width = 1800, height = 1000, res = 200)
PoreMeth2PairedExpQualityPlot(TableControl, TableTest)
dev.off()

# ===== 5. Run DMR detection =====
TableDMR <- PoreMeth2DMR(
  TableTest,
  TableControl,
  omega = 0.1,
  eta = 1e-5,
  FW = 3
)

# Save raw DMR results
fwrite(TableDMR, "DMR_downsampled_results.tsv", sep = "\t")

# ===== 6. Generate DMR summary statistics =====
stats_summary <- PoreMeth2DMRStatistics(
  TableDMR,
  Assembly = "hg38",
  BetaThr = 0.2,
  EntropyThr = 0.1,
  PValueThr = 0.05
)

# Convert to data frame & preserve row names
stats_summary_df <- as.data.frame(stats_summary)
stats_summary_df$Category <- rownames(stats_summary)
stats_summary_df <- stats_summary_df[, c("Category", setdiff(names(stats_summary_df), "Category"))]

# Save summary table
fwrite(stats_summary_df, "DMR_statistics_summary_ds.tsv", sep = "\t")

# ===== 7. Annotate DMR regions (using annotatr) =====
if (!require("annotatr")) install.packages("annotatr")
if (!require("GenomicRanges")) BiocManager::install("GenomicRanges")
if (!require("TxDb.Hsapiens.UCSC.hg38.knownGene")) {
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
}
if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
library(annotatr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# Load DMR results again to annotate
DMR_results <- fread("DMR_downsampled_results.tsv")

# Convert to GRanges
gr <- makeGRangesFromDataFrame(DMR_results,
                               seqnames.field = "chr",
                               start.field = "start",
                               end.field = "end")

# Build annotations (genes + CpG islands)
annotations <- build_annotations(genome = "hg38",
                                 annotations = c("hg38_basicgenes", "hg38_cpgs"))

# Annotate DMRs
annotated <- annotate_regions(
  regions = gr,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE
)

# Convert to data frame and save
AnnotatedTableDMR <- as.data.frame(annotated)
fwrite(AnnotatedTableDMR, "DMR_results_annotated_ds.tsv", sep = "\t")

# ===== 8. (Optional) Clean up environment =====
rm(TableControl, TableTest, TableDMR, stats_summary, annotated, gr)
gc()

cat("\n✅ PoreMeth2 analysis complete! All results saved in working directory.\n")
#############################################
