#####################
# LOAD LIBRARIES    #
#####################
library(dplyr)
library(DESeq2)
library(pheatmap)
library(scales)
library(rtracklayer)

########
# INIT #
########
workingDir   <- "/media/gbim/bigdata/folder" # Le chemin vers dossier dans lequel on travaille 
RDir         <- file.path(workingDir, "R")
countDir     <- file.path(workingDir, "count_mRNA")
conditionFiles <- c("conditions_CMT2A_vs_WT.tsv") # Le fichier condition avec .tsv à la fin

# Annotation file (GTF)
args <- c(
  workingDir,
  "/home/gbim/syno/Bioinfo_scripts/genomic_DB/GRCm39/Mus_musculus.GRCm39.108.chr.gtf"
)

########################
# GENE LISTS OF INTEREST
########################

geneList1 <- c("ABCD1", "ACAA1", "ACBD4", "ACBD5", "ACOT4", "ACOT6", "ACOT8", "ACOX1", "ACOX2",
               "ACOX3", "ACOXL", "ALDH3A2", "AMACR", "CRAT", "CROT", "DECR2", "ECI2", "EHHADH",
               "HACL1", "HAO2", "HSD17B4", "MLYCD", "NUDT19", "NUDT7", "PECR", "PHYH", "SCP2",
               "SLC25A17", "SLC27A2")

geneList2 <- c("ACAA2", "ACAD10", "ACAD11", "ACADL", "ACADM", "ACADS", "ACADVL", "ACBD6", "ACBD7",
               "ACOT1", "ACOT11", "ACOT12", "ACOT13", "ACOT2", "ACOT7", "ACOT9", "ACSF2", "ACSM3",
               "ACSM6", "DBI", "DECR1", "ECHS1", "ECI1", "HADH", "HADHA", "HADHB", "MCAT", "MCEE",
               "MECR", "MMAA", "MMUT", "NDUFAB1", "PCCA", "PCCB", "PCTP", "THEM4", "THEM5")

geneList3 <- c("ABCC1", "ABCD1", "ACAA1", "ACAA2", "ACACA", "ACACB", "ACAD10", "ACAD11", "ACADL",
               "ACADM", "ACADS", "ACADVL", "ACBD4", "ACBD5", "ACBD6", "ACBD7", "ACLY", "ACOT1",
               "ACOT11", "ACOT12", "ACOT13", "ACOT2", "ACOT4", "ACOT6", "ACOT7", "ACOT8", "ACOT9",
               "ACOX1", "ACOX2", "ACOX3", "ACOXL", "ACSBG1", "ACSBG2", "ACSF2", "ACSF3", "ACSL1",
               "ACSL3", "ACSL4", "ACSL5", "ACSL6", "ACSM3", "ACSM6", "AKR1C3", "ALDH3A2", "ALOX12",
               "ALOX12B", "ALOX15", "ALOX15B", "ALOX5", "ALOX5AP", "ALOXE3", "AMACR", "AWAT1",
               "CBR1", "CBR4", "CPT1A", "CPT1B", "CPT2", "CRAT", "CROT", "CYP1A1", "CYP1A2",
               "CYP1B1", "CYP2C19", "CYP2C8", "CYP2C9", "CYP2J2", "CYP2U1", "CYP4A11", "CYP4A22",
               "CYP4B1", "CYP4F11", "CYP4F2", "CYP4F22", "CYP4F3", "CYP4F8", "CYP8B1", "DBI",
               "DECR1", "DECR2", "DPEP1", "DPEP2", "DPEP3", "ECHS1", "ECI1", "ECI2", "EHHADH",
               "ELOVL1", "ELOVL2", "ELOVL3", "ELOVL4", "ELOVL5", "ELOVL6", "ELOVL7", "EPHX2",
               "FAAH", "FAAH2", "FADS1", "FADS2", "FASN", "GGT1", "GGT5", "GPX1", "GPX2", "GPX4",
               "HACD1", "HACD2", "HACD3", "HACD4", "HACL1", "HADH", "HADHA", "HADHB", "HAO2",
               "HPGD", "HPGDS", "HSD17B12", "HSD17B3", "HSD17B4", "HSD17B8", "HTD2", "LTA4H",
               "LTC4S", "MAPKAPK2", "MCAT", "MCEE", "MECR", "MID1IP1", "MLYCD", "MMAA", "MMUT",
               "MORC2", "NDUFAB1", "NUDT19", "NUDT7", "OLAH", "PCCA", "PCCB", "PCTP", "PECR",
               "PHYH", "PLA2G4A", "PON1", "PON2", "PON3", "PPARD", "PPT1", "PPT2", "PRKAA2",
               "PRKAB2", "PRKAG2", "PRXL2B", "PTGDS", "PTGES", "PTGES2", "PTGES3", "PTGIS",
               "PTGR1", "PTGR2", "PTGS1", "PTGS2", "RXRA", "SCD", "SCD5", "SCP2", "SLC22A5",
               "SLC25A17", "SLC25A20", "SLC27A2", "SLC27A3", "TBXAS1", "TECR", "TECRL", "THEM4",
               "THEM5", "THRSP")

geneList4 <- c("ABCD1", "ACAA1", "ACOT4", "ACOT8", "ACOX1", "DECR2", "ECI2", "EHHADH",
               "HSD17B4", "MLYCD", "SLC27A2")

geneList5 <- c("ABCB1", "ABCB4", "ABCB9", "ABCC5", "ABCC8", "ABCD1", "ABCD2", "ABCD3", "ACAA1",
               "ACOT8", "ACOX1", "ACSL1", "ACSL4", "ACSL5", "ALB", "ALDH1A1", "ALDH9A1", "ATXN1",
               "BCL10", "CACNA1B", "CADM1", "CAT", "CDK7", "CEL", "CLN6", "CLN8", "CNBP", "CRABP1",
               "CRABP2", "CRAT", "CTBP1", "CTPS1", "DHCR24", "DHRS3", "DIO1", "DLG4", "ECH1",
               "ECI2", "EHHADH", "ELOVL5", "EPHX2", "ERCC1", "ERCC3", "ESR2", "FABP6", "FADS1",
               "FDPS", "FIS1", "GNPAT", "GSTK1", "HAO2", "HMGCL", "HRAS", "HSD11B2", "HSD17B11",
               "HSD17B4", "HSD3B7", "IDE", "IDH1", "IDH2", "IDI1", "ISOC1", "ITGB1BP1", "LONP2",
               "MFN1", "MFN2", "MLYCD", "MSH2", "MVP", "NR1I2", "NUDT19", "PABPC1", "PEX11A",
               "PEX11B", "PEX13", "PEX14", "PEX2", "PEX5", "PEX6", "PRDX1", "PRDX5", "RDH11",
               "RETSAT", "RXRG", "SCGB1A1", "SCP2", "SEMA3C", "SERPINA6", "SIAH1", "SLC23A2",
               "SLC25A17", "SLC25A19", "SLC25A4", "SLC27A2", "SLC35B2", "SMARCC1", "SOD1", "SOD2",
               "STS", "SULT2B1", "TOP2A", "TSPO", "TTR", "UGT2B17", "VPS4B", "YWHAH")

# Named list of gene sets (for TSVs)
geneList_sets <- list(
  set1 = geneList1,
  set2 = geneList2,
  set3 = geneList3,
  set4 = geneList4,
  set5 = geneList5
)

# For the heatmap, we will use geneList5:
geneList <- geneList5

###################
# ADAPT GENE LIST #
###################
# NE PAS TOUCHER C'EST POUR CONVERTIR LES GÈNES HUMAINS EN GÈNES SOURIS

mouse_human_genes <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",
                              sep = "\t")

# separate human and mouse 
mouse <- split.data.frame(mouse_human_genes, mouse_human_genes$Common.Organism.Name)[[2]]
human <- split.data.frame(mouse_human_genes, mouse_human_genes$Common.Organism.Name)[[1]]

# remove some columns
mouse <- mouse[, c(1, 4)]
human <- human[, c(1, 4)]

# merge the 2 dataset  (note that the human list is longer than the mouse one)
mh_data <- merge.data.frame(mouse, human, by = "DB.Class.Key", all.y = TRUE)
colnames(mh_data) <- c("DB.class.Key", "Symbol.Mouse", "Symbol.Human")

# Function to convert one gene list to mouse symbols
convert_gene_list_to_mouse <- function(geneList, mh_data) {
  newGeneList <- c()
  for (gene in geneList) {
    if (gene %in% mh_data$Symbol.Mouse) {
      # Already mouse symbol
      newGeneList <- c(newGeneList, gene)
    } else {
      # Map human symbol to mouse
      newGeneName <- mh_data %>%
        dplyr::filter(Symbol.Human == gene) %>%
        dplyr::pull(Symbol.Mouse)
      if (length(newGeneName) >= 1) {
        newGeneList <- c(newGeneList, newGeneName[1])
      } else {
        warning(paste("Gene not found in mapping:", gene))
      }
    }
  }
  unique(newGeneList)
}

# Mouse-converted lists for ALL 5 sets
mouse_geneList_sets <- lapply(geneList_sets, convert_gene_list_to_mouse, mh_data = mh_data)

# For the heatmap we use geneList5 converted to mouse
newGeneList <- convert_gene_list_to_mouse(geneList, mh_data)

###########
# ANALYSE #
###########
DATA_DIR <- args[1]
GTF_FILE <- args[2]

gtf <- rtracklayer::import(GTF_FILE)
gtf_df_transcripts <- as.data.frame(gtf) %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::select(transcript_id, gene_name)

gtf_df_genes <- as.data.frame(gtf) %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::select(gene_id, gene_name)
gtf_df_genes <- gtf_df_genes[!duplicated(gtf_df_genes), ]

print("gtf file loaded")

for (conditionFile in conditionFiles) {
  
  CONDITION_FILE   <- file.path(workingDir, conditionFile)
  nameConditionFile <- tools::file_path_sans_ext(basename(CONDITION_FILE))
  
  OUTPUT_DIR <- paste0(DATA_DIR, "/R/", nameConditionFile)
  OUTPUT_DIR <- gsub("//", "/", OUTPUT_DIR)
  
  # Make sure output dir exists
  if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, recursive = TRUE)
  }
  
  COUNT_MRNA <- paste0(DATA_DIR, "/count_mRNA/")
  COUNT_MRNA <- gsub("//", "/", COUNT_MRNA)
  
  ########
  # Init #
  ########
  metaDataName <- CONDITION_FILE
  metaData <- read.csv(metaDataName, header = TRUE, sep = "\t", row.names = 1)
  c <- cbind(metaData, Sample = rownames(metaData))
  
  ########
  # mRNA #
  ########
  countsName_mRNA <- paste0(COUNT_MRNA, "/gene_count_matrix.tsv")
  countData_mRNA <- read.csv(countsName_mRNA, header = TRUE, sep = ",", row.names = 1)
  colnames(countData_mRNA) <- gsub("X", "", colnames(countData_mRNA))
  
  # Check sample consistency
  if (all(rownames(metaData) %in% colnames(countData_mRNA)) == TRUE) {
    countData_mRNA <- countData_mRNA[, rownames(metaData)]
    if (all(colnames(countData_mRNA) %in% rownames(metaData)) == FALSE) {
      cat("Attention, pas tous les échantillons présents dans la matrice sont dans le fichier conditions, seuls les échantillons du fichier conditions ont été faits")
    }
  } else {
    stop("Samples in condition file and count matrix do not match.")
  }
  
  ##############################
  # DESeq2 ANALYSE DIFFERENTIEL
  ##############################
  
  dds_mRNA <- DESeqDataSetFromMatrix(countData = countData_mRNA,
                                     colData   = metaData,
                                     design    = ~ condition)
  # Relevel condition: reference = last condition in metadata
  dds_mRNA$condition <- relevel(dds_mRNA$condition,
                                ref = metaData$condition[length(metaData$condition)])
  
  dds_mRNA <- DESeq(dds_mRNA)
  vsd_mRNA <- vst(dds_mRNA, blind = FALSE) #Normalization method is vst 
  
  ############################
  # DESeq2 RESULTS (ALL)     #
  ############################
  res <- results(dds_mRNA)
  res_ordered <- res[order(res$padj), ]
  
  res_df <- as.data.frame(res_ordered)
  res_df$gene_raw <- rownames(res_df)
  
  # Split "ENSEMBL|symbol" into ID + symbol (mouse gene name)
  res_df$gene_id   <- sub("\\|.*",  "", res_df$gene_raw)
  res_df$gene_name <- sub(".*\\|", "", res_df$gene_raw)
  
  res_df <- res_df[, c("gene_id", "gene_name",
                       setdiff(colnames(res_df), c("gene_id", "gene_name", "gene_raw")))]
  
  ############################################
  # SUBSET RESULTS FOR EACH GENE LIST (5x)   #
  ############################################
  for (set_name in names(mouse_geneList_sets)) {
    genes_of_interest <- mouse_geneList_sets[[set_name]]
    
    subset_df <- res_df[res_df$gene_name %in% genes_of_interest, ]
    subset_df <- subset_df[order(subset_df$padj), ]
    
    out_file <- file.path(
      OUTPUT_DIR,
      paste0("DESeq2_results_", nameConditionFile, "_", set_name, "_genes_of_interest.tsv")
    )
    
    write.table(subset_df,
                file      = out_file,
                sep       = "\t",
                quote     = FALSE,
                row.names = FALSE)
  }
  
  ##########################
  # SUPERVISED HEATMAP     #
  ##########################
  supervisedTable <- as.data.frame(assay(vsd_mRNA))
  supervisedTable$Genes <- gsub(".*\\|", "", rownames(supervisedTable))
  
  rowNbGeneList <- c()
  for (newGene in newGeneList) {
    nb <- which(supervisedTable$Genes == newGene)
    rowNbGeneList <- c(rowNbGeneList, nb)
  }
  
  # Colors for conditions
  mycolors <- hue_pal()(length(unique(metaData$condition)))
  names(mycolors) <- sort(unique(as.character(metaData$condition)))
  mycolors <- list(condition = mycolors)
  
  HM_supervised_mRNA <- pheatmap(assay(vsd_mRNA)[rowNbGeneList, ],
                                 show_rownames  = TRUE,
                                 scale          = "row",
                                 cluster_rows   = FALSE,
                                 cluster_cols   = FALSE,
                                 annotation_col = metaData,
                                 annotation_colors = mycolors)
  
  png(file.path(OUTPUT_DIR,
                paste0("HeatMap_mRNA_Genes_",
                       format(Sys.time(), "%Y-%m-%d-%H:%M:%S"),
                       ".png")),
      height = 1500, width = 750)
  print(HM_supervised_mRNA)
  dev.off()
}
