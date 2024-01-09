
## This code was written by Rita Bogorad (UGent, the PAE lab) with inputs from : 
## Dr. Willem Stock (UGent, the Phycology lab), Dr. Gust Bilcke (UGent, the PAE lab), and Luz Amadei Martinez (UGent, the PAE lab)
library("DESeq2") # version 1.26.0
library("tidyverse") # version 2.0.0
library("RColorBrewer") # version 1.1.3
library("pheatmap") # version 1.0.12
library("DEGreport") # version 1.22.0
library("tximport") # version 1.14.2
library("ggplot2") # version 3.4.3
library("ggrepel") # version 0.9.3
library("tibble") # version 3.2.1
library("dplyr") # version 1.1.3
library("vegan") # version 2.6.4
library("ggvegan") # version 0.1.0
library("clusterProfiler") # version 3.14.3
library("ggdendroplot") # version 0.1.0
library("VennDiagram") # version 1.7.3

### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.5

#
#### INPUTS ------------------------------------------------------------------------------------------------------  
#
# Loads custom KEGG annotations of the transcriptome
KEGGs <- read.table("only_best_hit_results_diamond_KEGG_SR_transcriptome.txt", sep="\t")
KEGGs <- KEGGs[,c(1:2)]; KEGGs$V2 <- as.character(KEGGs$V2)
KEGGs$V2 <- sapply(strsplit(KEGGs$V2, "\\|"), `[`, 2); colnames(KEGGs) <- c("ID", "KEGG")
# SqueezeMeta annotation file: Tamames & Puente-Sánchez 2019; https://doi.org/10.3389/fmicb.2018.03349
anno_KEGG <- read.table("keggfun2.txt", sep="\t", fill=TRUE, header = FALSE, quote = "")
anno_KEGG$V1 <- as.character(anno_KEGG$V1); anno_KEGG$V4 <- as.character(anno_KEGG$V4)
# Uploads official S.robusta transcriptome annotations: Osuna-Cruz et al 2020; https://doi.org/10.1038/s41467-020-17191-8
mt_anno <- read.table("Semro_gene_info_LATEST.txt", sep="\t", fill=TRUE, header = TRUE, quote = "")

all_signalP6_output <- as.data.frame(read_delim("prediction_results.txt", delim = "\t", skip = 1, show_col_types = FALSE))
all_signalP6_output$ID <- as.character(all_signalP6_output$ID)
all_signalP6_output$ID <- sapply(strsplit(all_signalP6_output$ID, " "), `[`, 1)
transcripts_with_signal_peptide <- all_signalP6_output$ID[all_signalP6_output$Prediction=="SP"]

# Uploads custom COG annotations of the transcriptome
COG_anno_dmnd <- read.table("only_best_hit_results_diamond_COG.txt", sep="\t")
COG_anno <- read.table("eukaryote_orfs.txt")
COG_anno$KEGG.ID <- as.character(gsub("\\*", "", COG_anno$KEGG.ID))
COG_anno <- unique(COG_anno[,c(8,13)])
COG_anno <- COG_anno[COG_anno$COGPATH != "", ]
COG_anno <- COG_anno[ ! grepl("Function unknown", COG_anno$COGPATH), ]
COG_anno <- COG_anno[ ! grepl("\\|", COG_anno$COGPATH), ]
COG_anno <- COG_anno[ ! grepl("General function prediction only", COG_anno$COGPATH), ]
COG_anno <- unique(COG_anno)
table(duplicated(COG_anno$KEGG.ID))

# Low the raw Salmon count data
data <- read.table("Srobusta_transcriptome_raw_salmon_counts.txt", sep="\t")
mean(unlist(data[rowSums(data[,1:4]) > 0,1:4]))
mean(unlist(data[rowSums(data[,5:9]) > 0,5:9]))
mean(unlist(data[rowSums(data[,9:14]) > 0,9:14]))

samples_df_no_outliers_seminvis_response <- data.frame("condition" = c(rep("AX", 4), rep("DPS", 5), rep("DVM", 5)), "Sample" = colnames(data)); rownames(samples_df_no_outliers_seminvis_response) <- colnames(data)

#
#### DE ANALYSIS ------------------------------------------------------------------------------------------------------
#
# Run DEseq2

dds <- DESeqDataSetFromMatrix(countData = data.matrix(data), colData = samples_df_no_outliers_seminvis_response, design = ~ condition)
# Remove low abundant orfs with the count under 10
dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds, normalized=TRUE)) >= 10
dds_filt <- dds[keep,]
table(keep)

## bacterial reads # input from squeezemeta summary table
bac_reads <- c(7824222,	6373524,	8115130,	7364736,	20802708,	37268410,	35534432,	22972202,	27839580,	12741174,	8173216,	35498834,	31718598,	21320156,	21598436,	17949020,	34480676,	40163210,	58756830,	65748872,	73894938)
names(bac_reads) <- c("AX2", "AX3",	"AX4",	"AX5",	"AX6",	"DPS1",	"DPS2",	"DPS3",	"DPS4",	"DPS5",	"DPS6",	"DVM1",	"DVM2",	"DVM4",	"DVM5",	"DVM6",	"PS4",	"PS5",	"PS6",	"VM1",	"VM2")
mean(bac_reads[grepl("DPS", names(bac_reads))]); sd(bac_reads[grepl("DPS", names(bac_reads))])
mean(bac_reads[grepl("DVM", names(bac_reads))]); sd(bac_reads[grepl("DVM", names(bac_reads))]) 

dds_filt$condition=factor(dds_filt$condition, levels=c('AX','DPS', "DVM"))
# Run DESeq2
dds2=DESeq(dds_filt)

## Total number of raw counts per sample
colSums(counts(dds2))
vst_mat <- assay(vst(dds2, blind=TRUE))


# Defining a contrast between DVM and AX conditions for differential expression analysis
contrast_AX_DVM <-  c("condition", "DVM", "AX")
# Performing DESeq2 analysis for the specified contrast
res_table_AX_DVM <- results(dds2, contrast=contrast_AX_DVM, alpha = padj.cutoff)
# Applying log2 fold change shrinkage for more accurate estimates using apeglm method
res_table_AX_DVM_shrink <- lfcShrink(dds2, coef="condition_DVM_vs_AX", res=res_table_AX_DVM, type="apeglm")
# Generating an MA plot using the shrunken log2 fold changes
plotMA(res_table_AX_DVM_shrink, ylim=c(-2,2))

# Converting the results table to a tibble for easier handling and filtering
res_table_AX_DVM_tb <- res_table_AX_DVM_shrink %>%
  data.frame() %>%
  rownames_to_column(var="Gene") %>% 
  as_tibble()

# Filtering significant results based on adjusted p-value and log2 fold change cutoff
sig_AX_DVM <- res_table_AX_DVM_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
# Counting the number of significant genes upregulated and downregulated
dim(sig_AX_DVM[sig_AX_DVM$log2FoldChange > 0.5,])
dim(sig_AX_DVM[sig_AX_DVM$log2FoldChange < 0.5,])
# Ordering the significant genes and classifying them as upregulated or downregulated
sig_AX_DVM_ordered <- sig_AX_DVM[order(sig_AX_DVM$log2FoldChange, sig_AX_DVM$padj),]
sig_AX_DVM_ordered <- as.data.frame(sig_AX_DVM_ordered)
sig_AX_DVM_ordered <- sig_AX_DVM_ordered[order(-sig_AX_DVM_ordered$log2FoldChange, sig_AX_DVM_ordered$padj),]
sig_AX_DVM_ordered$type <- ifelse(sig_AX_DVM_ordered$log2FoldChange > 0, "DVM_up", "DVM_down")
sig_AX_DVM_ordered$KEGG <-  unlist(lapply(sig_AX_DVM_ordered$Gene, function(x) KEGGs$KEGG[KEGGs$ID == x][1]))

# Repeating the process for the contrast between DPS and AX conditions
contrast_AX_DPS <-  c("condition", "DPS", "AX")
res_table_AX_DPS <- results(dds2, contrast=contrast_AX_DPS, alpha = padj.cutoff)
# Applying shrinkage to the DPS vs AX results
res_table_AX_DPS_shrink <- lfcShrink(dds2, coef="condition_DPS_vs_AX", res=res_table_AX_DPS, type="apeglm")
# Generating an MA plot for the DPS vs AX results
plotMA(res_table_AX_DPS_shrink, ylim=c(-2,2))

# Preparing the results table for DPS vs AX for analysis
res_table_AX_DPS_tb <- res_table_AX_DPS_shrink %>%
  data.frame() %>%
  rownames_to_column(var="Gene") %>% 
  as_tibble()
# Filtering significant genes for DPS vs AX
sig_AX_DPS <- res_table_AX_DPS_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

# Counting significant upregulated and downregulated genes in DPS vs AX
dim(sig_AX_DPS[sig_AX_DPS$log2FoldChange > 0.5,])
dim(sig_AX_DPS[sig_AX_DPS$log2FoldChange < 0.5,])

# Ordering and classifying significant genes in the DPS vs AX comparison
sig_AX_DPS_ordered <- sig_AX_DPS[order(sig_AX_DPS$log2FoldChange, sig_AX_DPS$padj),]
sig_AX_DPS_ordered <- as.data.frame(sig_AX_DPS_ordered)
sig_AX_DPS_ordered <- sig_AX_DPS_ordered[order(-sig_AX_DPS_ordered$log2FoldChange, sig_AX_DPS_ordered$padj),]
sig_AX_DPS_ordered$type <- ifelse(sig_AX_DPS_ordered$log2FoldChange > 0, "DPS_up", "DPS_down")
sig_AX_DPS_ordered$KEGG <-  unlist(lapply(sig_AX_DPS_ordered$Gene, function(x) KEGGs$KEGG[KEGGs$ID == x][1]))

shared_DE_AX <- Reduce(intersect, list(sig_AX_DVM_ordered$Gene[sig_AX_DVM_ordered$log2FoldChange < -0.5 & sig_AX_DVM_ordered$KEGG != ""], sig_AX_DPS_ordered$Gene[sig_AX_DPS_ordered$log2FoldChange < -0.5 & sig_AX_DPS_ordered$KEGG != ""]))

#
#### DE ANALYSIS (quality control and plots) ------------------------------------------------------------------------------------------------------
#
# perform hierarchical clustering
dist_vst_mat <- dist(t(vst_mat), method = "euclidean")
sampleDistMatrix <- as.matrix( dist_vst_mat )
colors <- colorRampPalette( brewer.pal(9, "Blues"))(255)
eucl_dist_heatmap <- pheatmap(sampleDistMatrix,
                              clustering_distance_rows = dist_vst_mat,
                              clustering_distance_cols = dist_vst_mat,
                              color = colors)

dist_vst_mat <- as.data.frame(as.matrix(dist_vst_mat))
mean(unlist(dist_vst_mat[1:4,1:4])[unlist(dist_vst_mat[1:4,1:4]) > 0 ])
sd(unlist(dist_vst_mat[1:4,1:4])[unlist(dist_vst_mat[1:4,1:4]) > 0 ])
mean(unlist(dist_vst_mat[5:9,5:9])[unlist(dist_vst_mat[5:9,5:9]) > 0 ])
sd(unlist(dist_vst_mat[5:9,5:9])[unlist(dist_vst_mat[5:9,5:9]) > 0 ])
mean(unlist(dist_vst_mat[10:14,10:14])[unlist(dist_vst_mat[10:14,10:14]) > 0 ])
sd(unlist(dist_vst_mat[10:14,10:14])[unlist(dist_vst_mat[10:14,10:14]) > 0 ])

## VENN DIAGRAM NUMBER 1
vst_mat_annotated <- vst_mat
data_filtered_annotated <- data[rownames(data) %in% rownames(vst_mat_annotated), ]
data_filtered_annotated$KEGG <-  unlist(lapply(rownames(data_filtered_annotated), function(x) KEGGs$KEGG[KEGGs$ID == x][1]))
data_filtered_annotated$ID <- rownames(data_filtered_annotated)
# Filtering data_filtered_annotated for non-NA KEGG values and non-zero row sums across the first 4 columns, then selecting columns 16 and 15
total_KEGG_AX <- data_filtered_annotated[!is.na(data_filtered_annotated$KEGG) & rowSums(data_filtered_annotated[,1:4]) > 0, c(16, 15)]
# Removing duplicate rows based on the KEGG column in the total_KEGG_AX dataframe
total_KEGG_AX <- total_KEGG_AX[!duplicated(total_KEGG_AX$KEGG), ]
# Similar to the first operation, but for a different set of columns (1:5) to filter for the DPS condition
total_KEGG_DPS <- data_filtered_annotated[!is.na(data_filtered_annotated$KEGG) & rowSums(data_filtered_annotated[,1:5]) > 0, c(16, 15)]
# Removing duplicates in the total_KEGG_DPS dataframe
total_KEGG_DPS <- total_KEGG_DPS[!duplicated(total_KEGG_DPS$KEGG), ]
# Filtering for the DVM condition by selecting non-NA KEGG values and non-zero row sums across columns 6 to 10
total_KEGG_DVM <- data_filtered_annotated[!is.na(data_filtered_annotated$KEGG) & rowSums(data_filtered_annotated[,6:10]) > 0, c(16, 15)]
# Removing duplicates in the total_KEGG_DVM dataframe
total_KEGG_DVM <- total_KEGG_DVM[!duplicated(total_KEGG_DVM$KEGG), ]

venn.diagram(
  x = list(unique(total_KEGG_AX$KEGG), unique(total_KEGG_DPS$KEGG), unique(total_KEGG_DVM$KEGG)),
  category.names = c("AX", "DPS","DVM"),
  main.fontfamily = "calibri",
  main = "Upregulated DE Transcripts",
  main.cex = 2, cat.cex = 2, inverted=FALSE, cat.pos = c(140,140,900),
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = as.character(c("#F8766D","#00BA38", "#619CFF")),
  
  # Numbers
  cex = 3,
  fontface = "bold",
  fontfamily = "calibri",
  filename = 'venn_diagramm_all_KEGGs.png',
  output=TRUE
)

## VENN DIAGRAM NUMBER 2 AND 3 - DE
# Extracting the names of genes upregulated in the DVM condition and omitting NA values
names_sigs_AX_DVM_upregulated <- na.omit(sig_AX_DVM_ordered$Gene[sig_AX_DVM_ordered$log2FoldChange > 0])

# Extracting the names of genes upregulated in the DPS condition and omitting NA values
names_sigs_AX_DPS_upregulated <- na.omit(sig_AX_DPS_ordered$Gene[sig_AX_DPS_ordered$log2FoldChange > 0])

# Creating a Venn diagram to visualize the overlap of upregulated genes between DVM and DPS conditions
venn.diagram(
  x = list(names_sigs_AX_DVM_upregulated, names_sigs_AX_DPS_upregulated),
  category.names = c("DVM", "DPS"),
  main.fontfamily = "calibri",
  main = "Upregulated DE Transcripts",
  main.cex = 2, cat.cex = 2, cat.pos = c(140, 900), inverted = TRUE,
  lwd = 2,
  lty = 'blank',
  fill = as.character(c("steelblue1", "#00BA38")),
  cex = 3,
  fontface = "bold",
  fontfamily = "calibri",
  filename = 'venn_diagramm_DVM_and_DPS_upregulated.png',
  output = TRUE
)

# Extracting the names of genes downregulated in the DVM condition and omitting NA values
names_sigs_AX_DVM_downregulated <- na.omit(sig_AX_DVM_ordered$Gene[sig_AX_DVM_ordered$log2FoldChange < 0])

# Extracting the names of genes downregulated in the DPS condition and omitting NA values
names_sigs_AX_DPS_downregulated <- na.omit(sig_AX_DPS_ordered$Gene[sig_AX_DPS_ordered$log2FoldChange < 0])

# Creating a Venn diagram to visualize the overlap of downregulated genes between DVM and DPS conditions
venn.diagram(
  x = list(names_sigs_AX_DVM_downregulated, names_sigs_AX_DPS_downregulated),
  category.names = c("DVM", "DPS"),
  main.fontfamily = "calibri",
  main = "Downregulated DE Transcripts",
  main.cex = 2, cat.cex = 2, cat.pos = c(140, 900), inverted = TRUE,
  lwd = 2,
  lty = 'blank',
  fill = as.character(c("steelblue1", "#00BA38")),
  cex = 3,
  fontface = "bold",
  fontfamily = "calibri",
  filename = 'venn_diagramm_DVM_and_DPS_downregulated.png',
  output = TRUE
)

## PCA
# Performing PCA on the transposed variance-stabilized transformed matrix
vst_mat.pca <- rda(t(vst_mat))

# Extracting PCA scores and adding the corresponding condition from samples_df_no_outliers_seminvis_response
vscores <- data.frame(vst_mat.pca$CA$u)
vscores$Condition <- samples_df_no_outliers_seminvis_response$condition

# Calculating the percentage of variance explained by the first two principal components
constrained_eigPCA <- eigenvals(vst_mat.pca) / sum(eigenvals(vst_mat.pca)) * 100

# Renaming columns in vscores to include the percentage of variance explained
colnames(vscores)[1] <- paste0(colnames(vscores)[1], " (", round(constrained_eigPCA[1]), "%)")
colnames(vscores)[2] <- paste0(colnames(vscores)[2], " (", round(constrained_eigPCA[2]), "%)")

# Creating a ggplot of the PCA results
ggplot(vscores, aes(x = `PC1 (19%)`, y = `PC2 (14%)`, col = Condition, shape = Condition)) + 
  geom_point(size = 5.5) +
  geom_text(data = vscores, aes(x = `PC1 (19%)`, y = `PC2 (14%)`, label = as.numeric(gsub("\\D", "", rownames(vscores)))), size = 2.5, color = "black") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0), legend.position = "bottom")

#
#### DIATOM SECRETION (interpreted from the SignalP and DeepLoc outputs) ------------------------------------------------------------------------------------------------------
#
# Loads SignalP and NetGPI outputs
# Converting the 'Protein_ID' and 'Localizations' columns in Deep_Loc_out to character type for processing
Deep_Loc_out$Protein_ID <- as.character(Deep_Loc_out$Protein_ID)
Deep_Loc_out$Localizations <- as.character(Deep_Loc_out$Localizations)

# Reading in the Net_GPI output file and converting relevant columns to character type
Net_GPI_out <- read.table("Net_GPI_all_sro_with_SP_prediction.txt", header=TRUE, sep="\t")
Net_GPI_out$ID <- as.character(Net_GPI_out$ID)
Net_GPI_out$Pred..GPI.Anchored <- as.character(Net_GPI_out$Pred..GPI.Anchored)

# Mapping SignalP predictions to the genes in sig_AX_DVM_ordered
sig_AX_DVM_ordered$SP <- unlist(lapply(sig_AX_DVM_ordered$Gene, function(x) all_signalP6_output$Prediction[all_signalP6_output$ID == x]))

# Reading in the DeepLoc output for DVM and converting relevant columns to character type
Deep_Loc_DVM_out <- read.table("DVM_DeepLoc.csv", sep=",", header=TRUE)
Deep_Loc_DVM_out$Protein_ID <- as.character(Deep_Loc_DVM_out$Protein_ID)
Deep_Loc_DVM_out$Localizations <- as.character(Deep_Loc_DVM_out$Localizations)

# Mapping DeepLoc localizations to the genes in sig_AX_DVM_ordered
sig_AX_DVM_ordered$DeepLoc <- unlist(lapply(sig_AX_DVM_ordered$Gene, function(x) {
  if (x %in% Deep_Loc_DVM_out$Protein_ID) {
    Deep_Loc_DVM_out$Localizations[Deep_Loc_DVM_out$Protein_ID == x]
  } else {
    NA
  }
}))

# Mapping Net_GPI predictions to the genes in sig_AX_DVM_ordered
sig_AX_DVM_ordered$Anch <- unlist(lapply(sig_AX_DVM_ordered$Gene, function(x) {
  if (x %in% Net_GPI_out$ID) {
    Net_GPI_out$Pred..GPI.Anchored[Net_GPI_out$ID == x]
  } else {
    NA
  }
}))

# Mapping SignalP predictions to the genes in sig_AX_DPS_ordered
sig_AX_DPS_ordered$SP <- unlist(lapply(sig_AX_DPS_ordered$Gene, function(x) all_signalP6_output$Prediction[all_signalP6_output$ID == x]))

# Reading in the DeepLoc output for DPS and converting relevant columns to character type
Deep_Loc_DPS_out <- read.table("DPS_DeepLoc.csv", sep=",", header=TRUE)
Deep_Loc_DPS_out$Protein_ID <- as.character(Deep_Loc_DPS_out$Protein_ID)
Deep_Loc_DPS_out$Localizations <- as.character(Deep_Loc_DPS_out$Localizations)

# Mapping DeepLoc localizations to the genes in sig_AX_DPS_ordered
sig_AX_DPS_ordered$DeepLoc <- unlist(lapply(sig_AX_DPS_ordered$Gene, function(x) {
  if (x %in% Deep_Loc_DPS_out$Protein_ID) {
    Deep_Loc_DPS_out$Localizations[Deep_Loc_DPS_out$Protein_ID == x]
  } else {
    NA
  }
}))

# Mapping Net_GPI predictions to the genes in sig_AX_DPS_ordered
sig_AX_DPS_ordered$Anch <- unlist(lapply(sig_AX_DPS_ordered$Gene, function(x) {
  if (x %in% Net_GPI_out$ID) {
    Net_GPI_out$Pred..GPI.Anchored[Net_GPI_out$ID == x]
  } else {
    NA
  }
}))



#
#### COG visualization of DE transcripts ------------------------------------------------------------------------------------------------------
#
# Preparing the DVM ordered dataset for figure creation
sig_AX_DVM_ordered_fig <- sig_AX_DVM_ordered
# Mapping COG pathways to KEGG IDs in sig_AX_DVM_ordered_fig
sig_AX_DVM_ordered_fig$COG <- unlist(lapply(sig_AX_DVM_ordered_fig$KEGG, function(x) COG_anno$COGPATH[COG_anno$KEGG.ID == x][1]))

# Converting COG column to character and removing NA entries
sig_AX_DVM_ordered_fig$COG <- as.character(sig_AX_DVM_ordered_fig$COG)
sig_AX_DVM_ordered_fig <- sig_AX_DVM_ordered_fig[!is.na(sig_AX_DVM_ordered_fig$COG), ]

# Converting COG column to a factor with levels sorted in ascending order of frequency
sig_AX_DVM_ordered_fig$COG <- factor(sig_AX_DVM_ordered_fig$COG, levels=names(table(sig_AX_DVM_ordered_fig$COG) %>% sort(decreasing = FALSE)))

# Recoding the 'type' column in sig_AX_DVM_ordered_fig
sig_AX_DVM_ordered_fig$type <- ifelse(sig_AX_DVM_ordered_fig$type == "DVM_up", "DVM", "AX")

# Repeating the process for the DPS ordered dataset
sig_AX_DPS_ordered_fig <- sig_AX_DPS_ordered
sig_AX_DPS_ordered_fig$COG <- unlist(lapply(sig_AX_DPS_ordered_fig$KEGG, function(x) COG_anno$COGPATH[COG_anno$KEGG.ID == x][1]))
sig_AX_DPS_ordered_fig$COG <- as.character(sig_AX_DPS_ordered_fig$COG)
sig_AX_DPS_ordered_fig <- sig_AX_DPS_ordered_fig[!is.na(sig_AX_DPS_ordered_fig$COG), ]
sig_AX_DPS_ordered_fig$COG <- factor(sig_AX_DPS_ordered_fig$COG, levels=names(table(sig_AX_DPS_ordered_fig$COG) %>% sort(decreasing = FALSE)))
sig_AX_DPS_ordered_fig$type <- ifelse(sig_AX_DPS_ordered_fig$type == "DPS_up", "DPS", "AX")

# Combining DVM and DPS datasets and adding subtype for diatoms
sig_AX_DVM_DPS_fig <- sig_AX_DVM_ordered_fig[sig_AX_DVM_ordered_fig$type == "DVM" | sig_AX_DVM_ordered_fig$Gene %in% shared_DE_AX,]
sig_AX_DVM_DPS_fig <- rbind(sig_AX_DVM_DPS_fig, sig_AX_DPS_ordered_fig[sig_AX_DPS_ordered_fig$type == "DPS", ])
sig_AX_DVM_DPS_fig$subtype <- "Diatom"

# Reading in bacterial COG data and setting subtype for bacteria
sig_bacteria_DPS_DVM_KEGG_fig <- readRDS("sig_DPS_DVM_KEGG_fig.rds")
sig_bacteria_DPS_DVM_KEGG_fig$subtype <- "Bacteria"

# Combining diatom and bacterial datasets
all_COG_fig <- rbind(sig_AX_DVM_DPS_fig[,c(1,3,11,7,12,13)], sig_bacteria_DPS_DVM_KEGG_fig[,c(1,3,9,10,12,13)])
all_COG_fig$subtype <- factor(all_COG_fig$subtype, levels=c("Diatom", "Bacteria"))

# Adjusting log2FoldChange values for visualization
all_COG_fig$log2FoldChange[all_COG_fig$log2FoldChange < 0] <- all_COG_fig$log2FoldChange[all_COG_fig$log2FoldChange < 0] * -1

# Creating a ggplot of differentially expressed transcripts categorized by COG function
COG_DE_diatom <- ggplot(all_COG_fig, aes(x = log2FoldChange, y = COG, color=type)) +
  geom_point()  + 
  facet_nested(.~subtype+type, scales = 'free') + 
  xlab("log2 fold change") + 
  ylab("") + 
  labs(fill = "") + theme_bw() + 
  guides(fill=guide_legend(ncol=5)) + theme(legend.position="none", 
                                            axis.text.x = element_text(angle = 0), 
                                            panel.grid = element_line(color = "white"), strip.background =element_rect(fill="white"),
                                            strip.text.y = element_text(angle = 0))  +
  scale_color_manual(values=setNames(c("#F8766D", "steelblue1", "#00BA38"), c("AX", "DVM", "DPS")))+ 
  scale_x_continuous(breaks= pretty_breaks()) +
  ggtitle("DE transcripts - COG function")

#
#### ENRICHMENT ANALYSIS ------------------------------------------------------------------------------------------------------
# 
# Performing gene set enrichment analysis for upregulated genes in DVM condition
enrich_DVM <- enricher(as.character(unique(sig_AX_DVM_ordered$Gene[sig_AX_DVM_ordered$log2FoldChange > 0])), TERM2GENE = data_filtered_annotated[, c(15:16)])
# Displaying the top rows of the enrichment results
head(enrich_DVM)
# Converting the enrichment results to a data frame for easier manipulation
enrich_DVM <- as.data.frame(enrich_DVM)
# Mapping KEGG pathway details to the enrichment results
enrich_DVM$KEGG_path <- unlist(lapply(enrich_DVM$ID, function(x) anno_KEGG$V4[anno_KEGG$V1 == x][1]))
enrich_DVM$KEGG_name <- unlist(lapply(enrich_DVM$ID, function(x) anno_KEGG$V3[anno_KEGG$V1 == x][1]))
# Selecting specific columns from the enrichment results and adding a new column for type
enrich_DVM_sort <- enrich_DVM[, c(1:7, 10, 11)]
enrich_DVM_sort$type <- "DVM"

# Performing gene set enrichment analysis for upregulated genes in DPS condition
enrich_DPS <- enricher(as.character(unique(sig_AX_DPS_ordered$Gene[sig_AX_DPS_ordered$log2FoldChange > 0])), TERM2GENE = data_filtered_annotated[, c(15:16)])
# Displaying the top rows of the enrichment results
head(enrich_DPS)

# Converting the enrichment results to a data frame
enrich_DPS <- as.data.frame(enrich_DPS)

# Mapping KEGG pathway details to the enrichment results
enrich_DPS$KEGG_path <- unlist(lapply(enrich_DPS$ID, function(x) anno_KEGG$V4[anno_KEGG$V1 == x][1]))
enrich_DPS$KEGG_name <- unlist(lapply(enrich_DPS$ID, function(x) anno_KEGG$V3[anno_KEGG$V1 == x][1]))

# Selecting specific columns from the enrichment results and cleaning the KEGG pathway names
enrich_DPS_sort <- enrich_DPS[, c(1:7, 10, 11)]
enrich_DPS_sort$KEGG_path <- gsub("Protein families: ", "", enrich_DPS_sort$KEGG_path)
enrich_DPS_sort$KEGG_path <- gsub("Brite Hierarchies; ", "", enrich_DPS_sort$KEGG_path)
enrich_DPS_sort$KEGG_path <- str_to_title(enrich_DPS_sort$KEGG_path)

# Combining enrichment results from DVM and DPS
enriched_DE <- rbind(enrich_DVM, enrich_DPS)

# Extracting unique transcripts from the combined enrichment results
vec_enriched_transcripts <- unique(unlist(strsplit(enriched_DE$geneID, "/")))


# downregulated SIG
enrich_AX_DVM <- enricher(as.character(unique(sig_AX_DVM_ordered$Gene[sig_AX_DVM_ordered$log2FoldChange<0])), TERM2GENE=data_filtered_annotated[,c(15:16)])
head(enrich_AX_DVM)
enrich_AX_DVM <- as.data.frame(enrich_AX_DVM)
enrich_AX_DVM$KEGG_name <- unlist(lapply(enrich_AX_DVM$ID, function(x) anno_KEGG$V3[anno_KEGG$V1 == x][1]))

enrich_AX_DPS <- enricher(as.character(unique(sig_AX_DPS_ordered$Gene[sig_AX_DPS_ordered$log2FoldChange<0])), TERM2GENE=data_filtered_annotated[,c(15:16)])
head(enrich_AX_DPS)

AX_intersection <- Reduce(intersect, list(unique(sig_AX_DVM_ordered$Gene[sig_AX_DVM_ordered$log2FoldChange<0]), unique(sig_AX_DPS_ordered$Gene[sig_AX_DPS_ordered$log2FoldChange<0])))
AX_intersect_vst <- vst_mat[rownames(vst_mat) %in% AX_intersection, ]
AX_intersect_vst_long <- reshape2::melt(AX_intersect_vst)
AX_intersect_vst_long$KEGG <-  unlist(lapply(AX_intersect_vst_long$Var1, function(x) sig_AX_DVM_ordered$KEGG[sig_AX_DVM_ordered$Gene == x][1]))
AX_intersect_vst_long$KEGG_name <-  unlist(lapply(AX_intersect_vst_long$KEGG, function(x) anno_KEGG$V3[anno_KEGG$V1 == x][1]))
AX_intersect_vst_long$condition <- "AX"
AX_intersect_vst_long$condition[grep("DPS", AX_intersect_vst_long$Var2)] <- "DPS"
AX_intersect_vst_long$condition[grep("DVM", AX_intersect_vst_long$Var2)] <- "DVM"
ggplot(AX_intersect_vst_long, aes(  Var2,Var1,fill= value)) + 
  geom_tile()   + facet_grid(.~condition, scales="free",space = "free", labeller = label_wrap_gen()) + theme_bw() +
  labs(y="Axenic Shared Upregulated Transcripts")  + 
  #guides(fill=guide_legend("VST Count"),  legend.key.width = unit(2, 'cm')) +
  theme(strip.background =element_rect(fill="white"), legend.position = "bottom", axis.text.y=element_text(size=7), strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.title.y = element_blank(),axis.title.x = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), panel.grid = element_line(color = "white"))  + # axis.text.y=element_blank(), axis.ticks.y=element_blank(),
  scale_fill_gradient2(midpoint = mean(AX_intersect_vst_long$value), low = "slateblue4", mid = "lemonchiffon", high = "orangered", space = "Lab" )




#
#### ENRICHMENT ANALYSIS - FIGURE ------------------------------------------------------------------------------------------------------
# 
# figure heatmap
# Subsetting data_filtered_annotated for transcripts in vec_enriched_transcripts
enriched_DE_vst <- data_filtered_annotated[rownames(data_filtered_annotated) %in% vec_enriched_transcripts, ]

# Reshaping the data to a long format for visualization
enriched_DE_vst_long <- reshape2::melt(enriched_DE_vst, id.vars = c("KEGG", "ID"))

# Assigning conditions based on the variable names in the reshaped data
enriched_DE_vst_long$condition <- "AX"
enriched_DE_vst_long$condition[grep("DPS", enriched_DE_vst_long$variable)] <- "DPS"
enriched_DE_vst_long$condition[grep("DVM", enriched_DE_vst_long$variable)] <- "DVM"

# Setting KEGG column as a factor and ordering by KEGG pathways
enriched_DE_vst_long$KEGG <- factor(enriched_DE_vst_long$KEGG, levels = c(unique(c(enrich_DPS$ID, enrich_DVM$ID))))
enriched_DE_vst_long <- enriched_DE_vst_long[order(enriched_DE_vst_long$KEGG),]

# Defining a function to scale values
scaleG <- function(vector){
  return(vector / max(vector))
}

# Applying the scaling function to the data
enriched_DE_vst_long_scaled <- enriched_DE_vst_long
for(i in 1:length(unique(enriched_DE_vst_long_scaled$ID))){
  enriched_DE_vst_long_scaled$value[enriched_DE_vst_long_scaled$ID == unique(enriched_DE_vst_long_scaled$ID)[i]] <- scaleG(enriched_DE_vst_long$value[enriched_DE_vst_long$ID == unique(enriched_DE_vst_long$ID)[i]])
}

# Creating a heatmap of the scaled data
seminavis_enrich <- ggplot(enriched_DE_vst_long_scaled, aes(  variable,ID,fill= value)) + 
  geom_tile()   + facet_grid(KEGG~condition, scales="free",space = "free", labeller = label_wrap_gen()) + theme_bw() +
  labs(y="Significant DE Transcript")  + 
  #guides(fill=guide_legend("VST Count"),  legend.key.width = unit(2, 'cm')) +
  theme(strip.background =element_rect(fill="white"), legend.position = "bottom", axis.text.y=element_text(size=7), strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.title.y = element_blank(),axis.title.x = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), panel.grid = element_line(color = "white"))  + # axis.text.y=element_blank(), axis.ticks.y=element_blank(),
  scale_fill_gradient2(midpoint = mean(enriched_DE_vst_long_scaled$value), low = "slateblue4", mid = "lemonchiffon", high = "orangered", space = "Lab" )

# re-coloring the panels
g <- ggplot_gtable(ggplot_build(seminavis_enrich))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("lightgoldenrod1","lightgoldenrod1","lightgoldenrod1","lightgoldenrod1", "lightgoldenrod1", "white", "white", "white", "white")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)



# figure annotation panel
tr_to_gene <- read.table("transcript_mapping.all_transcripts.sro.csv", sep=";")
colnames(tr_to_gene) <- c("transcript", "gene")
enriched_DE_vst_annotated <- enriched_DE_vst[,c(15,16)]
enriched_DE_vst_annotated$gene_id <-  unlist(lapply(enriched_DE_vst_annotated$ID, function(x) tr_to_gene$gene[tr_to_gene$transcript == x][1]))
enriched_DE_vst_annotated$gene_id <- as.character(enriched_DE_vst_annotated$gene_id)
enriched_DE_vst_annotated$KEGG_name <- unlist(lapply(enriched_DE_vst_annotated$KEGG, function(x) anno_KEGG$V3[anno_KEGG$V1 == x][1]))
enriched_DE_vst_annotated$mt_anno <- unlist(lapply(enriched_DE_vst_annotated$ID, function(x) mt_anno$description[mt_anno$Gene.ID == x][1]))
enriched_DE_vst_annotated$KEGG_name <- sapply(strsplit(sapply(strsplit(enriched_DE_vst_annotated$KEGG_name, " \\["), `[`, 1), ", other"), `[`, 1)
enriched_DE_vst_annotated$KEGG_name <- gsub("general ", "", enriched_DE_vst_annotated$KEGG_name)
enriched_DE_vst_annotated$ID <- factor(enriched_DE_vst_annotated$ID, levels=unique(enriched_DE_vst_long$ID))
enriched_DE_vst_annotated$mt_anno <- as.character(enriched_DE_vst_annotated$mt_anno)
enriched_DE_vst_annotated$mt_anno <- sapply(strsplit(enriched_DE_vst_annotated$mt_anno, "\\-19"), `[`, 1)
enriched_DE_vst_annotated$mt_anno <- sapply(strsplit(enriched_DE_vst_annotated$mt_anno, " 4"), `[`, 1)
enriched_DE_vst_annotated$mt_anno <- str_to_title(enriched_DE_vst_annotated$mt_anno)
enriched_DE_vst_annotated$mt_anno <- factor(enriched_DE_vst_annotated$mt_anno, levels=names(table(enriched_DE_vst_annotated$mt_anno)[order(table(enriched_DE_vst_annotated$mt_anno),decreasing=TRUE)]))
enriched_DE_vst_annotated$KEGG <- factor(enriched_DE_vst_annotated$KEGG, levels=c(unique(c(enrich_DPS$ID, enrich_DVM$ID))))
enriched_DE_vst_annotated <- enriched_DE_vst_annotated[order(enriched_DE_vst_annotated$KEGG),]
enriched_DE_vst_annotated$KEGG_name <- factor(enriched_DE_vst_annotated$KEGG_name, levels=unique(enriched_DE_vst_annotated$KEGG_name))

seminavis_annotation <- ggplot(enriched_DE_vst_annotated, aes(  mt_anno,ID)) + 
  geom_point()   + facet_grid(KEGG_name~., scales="free",space = "free", labeller = label_wrap_gen()) + theme_bw() +
  labs(y="Significant DE Transcript")  + 
  theme(strip.background =element_rect(fill="white"), legend.position = "bottom", axis.text.y=element_text(size=7), strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 35, hjust=1),axis.title.y = element_blank(),axis.title.x = element_blank(), axis.ticks.y=element_blank(), panel.grid = element_line(color = "white")) # axis.text.y=element_blank(), axis.ticks.y=element_blank(),

# re-coloring the panels
g2 <- ggplot_gtable(ggplot_build(seminavis_annotation))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("lightgoldenrod1","lightgoldenrod1","lightgoldenrod1","lightgoldenrod1", "lightgoldenrod1", "white", "white", "white", "white")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g2)

