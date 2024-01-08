
## This code was written by Rita Bogorad (UGent, the PAE lab) with inputs from : 
## Dr. Willem Stock (UGent, the Phycology lab), Dr. Gust Bilcke (UGent, the PAE lab), and Luz Amadei Martinez (UGent, the PAE lab)
library("tibble") # version 3.1.8
library("DESeq2") # version 1.26.0
library("vegan") # version 2.6-4
library("ggvegan") # version 0.1-0
library("pheatmap") # version 1.0.12
library("VennDiagram") # version 1.7.3
library("clusterProfiler") # version 3.14.3
library("ggh4x") # version 0.2.3
library("ggdendroplot") # version 0.1.0
library("pheatmap") # version 1.0.12

set.seed(10)

# Assignes consistent colors to the conditions
colors <- c("#F8766D","#00BA38", "#619CFF", "gold2", "skyblue"); names(colors) <- c("AX", "DPS", "DVM", "PS", "VM")
colors <- colors[order(names(colors))]
shape <- c(16,17,15, 18, 18); names(shape) <- c("AX", "DPS", "DVM", "PS", "VM")

### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.5

#
#### INPUTS ------------------------------------------------------------------------------------------------------  
#
# Loads the tab-separated file "bacterial_orfs.txt", a squeezemeta output, into a data frame
all_bacterial_data <- read.table("bacterial_orfs.txt", sep="\t")
# Replaces all instances of "*" in the "KEGG.ID" column of "all_bacterial_data" with an empty string
all_bacterial_data$KEGG.ID <- gsub("\\*", "", all_bacterial_data$KEGG.ID)
all_bacterial_data$KEGGFUN <- as.character(all_bacterial_data$KEGGFUN)
all_bacterial_data$KEGGPATH <- as.character(all_bacterial_data$KEGGPATH)
# Extracts columns 58 to 78 with the experiment raw counts from "all_bacterial_data"
bacterial_raw_count <- all_bacterial_data[,c(58:78)]
bacterial_raw_count <- as.data.frame(bacterial_raw_count)
# Removes all rows from "bacterial_raw_count" where the sum of the values in each row is equal to 0
bacterial_raw_count <- bacterial_raw_count[rowSums(bacterial_raw_count) !=0,]
colnames(bacterial_raw_count) <- gsub("Raw.read.count.", "", colnames(bacterial_raw_count))
# Removes the outliers AX6 and DPS6 
bacterial_raw_count_no_outlier <- bacterial_raw_count[,-c(1:5,11,17:21)]
# Creates the metadata file for the DE analysis
samples <- data.frame("condition" = c(rep("DPS", 5), rep("DVM", 5)), "Sample" = colnames(bacterial_raw_count_no_outlier)); rownames(samples) <- colnames(bacterial_raw_count_no_outlier)

#### summary statistics of the co-cultures
dim(bacterial_raw_count_no_outlier)
mean(colSums(bacterial_raw_count_no_outlier)[1:5]); sd(colSums(bacterial_raw_count_no_outlier)[1:5])
mean(colSums(bacterial_raw_count_no_outlier)[6:10]); sd(colSums(bacterial_raw_count_no_outlier)[6:10])

# Calculating the mean of counts for different conditions
mean(unlist(bacterial_raw_count_no_outlier[rowSums(bacterial_raw_count_no_outlier[,1:5]) > 0,1:5]))
mean(unlist(bacterial_raw_count_no_outlier[rowSums(bacterial_raw_count_no_outlier[,6:10]) > 0,6:10]))

#
#### PCA (co-culture and bacterial control samples) ------------------------------------------------------------------------------------------------------
#
# Creates a table without the AX control for the PCA (because this dataset is only for the analysis of the bacterial community)
bacterial_conrols_counts <- bacterial_raw_count[,-c(1:5,11)]
samples_bac_control <- data.frame("condition" = c(rep("DPS", 5), rep("DVM", 5), rep("PS", 3), rep("VM", 2)), "Sample" = colnames(bacterial_conrols_counts)); rownames(samples_bac_control) <- colnames(bacterial_conrols_counts)

# Creates the DESeq2 object to calculate VST counts, low abundant reads ( < 10 per sample) are NOT excluded
dds_bac_controls <- DESeqDataSetFromMatrix(countData = data.matrix(bacterial_conrols_counts), colData = samples_bac_control, design = ~ condition)
dds2_bac_controls=DESeq(dds_bac_controls)
vst_bac_controls_count <- vst(dds2_bac_controls, blind=T)
vst_bac_controls_mat <- assay(vst_bac_controls_count)

# Uses VST-normalized counts of the co-cultures and the bacterial controls to calculate PCA
vst_bac_controls_mat.pca <- rda(t(vst_bac_controls_mat))
vscores <- data.frame(vst_bac_controls_mat.pca$CA$u)
vscores$Condition <- samples_bac_control$condition

# Calculate the percentage of variance explained by each principal component
constrained_eigPCA<-eigenvals(vst_bac_controls_mat.pca)/ sum(eigenvals(vst_bac_controls_mat.pca))*100
colnames(vscores)[1] <- paste0(colnames(vscores)[1], " (", round(constrained_eigPCA[1]), "%)")
colnames(vscores)[2] <- paste0(colnames(vscores)[2], " (", round(constrained_eigPCA[2]), "%)")

# Plots calculated PCA
ggplot(vscores, aes(x = `PC1 (26%)`, y = `PC2 (13%)`, col = Condition, shape = Condition)) + 
  geom_point(size=5.5) + scale_color_manual(values = colors) + scale_shape_manual(values = shape) + 
  geom_text(data = vscores, aes(x = `PC1 (26%)`, y = `PC2 (13%)`, label = as.numeric(gsub("\\D", "", rownames(vscores)))), size=2.5, color="black") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0), legend.position = "bottom")


#
#### DE ANALYSIS (co-culture samples only) ------------------------------------------------------------------------------------------------------
#
dds <- DESeqDataSetFromMatrix(countData = data.matrix(bacterial_raw_count_no_outlier), colData = samples, design = ~ condition)
# Remove low abundant orfs (< 10 reads in the sum of all samples are excluded)
keep = rowSums(counts(dds)) >= 10
table(keep)
dds_filt = dds[keep,]

# Choose factor levels
dds_filt$condition=factor(dds_filt$condition, levels=c('DPS', "DVM"))
# Run DESeq2
dds2=DESeq(dds_filt)
resultsNames(dds2)
#results=results(dds2, name='condition_DVM_vs_DPS') # upregulated in DVM and downregulated in DPS
plotDispEsts(dds2)
# Generate a table with VST-normalized counts
vst_count <- vst(dds2, blind=T)
vst_mat <- assay(vst_count)
# Annotate the VST matrix with KEGG IDs
vst_mat_annotated <- vst_mat
vst_mat_annotated <- as.data.frame(vst_mat_annotated)
vst_mat_annotated$KEGG <-  unlist(lapply(rownames(vst_mat_annotated), function(x) all_bacterial_data$KEGG.ID[rownames(all_bacterial_data) == x][1]))
vst_mat_annotated$ID <- rownames(vst_mat_annotated)

# Calculates the percentage of annotated transcripts out of the total filtered transcripts
round(nrow(vst_mat_annotated[vst_mat_annotated$KEGG != "",])*100/nrow(vst_mat_annotated))

contrast_DVM_DPS <-  c("condition", "DVM", "DPS")
res_table_DVM_DPS <- results(dds2, contrast=contrast_DVM_DPS, alpha = padj.cutoff, parallel = TRUE)
res_table_DVM_DPS_tb <- res_table_DVM_DPS %>%
  data.frame() %>%
  rownames_to_column(var="Gene") %>% 
  as_tibble()

# shrinkage
res_table_DPS_DVM_shrink <- lfcShrink(dds2, coef="condition_DVM_vs_DPS", res=res_table_DVM_DPS, type="ashr")

# MA plot using shrunken fold changes
plotMA(res_table_DPS_DVM_shrink, ylim=c(-2,2))

res_table_DPS_DVM_shrink_tb <- res_table_DPS_DVM_shrink %>%
  data.frame() %>%
  rownames_to_column(var="Gene")



#
#### DE ANALYSIS (quality control and plots) ------------------------------------------------------------------------------------------------------
#
## PCA
# Uses VST-normalized counts of the co-cultures to calculate PCA - quality control step
vst_mat.pca <- rda(t(vst_mat))
vst_mat_vscores <- data.frame(vst_mat.pca$CA$u)
vst_mat_vscores$Condition <- samples$condition
constrained_eigPCA<-eigenvals(vst_mat.pca)/ sum(eigenvals(vst_mat.pca))*100
# Plots calculated PCA
ggplot(vst_mat_vscores, aes(x = PC1, y = PC2, col = Condition, shape = Condition)) + 
  geom_point(size=5.5) + scale_color_manual(values = colors) + scale_shape_manual(values = shape) + 
  geom_text(data = vst_mat_vscores, aes(x = PC1, y = PC2, label = as.numeric(gsub("\\D", "", rownames(vst_mat_vscores)))), size=2.5, color="black") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0), legend.position = "bottom")

## VENN DIAGRAM NUMBER 1
# Creates a new table from the Squeezemeta output which contains only the transcripts included in the DESeq2 analysis
all_bacterial_data_FILTERED <- all_bacterial_data[rownames(all_bacterial_data) %in% rownames(vst_count),]
# Extracts columns with the experiment TPM counts of the filtered transcripts and also their KEGG ids
bacterial_tpm <- all_bacterial_data_FILTERED[, c(21:25,27:31, 9)]
colnames(bacterial_tpm) <- gsub("TPM.", "", colnames(bacterial_tpm))
bacterial_tpm_KEGG <- bacterial_tpm[bacterial_tpm$KEGG.ID != "",]
bacterial_tpm_KEGG <- bacterial_tpm_KEGG[rowSums(bacterial_tpm_KEGG[,1:10]) > 0,]
# Create a Venn diagram visualizing the similarities between the total functional arsenal of both microbial communities
venn.diagram(
  x = list(unique(bacterial_tpm_KEGG$KEGG.ID[rowSums(bacterial_tpm_KEGG[,6:10]) > 0]), unique(bacterial_tpm_KEGG$KEGG.ID[rowSums(bacterial_tpm_KEGG[,1:5]) > 0])),
  category.names = c("DVM", "DPS"),
  main.fontfamily = "calibri",
  main = "Total expressed KEGGs",
  inverted=TRUE, main.cex = 2, cat.cex = 2, cat.pos = c(130,210),
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = as.character(c("steelblue1", "#00BA38")),
  disable.logging = TRUE, 
  # Numbers
  cex = 3,
  fontface = "bold",
  fontfamily = "calibri",
  filename = 'venn_diagramm_total_bac_KEGG_new.png',
  output=TRUE
)


## VENN DIAGRAM NUMBER 2
# Extract the results table, positive lfc correspond to transcripts upregulated in DVM, negative lfc to transcripts upregulated in DPS
sig_DPS_DVM <- res_table_DPS_DVM_shrink_tb[res_table_DPS_DVM_shrink_tb$padj < padj.cutoff & abs(res_table_DPS_DVM_shrink_tb$log2FoldChange) > lfc.cutoff , ]
# Add Phylum annotation to the DE transcripts
sig_DPS_DVM$Phylum <- unlist(lapply(sig_DPS_DVM$Gene, function(x) all_bacterial_data$Tax[rownames(all_bacterial_data) == x][1]))
sig_DPS_DVM$Phylum <- as.character(sig_DPS_DVM$Phylum)
sig_DPS_DVM$Phylum <- sapply(strsplit(sapply(strsplit(sig_DPS_DVM$Phylum, "p_"), `[`, 2), ";"), `[`, 1)
sig_DPS_DVM$Phylum[is.na(sig_DPS_DVM$Phylum)] <- "Unknown Phylum"

# Calculate the contribution of each phylum to DEs
round(table(sig_DPS_DVM$Phylum[sig_DPS_DVM$log2FoldChange > 0.5])*100/sum(table(sig_DPS_DVM$Phylum[sig_DPS_DVM$log2FoldChange > 0.5])))
round(table(sig_DPS_DVM$Phylum[sig_DPS_DVM$log2FoldChange < -0.5])*100/sum(table(sig_DPS_DVM$Phylum[sig_DPS_DVM$log2FoldChange < -0.5])))

# Simple calculation of the DE transcripts upregulated in DVM and DPS
nrow(sig_DPS_DVM[sig_DPS_DVM$log2FoldChange > 0,])
nrow(sig_DPS_DVM[sig_DPS_DVM$log2FoldChange < 0,])

# Add KEGG annotation to the DE transcripts
sig_DPS_DVM_KEGG <- sig_DPS_DVM
sig_DPS_DVM_KEGG$KEGG <- unlist(lapply(sig_DPS_DVM_KEGG$Gene, function(x) all_bacterial_data$KEGG.ID[rownames(all_bacterial_data) == x][1]))
sig_DPS_DVM_KEGG <- sig_DPS_DVM_KEGG[sig_DPS_DVM_KEGG$KEGG != "", ]
table(sig_DPS_DVM_KEGG$Phylum, sig_DPS_DVM_KEGG$type)
sig_DPS_DVM_KEGG$KEGG_type <- ifelse(sig_DPS_DVM_KEGG$KEGG %in% shared_KEGGs_DE, "unique", "shared")
table(sig_DPS_DVM_KEGG$KEGG_type, sig_DPS_DVM_KEGG$type)

## how many DEs have a KEGG ID
nrow(sig_DPS_DVM); nrow(sig_DPS_DVM_KEGG); nrow(sig_DPS_DVM) - nrow(sig_DPS_DVM_KEGG)
round(nrow(sig_DPS_DVM_KEGG)*100/nrow(sig_DPS_DVM), 1)

sig_DPS_DVM_KEGG$KEGG_path <- unlist(lapply(sig_DPS_DVM_KEGG$Gene, function(x) all_bacterial_data_FILTERED$KEGGPATH[rownames(all_bacterial_data_FILTERED) == x][1]))
sig_DPS_DVM_KEGG$type <- ifelse(sig_DPS_DVM_KEGG$log2FoldChange > 0, "DVM", "DPS")
sig_DPS_DVM_KEGG$name <- unlist(lapply(sig_DPS_DVM_KEGG$Gene, function(x) all_bacterial_data_FILTERED$KEGGFUN[rownames(all_bacterial_data_FILTERED) == x][1]))
# Get the shared KEGG IDs between "DPS" and "DVM" upregulated transcripts
shared_KEGGs_DE <- Reduce(intersect, list(unique(sig_DPS_DVM_KEGG$KEGG[sig_DPS_DVM_KEGG$type=="DPS"]), unique(sig_DPS_DVM_KEGG$KEGG[sig_DPS_DVM_KEGG$type=="DVM"])))
# Get the unique KEGG IDs for "DPS"
unique_DPS_KEGGs_DE <- Reduce(setdiff,  list(unique(sig_DPS_DVM_KEGG$KEGG[sig_DPS_DVM_KEGG$type=="DPS"]), unique(sig_DPS_DVM_KEGG$KEGG[sig_DPS_DVM_KEGG$type=="DVM"])))
# Get the unique KEGG IDs for "DVM"
unique_DVM_KEGGs_DE <- Reduce(setdiff,  list(unique(sig_DPS_DVM_KEGG$KEGG[sig_DPS_DVM_KEGG$type=="DVM"]), unique(sig_DPS_DVM_KEGG$KEGG[sig_DPS_DVM_KEGG$type=="DPS"])))
table(sig_DPS_DVM_KEGG$type)
sig_DPS_DVM_KEGG$gene_name <- unlist(lapply(sig_DPS_DVM_KEGG$Gene, function(x) all_bacterial_data_FILTERED$Gene.name[rownames(all_bacterial_data_FILTERED) == x][1]))
sig_DPS_DVM_KEGG$gene_name <- as.character(sig_DPS_DVM_KEGG$gene_name)

# Create a Venn diagram visualizing the similarities between the functional annotation of the DE transcripts
venn.diagram(
  x = list(unique(sig_DPS_DVM_KEGG$KEGG[sig_DPS_DVM_KEGG$type=="DPS"]), unique(sig_DPS_DVM_KEGG$KEGG[sig_DPS_DVM_KEGG$type=="DVM"])),
  #x = list(unique(sig_DPS_DVM_KEGG$KEGG[sig_DPS_DVM_KEGG$type=="DVM"]), unique(sig_DPS_DVM_KEGG$KEGG[sig_DPS_DVM_KEGG$type=="DPS"])),
  category.names = c("DPS", "DVM"),
  ext.line.lwd = 0, ext.pos=c(1,1),
  main.fontfamily = "calibri",
  main = "KEGGs of DE Transcripts",
  inverted=TRUE, main.cex = 2, cat.cex = 2, cat.pos = c(40,320), cat.dist=0.07,
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = as.character(c("#00BA38", "steelblue1")),
  disable.logging = TRUE, 
  # Numbers
  cex = 3,
  fontface = "bold",
  fontfamily = "calibri",
  filename = 'venn_diagramm_DVM_and_DPS_bac_KEGG.png',
  output=TRUE
)

## HEATMAPS
# Calculating Euclidean distance matrix from transposed vst_mat
dist_vst_mat <- dist(t(vst_mat), method = "euclidean")

# Converting the distance matrix to a matrix format for heatmap
sampleDistMatrix <- as.matrix(dist_vst_mat)

# Defining a color palette for the heatmap
heatmap_colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)

# Creating a heatmap of the distance matrix using pheatmap
eucl_dist_heatmap <- pheatmap(sampleDistMatrix,
                              clustering_distance_rows = dist_vst_mat,
                              clustering_distance_cols = dist_vst_mat,
                              color = heatmap_colors)

# Converting the distance matrix to a data frame for further analysis
dist_vst_mat <- as.data.frame(as.matrix(dist_vst_mat))

# Calculating the mean and standard deviation of distances within the first 5 samples, excluding zero distances
mean(unlist(dist_vst_mat[1:5,1:5])[unlist(dist_vst_mat[1:5,1:5]) > 0 ])
sd(unlist(dist_vst_mat[1:5,1:5])[unlist(dist_vst_mat[1:5,1:5]) > 0 ])

# Calculating the mean and standard deviation of distances within the next 5 samples, excluding zero distances
mean(unlist(dist_vst_mat[6:10,6:10])[unlist(dist_vst_mat[6:10,6:10]) > 0 ])
sd(unlist(dist_vst_mat[6:10,6:10])[unlist(dist_vst_mat[6:10,6:10]) > 0 ])

# Conducting a Wilcoxon test to compare the distance distributions between the first 5 and the next 5 samples
wilcox.test(unlist(dist_vst_mat[1:5,1:5])[unlist(dist_vst_mat[1:5,1:5]) > 0 ],
                              unlist(dist_vst_mat[6:10,6:10])[unlist(dist_vst_mat[6:10,6:10]) > 0 ],
                              p.adjust.methods="BH", exact=FALSE)$p.value

# Performing a permutational multivariate analysis of variance using Euclidean distance
vegan::adonis2(as.matrix(t(vst_mat)) ~ condition, data=samples, permutations=999, method="euclidean")

### Compute pairwise correlation values
vst_cor <- cor(vst_mat) 
### Plot heatmap
pheatmap(vst_cor, annotation = samples)


## VOLCANO PLOT 
res_table_DPS_DVM_shrink_tb <- as.data.frame(res_table_DPS_DVM_shrink_tb)
# Removes all rows from "res_table_DPS_DVM_shrink_tb" where the "padj" column has a missing value
res_table_DPS_DVM_shrink_tb <- res_table_DPS_DVM_shrink_tb[ ! is.na(res_table_DPS_DVM_shrink_tb$padj) , ]
# Adds a new column named "Phylum" and assigns to each row the corresponding taxonomic information
res_table_DPS_DVM_shrink_tb$Phylum <- unlist(lapply(res_table_DPS_DVM_shrink_tb$Gene, function(x) all_bacterial_data$Tax[rownames(all_bacterial_data) == x][1]))
res_table_DPS_DVM_shrink_tb$Phylum <- as.character(res_table_DPS_DVM_shrink_tb$Phylum)
res_table_DPS_DVM_shrink_tb$Phylum <- sapply(strsplit(sapply(strsplit(res_table_DPS_DVM_shrink_tb$Phylum, "p_"), `[`, 2), ";"), `[`, 1)
# Replaces any missing values in the "Phylum" column with the string "Unknown Phylum"
res_table_DPS_DVM_shrink_tb$Phylum[is.na(res_table_DPS_DVM_shrink_tb$Phylum)] <- "Unknown Phylum"
# This line assigns "non_DE" to "Phylum" for rows where either absolute log2FoldChange is less than 0.5 or "padj" is greater than 0.05
res_table_DPS_DVM_shrink_tb$Phylum[abs(res_table_DPS_DVM_shrink_tb$log2FoldChange) < lfc.cutoff | res_table_DPS_DVM_shrink_tb$padj > padj.cutoff ] <- "non_DE"
res_table_DPS_DVM_shrink_tb$col <- ifelse(res_table_DPS_DVM_shrink_tb$Phylum %in% c("Proteobacteria", "Planctomycetes", "Bacteroidetes", "Firmicutes", "Actinobacteria", "Cyanobacteria", "non_DE"), res_table_DPS_DVM_shrink_tb$Phylum, "Other")
# Renames according to the SILVA standard
res_table_DPS_DVM_shrink_tb$col[res_table_DPS_DVM_shrink_tb$col=="Planctomycetes"] <- "Planctomycetota"
res_table_DPS_DVM_shrink_tb$col[res_table_DPS_DVM_shrink_tb$col=="Bacteroidetes"] <- "Bacteroidota"
res_table_DPS_DVM_shrink_tb$col[res_table_DPS_DVM_shrink_tb$col=="Actinobacteria"] <- "Actinobacteriota"
# Consistent colors
RNA_colors_bacteria_class <- c("black", "grey", "navy", "#FF6347","#2ed507", "#CDAD00", "yellow1", "#BFEFFF")
names(RNA_colors_bacteria_class) <- unique(res_table_DPS_DVM_shrink_tb$col)

## Volcano plot
ggplot(res_table_DPS_DVM_shrink_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color=col), size=.5)  +  
  scale_color_manual(values = RNA_colors_bacteria_class) +
  geom_hline(yintercept = -log10(padj.cutoff)) +
  geom_vline(xintercept = lfc.cutoff) +
  geom_vline(xintercept = (-lfc.cutoff)) +
  geom_text(y=2,x=-20, label="p-value < 0.05", check_overlap = TRUE) +
  geom_text(y=2,x=20, label="p-value < 0.05", check_overlap = TRUE) +
  geom_text(y=15,x=-5, label="l2fc < - 0.5", check_overlap = TRUE) +
  geom_text(y=15,x=5, label="l2fc > 0.5", check_overlap = TRUE) +
  geom_text(y=17,x=-15, label="DPS", check_overlap = TRUE, color="#00BA38", size=8) +
  geom_text(y=17,x=15, label="DVM", check_overlap = TRUE, color="#619CFF", size=8)+ theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))+
  guides(color = guide_legend(override.aes = list(size = 4, shape = 15))) + labs(color = "Phylum") +
#  xlim(-20, 30)  +  ylim(-2, 20)  +
  xlab("log2 fold change") +  
  ylab("-log10 adjusted p-value")#




#
######## ENRICHMENT ANALYSIS ------------------------------------------------------------------------------------------------------
#
# Find enriched KEGG ids from the upregulated in DVM tranacripts against the custom background
enrich_DVM <- enricher(as.character(unique(sig_DPS_DVM$Gene[sig_DPS_DVM$log2FoldChange>0])), TERM2GENE=vst_mat_annotated[,c(11:12)])
head(enrich_DVM)
enrich_DVM <- as.data.frame(enrich_DVM)

# Annotate the enriched KEGG ids
enrich_DVM$KEGG_path <- unlist(lapply(enrich_DVM$ID, function(x) all_bacterial_data_FILTERED$KEGGPATH[all_bacterial_data_FILTERED$KEGG.ID == x][1]))
enrich_DVM$KEGG_name <- unlist(lapply(enrich_DVM$ID, function(x) all_bacterial_data_FILTERED$KEGGFUN[all_bacterial_data_FILTERED$KEGG.ID == x][1]))
enrich_DVM$type <- "DVM"

enrich_DPS <- enricher(as.character(unique(sig_DPS_DVM$Gene[sig_DPS_DVM$log2FoldChange<0])), TERM2GENE=vst_mat_annotated[,c(11:12)])
head(enrich_DPS)
enrich_DPS <- as.data.frame(enrich_DPS)


#
#### BACTERIAL SECRETION (interpreted from the SignalP and NetGPI outputs) ------------------------------------------------------------------------------------------------------
#
# Define a vector of KEGG IDs related to secretion pathways
secretion_KEGGs <- c("K02453", "K02454", "K02455", "K02456", "K02457", "K03070", "K03071", "K03072", "K03073", "K03074", "K03075", "K03076", "K03106", "K03110", "K03116", "K03117", "K03118", "K03210", "K03217", "K11892", "K11903", "K11904", "K11912", "K12257", "K12340")

# Filter sig_DPS_DVM_KEGG for positive log2FoldChange and KEGG IDs in secretion_KEGGs
secretion_DVM <- sig_DPS_DVM_KEGG[sig_DPS_DVM_KEGG$log2FoldChange > 0 & sig_DPS_DVM_KEGG$KEGG %in% secretion_KEGGs, ]

# Calculate and round the percentage of each Phylum within secretion_DVM
round(table(secretion_DVM$Phylum) * 100 / sum(table(secretion_DVM$Phylum)), 2)

# Define KEGG IDs for type 2 secretion systems
type2_sec <- c("K02453", "K02455", "K02456", "K02457", "K02454", "K03072", "K12257", "K03073", "K03074", "K03075", "K03076", "K03210", "K03217", "K03070", "K03110", "K03071", "K03106", "K03116", "K03117","K03118")

# Filter for type 2 secretion systems in sig_DPS_DVM_KEGG
secretion_type2_DVM <- sig_DPS_DVM_KEGG[sig_DPS_DVM_KEGG$log2FoldChange > 0 & sig_DPS_DVM_KEGG$KEGG %in% type2_sec, ]
round(table(secretion_type2_DVM$Phylum) * 100 / sum(table(secretion_type2_DVM$Phylum)), 2)

# Define KEGG IDs for type 6 secretion systems
type6_sec <- c("K11904", "K11903", "K11892", "K11912")

# Filter for type 6 secretion systems in sig_DPS_DVM_KEGG
secretion_type6_DVM <- sig_DPS_DVM_KEGG[sig_DPS_DVM_KEGG$log2FoldChange > 0 & sig_DPS_DVM_KEGG$KEGG %in% type6_sec, ]
round(table(secretion_type6_DVM$Phylum) * 100 / sum(table(secretion_type6_DVM$Phylum)), 2)

# Reading SignalP prediction data for DPS
DPS_signalP <- read.table("DPS_output_protein_type.txt", sep="\t", header = TRUE)
table(DPS_signalP$Prediction)

# Convert IDs to character for processing
DPS_signalP$ID <- as.character(DPS_signalP$ID)

# Reading NetGPI data for DPI
Net_GPI_DPS <- read.table("NetGPI_DPS_with_SP.txt", sep="\t", header=TRUE)
Net_GPI_DPS$ID <- as.character(Net_GPI_DPS$ID); Net_GPI_DPS$Pred..GPI.Anchored <- as.character(Net_GPI_DPS$Pred..GPI.Anchored)

# Map KEGG IDs and GPI Anchored predictions to SignalP DPS data
DPS_signalP$KEGG <- unlist(lapply(DPS_signalP$ID, function(x) all_bacterial_data_FILTERED$KEGG.ID[rownames(all_bacterial_data_FILTERED) == x][1]))
DPS_signalP$Anch <- unlist(lapply(DPS_signalP$ID, function(x) Net_GPI_DPS$Pred..GPI.Anchored[Net_GPI_DPS$ID == x][1]))

# Identify differentially expressed (DE) transcripts in SignalP data that are excreted in DPS
DE_enriched_Signal_EXCRETED_DPS <- enrich_DPS$ID[enrich_DPS$ID %in% DPS_signalP$KEGG[DPS_signalP$KEGG != "" & DPS_signalP$Prediction != "OTHER"] ]
unique(all_bacterial_data_FILTERED[all_bacterial_data_FILTERED$KEGG.ID %in% DE_enriched_Signal_EXCRETED_DPS, c(9, 10)])

# Prepare an empty dataframe for DVM SignalP data
DVM_signalP <- data_frame(ID=character(), Prediction=character(), SP.Sec.SPI.=numeric(), TAT.Tat.SPI.=numeric(), LIPO.Sec.SPII.=numeric(), OTHER=numeric(), CS.Position=character())

# List files in the specified directory and remove the first entry (usually '.' or '..')
files <- list.files("./signalP/")[-1]

# Loop through files to read and combine SignalP data for DVM
for(i in 1:4){
  temp <- read.table(paste0("./signalP/", files[i]), sep="\t", header = FALSE)
  DVM_signalP <- rbind(DVM_signalP, temp)
}

# Set column names for DVM_signalP dataframe
colnames(DVM_signalP) <- colnames(DPS_signalP)[1:7]
table(DVM_signalP$Prediction)

# Reading NetGPI data for DVM
Net_GPI_DVM <- read.table("NetGPI_DVM_with_SP.txt", sep="\t", header=TRUE)
Net_GPI_DVM$ID <- as.character(Net_GPI_DVM$ID)
Net_GPI_DVM$Pred..GPI.Anchored <- as.character(Net_GPI_DVM$Pred..GPI.Anchored)

# Identify unique KEGG IDs associated with genes in DVM_signalP predicted as 'SP(Sec/SPI)' and having positive log2FoldChange
unique(sig_DPS_DVM_KEGG$KEGG[sig_DPS_DVM_KEGG$Gene %in% DVM_signalP$ID[DVM_signalP$Prediction == "SP(Sec/SPI)"] & sig_DPS_DVM_KEGG$log2FoldChange > 0])

# Identify unique KEGG IDs associated with genes in DVM_signalP predicted as 'TAT(Tat/SPI)' and having positive log2FoldChange
unique(sig_DPS_DVM_KEGG$KEGG[sig_DPS_DVM_KEGG$Gene %in% DVM_signalP$ID[DVM_signalP$Prediction == "TAT(Tat/SPI)"] & sig_DPS_DVM_KEGG$log2FoldChange > 0])

# Map KEGG IDs to DVM_signalP data based on the 'ID' column
DVM_signalP$KEGG <- unlist(lapply(DVM_signalP$ID, function(x) all_bacterial_data_FILTERED$KEGG.ID[rownames(all_bacterial_data_FILTERED) == x][1]))

# Map GPI Anchored predictions to DVM_signalP data based on the 'ID' column
DVM_signalP$Anch <- unlist(lapply(DVM_signalP$ID, function(x) Net_GPI_DVM$Pred..GPI.Anchored[Net_GPI_DVM$ID == x][1]))

# Identify DE enriched transcripts in DVM_signalP that are excreted (excluding 'OTHER' predictions)
DE_enriched_Signal_EXCRETED_DVM <- enrich_DVM$ID[enrich_DVM$ID %in% DVM_signalP$KEGG[DVM_signalP$KEGG != "" & DVM_signalP$Prediction != "OTHER"]]

# Filter sig_DPS_DVM_KEGG for KEGG IDs in DE_enriched_Signal_EXCRETED_DVM and positive log2FoldChange
sig_DVM_excreted <- sig_DPS_DVM_KEGG[sig_DPS_DVM_KEGG$KEGG %in% DE_enriched_Signal_EXCRETED_DVM & sig_DPS_DVM_KEGG$log2FoldChange > 0 , ]

# Generate a table of the distribution of 'Phylum' in excreted transcripts
table(sig_DVM_excreted$Phylum)

# Identify unique pairs of Phylum and corresponding KEGG function for excreted transcripts
unique(all_bacterial_data_FILTERED[all_bacterial_data_FILTERED$KEGG.ID %in% DE_enriched_Signal_EXCRETED_DVM, c(9, 10)])

# Identify DE enrichedtranscripts in DVM_signalP that are GPI-Anchored and belong to either 'SP(Sec/SPI)' or 'TAT(Tat/SPI)' predictions
DE_enriched_Signal_ANCHORED_DVM <- enrich_DVM$ID[enrich_DVM$ID %in% DVM_signalP$KEGG[DVM_signalP$Prediction %in% c("SP(Sec/SPI)", "TAT(Tat/SPI)")] & enrich_DVM$ID %in% DVM_signalP$KEGG[DVM_signalP$Anch == "GPI-Anchored"]]

# Identify DE enriched transcripts in DVM_signalP predicted as 'LIPO(Sec/SPII)'
DE_enriched_LIPO_DVM <- enrich_DVM$ID[enrich_DVM$ID %in% DVM_signalP$KEGG[DVM_signalP$Prediction == "LIPO(Sec/SPII)"]]

# Identify unique KEGG functions associated with DE LIPO transcripts
unique(all_bacterial_data_FILTERED$KEGGFUN[all_bacterial_data_FILTERED$KEGG.ID %in% DE_enriched_LIPO_DVM])


#
######## THE ANALYSIS OF ALL DE FUNCTIONS ------------------------------------------------------------------------------------------------------
#
# Duplicate the 'sig_DPS_DVM_KEGG' dataframe for figure creation purposes
sig_DPS_DVM_KEGG_fig <- sig_DPS_DVM_KEGG

# Map COG (Clusters of Orthologous Groups) pathways to genes in 'sig_DPS_DVM_KEGG_fig'
# This is done by applying a function to each 'Gene' in the dataframe that finds the corresponding COG pathway in 'all_bacterial_data'
sig_DPS_DVM_KEGG_fig$COG <- unlist(lapply(sig_DPS_DVM_KEGG_fig$Gene, function(x) all_bacterial_data$COGPATH[rownames(all_bacterial_data) == x][1]))

# Check for genes without an assigned COG pathway (empty strings)
table(sig_DPS_DVM_KEGG_fig$COG == "")

# Convert the 'COG' column to character type
sig_DPS_DVM_KEGG_fig$COG <- as.character(sig_DPS_DVM_KEGG_fig$COG)

# Remove rows where the COG pathway is not specified (empty string)
sig_DPS_DVM_KEGG_fig <- sig_DPS_DVM_KEGG_fig[sig_DPS_DVM_KEGG_fig$COG != "", ]

# Filter out rows where COG contains a pipe (|) character, indicating multiple annotations
sig_DPS_DVM_KEGG_fig <- sig_DPS_DVM_KEGG_fig[ ! grepl("\\|", sig_DPS_DVM_KEGG_fig$COG), ]

# Exclude rows with 'Function unknown' in the COG column
sig_DPS_DVM_KEGG_fig <- sig_DPS_DVM_KEGG_fig[ ! grepl("Function unknown", sig_DPS_DVM_KEGG_fig$COG), ]

# Filter out general function predictions from the COG column
sig_DPS_DVM_KEGG_fig <- sig_DPS_DVM_KEGG_fig[sig_DPS_DVM_KEGG_fig$COG != "General function prediction only", ]

# Display frequency tables of COG pathways for 'DVM' and 'DPS' types, sorted in ascending order
table(sig_DPS_DVM_KEGG_fig$COG[sig_DPS_DVM_KEGG_fig$type  == "DVM"]) %>% sort(decreasing = FALSE)
table(sig_DPS_DVM_KEGG_fig$COG[sig_DPS_DVM_KEGG_fig$type  == "DPS"]) %>% sort(decreasing = FALSE)

# Convert the 'COG' column to a factor with levels ordered based on frequency in ascending order
sig_DPS_DVM_KEGG_fig$COG <- factor(sig_DPS_DVM_KEGG_fig$COG, levels=names(table(sig_DPS_DVM_KEGG_fig$COG) %>% sort(decreasing = FALSE)))

# Create a ggplot of DE transcripts categorized by COG function
# Points are colored by type ('DPS' or 'DVM') and positioned by log2 fold change
COG_DE_bacteria <- ggplot(sig_DPS_DVM_KEGG_fig, aes(x = log2FoldChange, y = COG, color=type)) +
  geom_point()  + 
  facet_nested(.~type, scales = 'free', space = "free") + 
  xlab("log2 fold change") + 
  ylab("") + 
  labs(fill = "") + theme_bw() + 
  guides(fill=guide_legend(ncol=5)) + theme(legend.position="none", 
                                            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1), 
                                            panel.grid = element_line(color = "white"), strip.background =element_rect(fill="white"),
                                            strip.text.y = element_text(angle = 0))  +
  scale_color_manual(values=setNames(c("#00BA38", "steelblue1"), c("DPS", "DVM"))) +
  ggtitle("Microbiome DE transcripts - COG function")

# Create contingency tables for 'Day' and 'Night' types, comparing COG pathways and Kingdom classification
table(sig_DPS_DVM_KEGG_fig$COG[sig_DPS_DVM_KEGG_fig$type == "Day"], sig_DPS_DVM_KEGG_fig$Kingdom[sig_DPS_DVM_KEGG_fig$type == "Day"])
table(sig_DPS_DVM_KEGG_fig$COG[sig_DPS_DVM_KEGG_fig$type == "Night"], sig_DPS_DVM_KEGG_fig$Kingdom[sig_DPS_DVM_KEGG_fig$type == "Night"])


#
######## MARKER FUNCTIONS ------------------------------------------------------------------------------------------------------
#
# Upload the list of used markers
marker_genes <- read.table("marker_gene_list.csv", sep=",", header=TRUE)
rownames(marker_genes) <- marker_genes$X; marker_genes$X <- NULL
marker_genes$KEGG <- as.character(marker_genes$KEGG)
marker_genes$Step <- as.character(marker_genes$Step)
# Edit the titles to show the number of KEGGs in each tested functional category
for(i in 1:length(unique(marker_genes$Step))){
  str <- paste0(unique(marker_genes$Step)[i], " (n = ", length(unique(marker_genes$KEGG[marker_genes$Step == unique(marker_genes$Step)[i]])), ")")
  marker_genes$Step[marker_genes$Step ==unique(marker_genes$Step)[i]] <- rep(str, length(marker_genes$Step[marker_genes$Step == unique(marker_genes$Step)[i]]))
}
## find out how many marker genes are in the data
round(length(marker_genes$KEGG[marker_genes$KEGG %in% bacterial_tpm$KEGG.ID])*100/length(marker_genes$KEGG))
round(length(marker_genes$KEGG[  !marker_genes$KEGG %in% bacterial_tpm$KEGG.ID])*100/length(marker_genes$KEGG))
unique(marker_genes$Step[  !marker_genes$KEGG %in% bacterial_tpm$KEGG.ID])

# Subset TPM counts of the transcripts coding for the marker KEGG ids
marker_genes_tpm <- bacterial_tpm_KEGG[bacterial_tpm_KEGG$KEGG.ID %in% marker_genes$KEGG,]
# Add annotation to the marker transcripts
marker_genes_tpm$Cycle <- unlist(lapply(marker_genes_tpm$KEGG.ID, function(x) marker_genes$Cycle[marker_genes$KEGG == x]))
marker_genes_tpm$Step <- unlist(lapply(marker_genes_tpm$KEGG.ID, function(x) marker_genes$Step[marker_genes$KEGG == x]))

# Assessing the homogeneity of multivariate dispersions for different conditions
betadisper(vegdist(as.matrix(t(marker_genes_tpm[,c(1:10)]))), samples$condition)

# Conducting a multivariate analysis of variance using Bray-Curtis dissimilarity on a transposed subset of marker_genes_tpm
vegan::adonis2(as.matrix(t(marker_genes_tpm[,c(1:10)])) ~ condition, data=samples, permutations=999, method="bray")

# Performing the same multivariate analysis of variance as above but using Euclidean distance
vegan::adonis2(as.matrix(t(marker_genes_tpm[,c(1:10)])) ~ condition, data=samples, permutations=999, method="euclidean")

# Subsetting the vst_mat to include only rows that are in marker_genes_tpm
vst_mat_marker <- vst_mat[rownames(vst_mat) %in% rownames(marker_genes_tpm), ]

# Conducting a multivariate analysis of variance on the variance-stabilized transformed data using Euclidean distance
vegan::adonis2(as.matrix(t(vst_mat_marker)) ~ condition, data=samples, permutations=999, method="euclidean")

# Add taxa of the marker genes
marker_genes_tpm$Phylum <- unlist(lapply(rownames(marker_genes_tpm), function(x) all_bacterial_data_FILTERED$Tax[rownames(all_bacterial_data_FILTERED) == x]))
marker_genes_tpm$Phylum <- sapply(strsplit(sapply(strsplit(marker_genes_tpm$Phylum, "p_"), `[`, 2), ";"), `[`, 1)

# Adding transcript names as a new column in marker_genes_tpm dataframe
marker_genes_tpm$transcript <- rownames(marker_genes_tpm)

# Calculating distances between columns (transcripts) and rounding off the diagonal values of the distance matrix
distances <- dist(t(marker_genes_tpm[,1:10]))
round(diag(as.matrix(distances)[,-1]),1)

# Performing hierarchical clustering on the distance matrix and manually adjusting the order of clusters
colclus <- hclust(distances)
colclus$order <- c(4,5,9,6,3,7,8,10,1,2)

# Reshaping the marker_genes_tpm data to a long format for easier analysis and plotting
marker_genes_tpm_long <- reshape2::melt(marker_genes_tpm, id.vars=c("Cycle", "Step", "KEGG.ID", "Phylum", "transcript"))

# Assigning types to the data based on the presence of 'DPS' or 'DVM' in the variable names
marker_genes_tpm_long$type[grepl("DPS", marker_genes_tpm_long$variable)] <- "DPS"
marker_genes_tpm_long$type[grepl("DVM", marker_genes_tpm_long$variable)] <- "DVM"

# Creating contingency tables to compare the transcripts in marker_genes_tpm_long with those in sig_DPS_DVM
table(unique(marker_genes_tpm_long$transcript) %in% sig_DPS_DVM$Gene)
table(unique(marker_genes_tpm_long$transcript) %in% sig_DPS_DVM$Gene[sig_DPS_DVM$log2FoldChange < 0]) # For DE DPS functional genes
unique(marker_genes_tpm_long$KEGG.ID[marker_genes_tpm_long$transcript %in% sig_DPS_DVM$Gene[sig_DPS_DVM$log2FoldChange < 0]])
unique(marker_genes_tpm_long$Phylum[marker_genes_tpm_long$transcript %in% sig_DPS_DVM$Gene[sig_DPS_DVM$log2FoldChange < 0]])

# Repeating the process for DE DVM functional genes
table(unique(marker_genes_tpm_long$transcript) %in% sig_DPS_DVM$Gene[sig_DPS_DVM$log2FoldChange > 0])
unique(marker_genes_tpm_long$KEGG.ID[marker_genes_tpm_long$transcript %in% sig_DPS_DVM$Gene[sig_DPS_DVM$log2FoldChange > 0]])
unique(marker_genes_tpm_long$Phylum[marker_genes_tpm_long$transcript %in% sig_DPS_DVM$Gene[sig_DPS_DVM$log2FoldChange > 0]])

# Checking if all KEGG IDs coded by DE DPS are also included in DE DVM KEGG IDs to infer functional redundancy
all(unique(marker_genes_tpm_long$KEGG.ID[marker_genes_tpm_long$transcript %in% sig_DPS_DVM$Gene[sig_DPS_DVM$log2FoldChange < 0]]) %in% unique(marker_genes_tpm_long$KEGG.ID[marker_genes_tpm_long$transcript %in% sig_DPS_DVM$Gene[sig_DPS_DVM$log2FoldChange > 0]]))

# Aggregating the data based on certain conditions and computing the mean
marker_genes_tpm_long_aggregated <- aggregate(marker_genes_tpm_long$value, by=list(marker_genes_tpm_long$Cycle, marker_genes_tpm_long$Step, marker_genes_tpm_long$variable, marker_genes_tpm_long$type), mean)
colnames(marker_genes_tpm_long_aggregated) <- colnames(marker_genes_tpm_long)[c(1,2,6, 8, 7)]

# Converting 'Cycle' to title case and 'Step' to character type
marker_genes_tpm_long_aggregated$Cycle <- str_to_title(marker_genes_tpm_long_aggregated$Cycle)
marker_genes_tpm_long_aggregated$Step <- as.character(marker_genes_tpm_long_aggregated$Step)

# Applying log10 transformation to the aggregated value column
marker_genes_tpm_long_aggregated$value <- log10(marker_genes_tpm_long_aggregated$value)

# Setting the factor levels for the 'variable' column based on a specific order
marker_genes_tpm_long_aggregated$variable <- factor(marker_genes_tpm_long_aggregated$variable, levels=unique(marker_genes_tpm_long_aggregated$variable)[c(4,5,9,6,3,7,8,10,1,2)])

# Plotting a dendrogram using ggplot2 along with customized theme settings
ggplot() + geom_dendro(colclus) +  theme_hm()  + theme( axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

# Assigning colors to the replicates and then generating a heat map using ggplot2
rep_colors <- rep("col", 10)
names(rep_colors) <- unique(marker_genes_tpm_long_aggregated$variable)[c(4,5,9,6,3,7,8,10,1,2)]
rep_colors <- ifelse(grepl("DPS", names(rep_colors)),  "#00BA38", "#619CFF")
names(rep_colors) <- unique(marker_genes_tpm_long_aggregated$variable)[c(4,5,9,6,3,7,8,10,1,2)]

ggplot() +
  geom_tile(data=marker_genes_tpm_long_aggregated, aes(x=variable, y=Step, fill=value)) +  theme_hm() +  
  labs(y="Significant DE Transcript", fill = "log-transformed mean TPM")  +  
  scale_fill_gradient(na.value = "white") + 
  facet_grid(Cycle~., scales="free", space = "free") + theme_bw() +
  theme(legend.position = "bottom", strip.text.y = element_text(angle = 0), axis.text.x = element_text(color = rep_colors, size = 10, angle = 90, vjust = 0.5, hjust=1), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid = element_line(color = "white"), strip.background =element_rect(fill="white"))

# Further dendrogram plotting for a specific step "Biotin biosynthesis (n = 1)"
ggplot(marker_genes_tpm_long_aggregated[marker_genes_tpm_long_aggregated$Step=="Biotin biosynthesis (n = 1)", ], aes(variable, value)) + geom_dendro(colclus)