
library("phyloseq") # version 1.30.0
library("ggh4x") # version 0.2.3
library("ggpubr") # version 0.5.0
library("vegan") # version 2.6.4
library("ggvegan") # version 0.1.0
library("plotly") # version 4.10.1
library("stringr") # version 1.5.0
library("reshape2") # version 1.4.4
library("dplyr") # version 1.0.10
library("eoffice") # version 0.2.2
library("seqinr") # version 4.2.30
library("forcats") # version 1.0.0
library("tibble") # version 3.2.1
library("tidyr") # version 1.3.0

set.seed(1)

# RNAseq Dataset ------------------------------------------------------------------------------------

## Data Import and Initial Processing

# Import SQM-derived table representing filtered RNAseq taxa data
all_RNAseq_bacterial_data_FILTERED <- read.table(
  "all_bacterial_data_FILTERED.txt", 
  sep = "\t"
)

# Select relevant columns
all_RNAseq_bacterial_data_FILTERED <- all_RNAseq_bacterial_data_FILTERED[, c(16:19, 21:25, 27:36, 8)]

# Convert Tax column to character
all_RNAseq_bacterial_data_FILTERED$Tax <- as.character(all_RNAseq_bacterial_data_FILTERED$Tax)

# Add columns for Kingdom and Phylum, initialize with 'Unknown'
all_RNAseq_bacterial_data_FILTERED[, c("Kingdom", "Phylum")] <- "Unknown"

# Extract taxonomic information from Tax column
for(i in 1:nrow(all_RNAseq_bacterial_data_FILTERED)) {
  # Split taxonomic information
  vec <- str_split(all_RNAseq_bacterial_data_FILTERED$Tax[i], ";")[[1]]
  
  # Extract Kingdom and Phylum
  all_RNAseq_bacterial_data_FILTERED$Kingdom[i] <- str_replace(str_subset(vec, "^k_"), "k_", "")
  all_RNAseq_bacterial_data_FILTERED$Phylum[i] <- ifelse(
    any(grepl("p_", vec)),
    str_replace(str_subset(vec, "^p_"), "p_", ""),
    "Unknown"
  )
}

# Assign rownames to ID column
all_RNAseq_bacterial_data_FILTERED$ID <- rownames(all_RNAseq_bacterial_data_FILTERED)

# Data Transformation

# Melt data and rename columns
all_RNAseq_bacterial_data_FILTERED_long <- melt(
  all_RNAseq_bacterial_data_FILTERED[, -20], 
  id.vars = c("Kingdom", "Phylum", "ID")
) %>%
  mutate(value = as.numeric(value))
colnames(all_RNAseq_bacterial_data_FILTERED_long)[4] <- "Sample"

# Transform counts to relative abundance
all_RNAseq_long_rel_abund <- all_RNAseq_bacterial_data_FILTERED_long %>%
  group_by(Sample) %>%
  mutate(Abundance = value / sum(value) * 100) %>%
  ungroup()

# --- Step 3: Data Aggregation and Annotation ---

# Aggregate data by Kingdom, Phylum, and Sample
all_RNAseq_long_aggregated_rel_abund <- aggregate(
  Abundance ~ Kingdom + Phylum + Sample,
  data = all_RNAseq_long_rel_abund[,-c(3,5)],
  sum
)

# Rename low abundant phyla to 'Others'
all_RNAseq_long_aggregated_rel_abund_others <- all_RNAseq_long_aggregated_rel_abund
all_RNAseq_long_aggregated_rel_abund_others[all_RNAseq_long_aggregated_rel_abund_others$Abundance < 1, c(1,2)] <- "Others"

# Aggregate data by Kingdom, Phylum, and Sample again
all_RNAseq_long_aggregated_rel_abund_others <- aggregate(
  Abundance ~ Kingdom + Phylum + Sample,
  data = all_RNAseq_long_aggregated_rel_abund_others,
  sum
)

# Annotate the data
all_RNAseq_long_aggregated_rel_abund_others$Sample <- gsub("TPM.", "", all_RNAseq_long_aggregated_rel_abund_others$Sample)
all_RNAseq_long_aggregated_rel_abund_others$treatment <- gsub('[[:digit:]]+', '', all_RNAseq_long_aggregated_rel_abund_others$Sample)

# Subset data for specific treatments (DPS, DVM)
experiment_RNAseq <- all_RNAseq_long_aggregated_rel_abund_others[all_RNAseq_long_aggregated_rel_abund_others$treatment %in% c("DPS", "DVM"), ]

# Visualization

# Create bar plot for RNAseq data
ggplot(experiment_RNAseq, aes(x = Sample, y = Abundance, fill = fct_reorder(Phylum, Abundance, .desc = TRUE))) +
  geom_bar(stat = "identity", color = "black") +
  labs(fill = "Phylum") + theme_bw() +
  ylab("Abundance (%)") + xlab("") +
  guides(fill = guide_legend(ncol = 7)) +
  facet_nested(.~treatment, scales = "free") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.grid = element_line(color = "white"))

# Data Analysis

# Calculate mean phyla abundances

# Shannon diversity calculation
RNAseq_wide_matrix <- experiment_RNAseq[, -c(1, 5)] %>%
  pivot_wider(names_from = Sample, values_from = Abundance) %>%
  column_to_rownames(var = "Phylum") %>%
  as.matrix()
RNAseq_wide_matrix[is.na(RNAseq_wide_matrix)] <- 0
RNAseq_wide_matrix <- t(RNAseq_wide_matrix)

# Shannon diversity indices
RNA_seq_shannon_div <- diversity(RNAseq_wide_matrix, index = "shannon")
mean(RNA_seq_shannon_div[1:5]); sd(RNA_seq_shannon_div[1:5])
mean(RNA_seq_shannon_div[6:10]); sd(RNA_seq_shannon_div[6:10])

# Perform Wilcoxon rank-sum test to compare the Shannon diversity indices between two groups of samples (DPS and DVM)
wilcox.test(RNA_seq_shannon_div[1:5], RNA_seq_shannon_div[6:10], p.adjust.methods="BH", exact=FALSE)$p.value






# Amplicon Dataset ------------------------------------------------------------------------------
## RNA stands for 16S cDNA, DNA stands for 16S rRNA. RNA/DNA terms used relate to the sequenced fraction of the extracted in parallel material from the same sample

## upload dada2-generated output as phyloseq object
RNA <- readRDS("ASV_physeq_Prok.rds")
DNA <- readRDS("physeq_DNA.rds")

## initial taxonomy check - caution on Order chloroplast and Family Mitochondria
##plot_bar(DNA, fill="Phylum")
##plot_bar(RNA, fill="Order")

# Analysis of the bacterial controls ------------------------------------------------------------------------------
# Subsetting the RNA dataset to include only bacterial control samples
RNA_bac_controls <- subset_samples(RNA, sample_names(RNA) %in% c("PS4", "PS5", "PS6", "VM1", "VM2"))

# Filtering out taxa with zero sums across all samples
RNA_bac_controls <- subset_taxa(RNA_bac_controls, taxa_sums(RNA_bac_controls) > 0)

# Calculating Bray-Curtis distance matrix on raw counts
RNA_bac_controls_bray <- phyloseq::distance(RNA_bac_controls, method = "bray")

# Transforming raw counts to relative abundance
RNA_bac_controls_rel_abund <- transform_sample_counts(RNA_bac_controls, function(x) 100 * x / sum(x))

# Melting the phyloseq object into a data frame
RNA_bac_controls_rel_abund_melt <- psmelt(RNA_bac_controls_rel_abund)

# Assign ASV names and row names
RNA_bac_controls_rel_abund_melt <- RNA_bac_controls_rel_abund_melt %>%
  group_by(OTU) %>%
  mutate(ASV = paste0("ASV", cur_group_id())) %>%
  ungroup()

RNA_bac_controls_rel_abund_melt <- as.data.frame(RNA_bac_controls_rel_abund_melt)
# Renaming unknown taxa as 'Unknown'
for (i in 6:12) {
  RNA_bac_controls_rel_abund_melt[, i] <- as.character(RNA_bac_controls_rel_abund_melt[, i])
  RNA_bac_controls_rel_abund_melt[is.na(RNA_bac_controls_rel_abund_melt[, i]), i] <- "Unknown"
}

# Assign 'Others' to low-abundance taxa (treshold of 1%)
RNA_bac_controls_rel_abund_melt[RNA_bac_controls_rel_abund_melt$Abundance < 1, 6:12] <- "Others"

# Shorten overly long genus names
RNA_bac_controls_rel_abund_melt$Genus[RNA_bac_controls_rel_abund_melt$Genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "A-N-P-R"

# Add a column to specify the type of data
RNA_bac_controls_rel_abund_melt$type <- "16S cDNA"

# Visualization of relative abundance data of the bacterial controls with Genus as labels
ggplot(RNA_bac_controls_rel_abund_melt, aes(x = Sample, y = Abundance, fill = fct_reorder(Family, Abundance, .desc = TRUE))) +
  geom_bar(stat = "identity", color = "black") +
  facet_nested(.~treatment + type, scale = 'free_x') +
  labs(fill = "Phylum") +
  theme_bw() +
  ylab("Abundance (%)") +
  xlab("") +
  guides(fill = guide_legend(ncol = 4)) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.grid = element_line(color = "white")) +
  geom_text(aes(label = ifelse(Abundance > 5, as.character(Genus), '')), size = 3, position = position_stack(vjust = 0.6), show.legend = FALSE)

# Analyses on unique ASVs and taxa in the bacterial controls
length(unique(RNA_bac_controls_rel_abund_melt$ASV[RNA_bac_controls_rel_abund_melt$Class != "Others"]))
unique(RNA_bac_controls_rel_abund_melt$ASV[RNA_bac_controls_rel_abund_melt$Abundance > 0 & RNA_bac_controls_rel_abund_melt$treatment == "PS"])
unique(RNA_bac_controls_rel_abund_melt$ASV[RNA_bac_controls_rel_abund_melt$Abundance > 0 & RNA_bac_controls_rel_abund_melt$treatment == "VM"])
unique(RNA_bac_controls_rel_abund_melt$Genus[RNA_bac_controls_rel_abund_melt$Family ==  "Rhodobacteraceae" & RNA_bac_controls_rel_abund_melt$treatment == "PS"])
unique(RNA_bac_controls_rel_abund_melt$Species[RNA_bac_controls_rel_abund_melt$Family ==  "Rhodobacteraceae" & RNA_bac_controls_rel_abund_melt$treatment == "PS"])

any(RNA_bac_controls_rel_abund_melt$Family == "Mitochondria")
any(RNA_bac_controls_rel_abund_melt$Order == "Chloroplast")






# Analysis of the co-culture bacteria ------------------------------------------------------------------------------
## select only the treatments needed in the data analysis (N stands for "New" sequencing run of the same ("O") samples in 16S rRNA (cDNA))
DNA_experiment <- subset_samples(DNA, run %in% c("O", "env"))
RNA_experiment <- subset_samples(RNA, treatment!="control")

# Data Preparation for Amplicon Analysis in an Experiment

sample_data(RNA_experiment)$sample <- rownames(sample_data(RNA_experiment))
sample_data(DNA_experiment)$type <- "DNA"
sample_data(DNA_experiment) <- sample_data(DNA_experiment)[,c(4,3,1)]
## unite the object of RNA and DNA together
experiment_amplicon <- merge_phyloseq(DNA_experiment, RNA_experiment)
## Subsets only the co-culture samples without the outliers or redundant technical replicates with lower number of reads
experiment_amplicon <- subset_samples(experiment_amplicon,  ! sample_names(experiment_amplicon) %in% c("DVM1-O", "DVM1C-O", "PS4-O", "PS4C-O", "AX6", "AX6-B", "DVM4-B", "DVM4-C"))
## rename mislabelled sample
sample_names(experiment_amplicon)[sample_names(experiment_amplicon)=="D2P5-O"] <- "DP5-O"

## clean the experiment data from the Chloroplast and Mitochondria reads
experiment_amplicon_no_chloro_mito <- subset_taxa(experiment_amplicon, Order != "Chloroplast")
experiment_amplicon_no_chloro_mito <- subset_taxa(experiment_amplicon_no_chloro_mito, Family != "Mitochondria")
## check that they are removed
"Mitochondria" %in% as.character(unique(tax_table(experiment_amplicon_no_chloro_mito)[,5]))
"Chloroplast" %in% as.character(unique(tax_table(experiment_amplicon_no_chloro_mito)[,4]))
chloro <- subset_taxa(experiment_amplicon, Order == "Chloroplast")

# Display the proportion of Chloroplast and Mitochondria ASVs from the total per group
round(colSums(otu_table(chloro))*100/colSums(otu_table(experiment_amplicon)),2)
mean(round(colSums(otu_table(chloro))*100/colSums(otu_table(experiment_amplicon)),2)[1:4]) ## AX chloroplast content
mean(round(colSums(otu_table(chloro))*100/colSums(otu_table(experiment_amplicon)),2)[7:11]) ## DPS DNA chloroplast content
mean(round(colSums(otu_table(chloro))*100/colSums(otu_table(experiment_amplicon)),2)[12:16]) ## DVM DNA chloroplast content


# Remove taxa with zero sums across all samples, thus retaining only the taxa that have counts greater than zero
experiment_amplicon_no_chloro_mito <- subset_taxa(experiment_amplicon_no_chloro_mito, 
                                                  taxa_sums(experiment_amplicon_no_chloro_mito) > 0)

# Subset the samples based on 'treatment' values of "DPS" or "DVM" and the inocula samples ("env")
DPS_DVM_env <- subset_samples(experiment_amplicon_no_chloro_mito, 
                              treatment %in% c("DPS", "DVM") | 
                                grepl("env", sample_names(experiment_amplicon_no_chloro_mito)))

# Remove outliers by sample names "DPS6" and "DP5-O" from the subset
DPS_DVM_env_no_outlier <- subset_samples(DPS_DVM_env, 
                                         !sample_names(DPS_DVM_env) %in% c("DPS6", "DP5-O"))

# Fetch row names of the ASV table of the filtered dataset
seqs <- rownames(otu_table(DPS_DVM_env_no_outlier))

# Rename row names to start with "ASV" followed by a sequence number
names(seqs) <- paste0("ASV", seq(1:nrow(otu_table(DPS_DVM_env_no_outlier))))

# Apply the renamed taxa names to the OTU table.
taxa_names(DPS_DVM_env_no_outlier) <- names(seqs)

# Update 'sample' column in sample data. Concatenate 'type' and 'sample' columns
sample_data(DPS_DVM_env_no_outlier)$sample <- paste0(sample_data(DPS_DVM_env_no_outlier)$type, "_", 
                                                     sample_data(DPS_DVM_env_no_outlier)$sample)

# If the treatment is either "PS" or "VM", prepend "env_" to the 'sample' name, this is done to label the inocula for plotting
sample_data(DPS_DVM_env_no_outlier)$sample[sample_data(DPS_DVM_env_no_outlier)$treatment %in% c("PS", "VM")] <- 
  paste0("env_", sample_data(DPS_DVM_env_no_outlier)$sample[sample_data(DPS_DVM_env_no_outlier)$treatment %in% c("PS", "VM")])

# Remove any "B" from the 'sample' names (artifact of the labelling)
sample_data(DPS_DVM_env_no_outlier)$sample <- gsub("B", "", sample_data(DPS_DVM_env_no_outlier)$sample)




# 16S rRNA (DNA fraction) - Community composition ------------------------------------------------------------------------------------

# Calculate bray curtis distance matrix on raw counts for DNA samples in "DVM" and "DPS" treatments
DNA_experiment_bray <- phyloseq::distance(subset_samples(DPS_DVM_env_no_outlier, type=="DNA" & treatment %in% c("DVM", "DPS")), method = "bray")

# Convert sample_data to a data frame for further analyses
sampledf <- data.frame(sample_data(subset_samples(DPS_DVM_env_no_outlier, type=="DNA" & treatment %in% c("DVM", "DPS"))))

# Perform homogeneity of group dispersions test (betadisper)
beta_DNA <- betadisper(DNA_experiment_bray, sampledf$treatment)
# Conduct permutation test for betadisper
permutest(beta_DNA)

# Perform Adonis test to assess significant differences among treatment groups
adonis2(DNA_experiment_bray ~ treatment, data = sampledf)

# Convert the Bray-Curtis distance matrix to a data frame
DNA_experiment_bray_df <- as.data.frame(as.matrix(DNA_experiment_bray))

# Calculate the mean Bray-Curtis distance for DPS
mean(as.numeric(unlist(DNA_experiment_bray_df[1:5, 1:5])))

# Calculate the mean Bray-Curtis distance for DVM
mean(as.numeric(unlist(DNA_experiment_bray_df[6:10, 6:10])))


# 16S cDNA (RNA fraction) - Active comunity composition ------------------------------------------------------------------------------------
# Calculate the Bray-Curtis distance matrix for RNA samples in "DVM" and "DPS" treatments
RNA_experiment_bray <- phyloseq::distance(subset_samples(DPS_DVM_env_no_outlier, type=="RNA" & treatment %in% c("DVM", "DPS")), method = "bray")

# Convert sample_data to a data frame for subsequent analyses
sampledf_RNA <- data.frame(sample_data(subset_samples(DPS_DVM_env_no_outlier, type=="RNA" & treatment %in% c("DVM", "DPS"))))

# Perform homogeneity of group dispersions test (betadisper)
beta_RNA <- betadisper(RNA_experiment_bray, sampledf_RNA$treatment)
# Run permutation test for betadisper
permutest(beta_RNA)

# Conduct Adonis test to determine significant differences between treatment groups
adonis2(RNA_experiment_bray ~ treatment, data = sampledf_RNA, permutations = 999)

# Convert the Bray-Curtis distance matrix to a data frame
RNA_experiment_bray_df <- as.data.frame(as.matrix(RNA_experiment_bray))

# Calculate the mean Bray-Curtis distance for DPS
mean(as.numeric(unlist(RNA_experiment_bray_df[1:5, 1:5])))

# Calculate the mean Bray-Curtis distance for DVM
mean(as.numeric(unlist(RNA_experiment_bray_df[6:10, 6:10])))

## Diversity
## Rarefaction functions should be used for observed counts. Removing rare species before rarefaction can also give biased results. ?rarefy
DPS_DVM_env_no_outlier_rarefied <- rarefy_even_depth(DPS_DVM_env_no_outlier, sample.size = min(sample_sums(DPS_DVM_env_no_outlier)))

# Melt the rarefied phyloseq object and extract relevant columns
amplicon_wide_matrix_rarefied <- psmelt(DPS_DVM_env_no_outlier_rarefied) %>%
  select(OTU, sample, Abundance) %>%
  pivot_wider(names_from = sample, values_from = Abundance) %>%
  column_to_rownames(var = "OTU") %>%
  as.matrix()

# Transpose the matrix for downstream analyses
amplicon_wide_matrix_rarefied <- t(amplicon_wide_matrix_rarefied)

# Rarefy the amplicon-wide matrix
amplicon_wide_matrix_rarefied_rarefied <- rarefy(amplicon_wide_matrix_rarefied, min(rowSums(amplicon_wide_matrix_rarefied)))

# Generate rarefaction curve
rare_curve <- rarecurve(amplicon_wide_matrix_rarefied, step = 20, sample =  min(rowSums(amplicon_wide_matrix_rarefied)), col = "blue", cex = 0.6, tidy = TRUE)
# Extract treatment information from Site
rare_curve$treatment <- gsub('[[:digit:]]+', '', rare_curve$Site)
rare_curve$treatment <- gsub('[DR]NA_', '', rare_curve$treatment)
# Plot rarefaction curve using ggplot2
rare_curve_fig <- ggplot(rare_curve, aes(x = Sample, y = Species, group=Site, color=treatment)) +
  geom_line() + theme_minimal()

# Compute Shannon diversity indices
shannon_div <- diversity(amplicon_wide_matrix_rarefied, index = "shannon")

# Calculate mean and standard deviation of Shannon diversity indices for various treatments
mean(shannon_div[grepl("env_DNA_VM", names(shannon_div))]); sd(shannon_div[grepl("env_DNA_VM", names(shannon_div))])
mean(shannon_div[grepl("env_DNA_PS", names(shannon_div))]); sd(shannon_div[grepl("env_DNA_PS", names(shannon_div))])
mean(shannon_div[grepl("DNA_DVM", names(shannon_div))]); sd(shannon_div[grepl("DNA_DVM", names(shannon_div))])
mean(shannon_div[grepl("DNA_DPS", names(shannon_div))]); sd(shannon_div[grepl("DNA_DPS", names(shannon_div))])
mean(shannon_div[grepl("RNA_DVM", names(shannon_div))]); sd(shannon_div[grepl("RNA_DVM", names(shannon_div))])
mean(shannon_div[grepl("RNA_DPS", names(shannon_div))]); sd(shannon_div[grepl("RNA_DPS", names(shannon_div))])

wilcox.test(shannon_div[grepl("env_DNA_VM", names(shannon_div))], shannon_div[grepl("env_DNA_PS", names(shannon_div))], p.adjust.methods="BH", exact=FALSE)$p.value
wilcox.test(shannon_div[grepl("RNA_DVM", names(shannon_div))], shannon_div[grepl("RNA_DPS", names(shannon_div))], p.adjust.methods="BH", exact=FALSE)$p.value
wilcox.test(shannon_div[grepl("DNA_DVM", names(shannon_div))], shannon_div[grepl("DNA_DPS", names(shannon_div))], p.adjust.methods="BH", exact=FALSE)$p.value






# Mean Relative Abundance Calculations per group ------------------------------------------------------------------------------------

# Melt the phyloseq object to a data frame
DPS_DVM_env_melt <- psmelt(DPS_DVM_env_no_outlier)

# Rename missing or unknown taxonomic levels to 'Unknown'
for(i in 7:13){
  DPS_DVM_env_melt[,i] <- as.character(DPS_DVM_env_melt[,i])
  DPS_DVM_env_melt[is.na(DPS_DVM_env_melt[,i]),i] <- "Unknown"
}

# Retrieve unique OTUs for a specific condition and family
PS_Rhodobacteraceae_ASVs <- unique(DPS_DVM_env_melt$OTU[DPS_DVM_env_melt$type == "DNA" & DPS_DVM_env_melt$treatment == "PS" & DPS_DVM_env_melt$Family == "Rhodobacteraceae" & DPS_DVM_env_melt$Abundance > 0])

# Convert counts to relative abundance
DPS_DVM_env_relative_abundance_melt <- DPS_DVM_env_melt %>% 
  group_by(sample) %>% 
  mutate(Abundance = Abundance/sum(Abundance) * 100) %>%
  ungroup()

# Rename low abundance (<1%) taxa to 'Others'
DPS_DVM_env_relative_abundance_melt <- DPS_DVM_env_relative_abundance_melt[DPS_DVM_env_relative_abundance_melt$Abundance > 0, ]
DPS_DVM_env_relative_abundance_melt[DPS_DVM_env_relative_abundance_melt$Abundance < 1, c(7:13)] <- "Others"

## Calculations of mean relative abundances per group
## Calculate the mean phyla abundances per group
mean_phyla_amplicon_DPS_DNA <- DPS_DVM_env_relative_abundance_melt %>% 
  filter(treatment == "DPS" & type == "DNA") %>% 
  group_by(Phylum) %>% 
  summarize(mean = sum(Abundance)/5) %>% 
  arrange(desc(mean))

mean_phyla_amplicon_DVM_DNA <- DPS_DVM_env_relative_abundance_melt %>% 
  filter(treatment == "DVM" & type == "DNA") %>% 
  group_by(Phylum) %>% 
  summarize(mean = sum(Abundance)/5) %>% 
  arrange(desc(mean))

mean_Class_amplicon_DPS_DNA <- DPS_DVM_env_relative_abundance_melt %>% 
  filter(treatment == "DPS" & type == "DNA") %>% 
  group_by(Class) %>% 
  summarize(mean = sum(Abundance)/5) %>% 
  arrange(desc(mean))

mean_Class_amplicon_DVM_DNA <- DPS_DVM_env_relative_abundance_melt %>% 
  filter(treatment == "DVM" & type == "DNA") %>% 
  group_by(Class) %>% 
  summarize(mean = sum(Abundance)/5) %>% 
  arrange(desc(mean))


mean_Genus_amplicon_DPS_DNA <- DPS_DVM_env_relative_abundance_melt %>% 
  filter(treatment == "DPS" & type == "DNA") %>% 
  group_by(Family, Genus) %>% 
  summarize(mean = sum(Abundance)/5) %>% 
  arrange(desc(mean))

mean_Genus_amplicon_DVM_DNA <- DPS_DVM_env_relative_abundance_melt %>% 
  filter(treatment == "DVM" & type == "DNA") %>% 
  group_by(Family, Genus) %>% 
  summarize(mean = sum(Abundance)/5) %>% 
  arrange(desc(mean))

## Calculate the mean phyla abundances per group - bacterial activity
mean_phyla_amplicon_DPS_RNA <- DPS_DVM_env_relative_abundance_melt %>% 
  filter(treatment == "DPS" & type == "RNA") %>% 
  group_by(Phylum) %>% 
  summarize(mean = sum(Abundance)/5) %>% 
  arrange(desc(mean))

mean_phyla_amplicon_DVM_RNA <- DPS_DVM_env_relative_abundance_melt %>% 
  filter(treatment == "DVM" & type == "RNA") %>% 
  group_by(Phylum) %>% 
  summarize(mean = sum(Abundance)/5) %>% 
  arrange(desc(mean))

mean_Class_amplicon_DPS_RNA <- DPS_DVM_env_relative_abundance_melt %>% 
  filter(treatment == "DPS" & type == "RNA") %>% 
  group_by(Class) %>% 
  summarize(mean = sum(Abundance)/5) %>% 
  arrange(desc(mean))

mean_Class_amplicon_DVM_RNA <- DPS_DVM_env_relative_abundance_melt %>% 
  filter(treatment == "DVM" & type == "RNA") %>% 
  group_by(Class) %>% 
  summarize(mean = sum(Abundance)/5) %>% 
  arrange(desc(mean))


mean_Genus_amplicon_DPS_RNA <- DPS_DVM_env_relative_abundance_melt %>% 
  filter(treatment == "DPS" & type == "RNA") %>% 
  group_by(Family, Genus) %>% 
  summarize(mean = sum(Abundance)/5) %>% 
  arrange(desc(mean))

mean_Genus_amplicon_DVM_RNA <- DPS_DVM_env_relative_abundance_melt %>% 
  filter(treatment == "DVM" & type == "RNA") %>% 
  group_by(Family, Genus) %>% 
  summarize(mean = sum(Abundance)/5) %>% 
  arrange(desc(mean))



# RNA & DNA - Shared ASVs validative cross-contamination test ------------------------------------------------------------------------------------

shared_ASVs_RNA <- Reduce(intersect, list(unique(DPS_DVM_env_relative_abundance_melt$OTU[DPS_DVM_env_relative_abundance_melt$type=="RNA" & DPS_DVM_env_relative_abundance_melt$treatment=="DPS" & DPS_DVM_env_relative_abundance_melt$Phylum != "Others"]), unique(DPS_DVM_env_relative_abundance_melt$OTU[DPS_DVM_env_relative_abundance_melt$type=="RNA" & DPS_DVM_env_relative_abundance_melt$treatment=="DVM" & DPS_DVM_env_relative_abundance_melt$Phylum != "Others"])))
mean(DPS_DVM_env_relative_abundance_melt$Abundance[DPS_DVM_env_relative_abundance_melt$type=="RNA" & DPS_DVM_env_relative_abundance_melt$treatment=="DPS" & DPS_DVM_env_relative_abundance_melt$Phylum != "Others" & DPS_DVM_env_relative_abundance_melt$OTU %in% shared_ASVs_RNA])
mean(DPS_DVM_env_relative_abundance_melt$Abundance[DPS_DVM_env_relative_abundance_melt$type=="RNA" & DPS_DVM_env_relative_abundance_melt$treatment=="DVM" & DPS_DVM_env_relative_abundance_melt$Phylum != "Others" & DPS_DVM_env_relative_abundance_melt$OTU %in% shared_ASVs_RNA])
unique(DPS_DVM_env_relative_abundance_melt$Genus[DPS_DVM_env_relative_abundance_melt$type=="RNA" & DPS_DVM_env_relative_abundance_melt$treatment=="DVM" & DPS_DVM_env_relative_abundance_melt$Phylum != "Others" & DPS_DVM_env_relative_abundance_melt$OTU %in% shared_ASVs_RNA])

shared_ASVs_DNA <- Reduce(intersect, list(unique(DPS_DVM_env_relative_abundance_melt$OTU[DPS_DVM_env_relative_abundance_melt$type=="DNA" & DPS_DVM_env_relative_abundance_melt$treatment=="DPS" & DPS_DVM_env_relative_abundance_melt$Phylum != "Others"]), unique(DPS_DVM_env_relative_abundance_melt$OTU[DPS_DVM_env_relative_abundance_melt$type=="DNA" & DPS_DVM_env_relative_abundance_melt$treatment=="DVM" & DPS_DVM_env_relative_abundance_melt$Phylum != "Others"])))
mean(DPS_DVM_env_relative_abundance_melt$Abundance[DPS_DVM_env_relative_abundance_melt$type=="DNA" & DPS_DVM_env_relative_abundance_melt$treatment=="DPS" & DPS_DVM_env_relative_abundance_melt$Phylum != "Others" & DPS_DVM_env_relative_abundance_melt$OTU %in% shared_ASVs_DNA])
mean(DPS_DVM_env_relative_abundance_melt$Abundance[DPS_DVM_env_relative_abundance_melt$type=="DNA" & DPS_DVM_env_relative_abundance_melt$treatment=="DVM" & DPS_DVM_env_relative_abundance_melt$Phylum != "Others" & DPS_DVM_env_relative_abundance_melt$OTU %in% shared_ASVs_DNA])
unique(DPS_DVM_env_relative_abundance_melt$Genus[DPS_DVM_env_relative_abundance_melt$type=="DNA" & DPS_DVM_env_relative_abundance_melt$treatment=="DVM" & DPS_DVM_env_relative_abundance_melt$Phylum != "Others" & DPS_DVM_env_relative_abundance_melt$OTU %in% shared_ASVs_DNA])


shared_DNA_DVM_VM <- Reduce(intersect, list(unique(DPS_DVM_env_relative_abundance_melt$OTU[DPS_DVM_env_relative_abundance_melt$type=="DNA" & DPS_DVM_env_relative_abundance_melt$treatment=="VM"]), unique(DPS_DVM_env_relative_abundance_melt$OTU[DPS_DVM_env_relative_abundance_melt$type=="DNA" & DPS_DVM_env_relative_abundance_melt$treatment=="DVM" & DPS_DVM_env_relative_abundance_melt$Phylum != "Others"])))
mean(DPS_DVM_env_relative_abundance_melt$Abundance[DPS_DVM_env_relative_abundance_melt$type=="DNA" & DPS_DVM_env_relative_abundance_melt$treatment=="DVM" & DPS_DVM_env_relative_abundance_melt$Phylum != "Others" & DPS_DVM_env_relative_abundance_melt$OTU %in% shared_DNA_DVM_VM])
unique(DPS_DVM_env_relative_abundance_melt$Genus[DPS_DVM_env_relative_abundance_melt$type=="DNA" & DPS_DVM_env_relative_abundance_melt$treatment=="DVM" & DPS_DVM_env_relative_abundance_melt$Phylum != "Others" & DPS_DVM_env_relative_abundance_melt$OTU %in% shared_DNA_DVM_VM])
mean(DPS_DVM_env_relative_abundance_melt$Abundance[DPS_DVM_env_relative_abundance_melt$type=="DNA" & DPS_DVM_env_relative_abundance_melt$treatment=="VM" & DPS_DVM_env_relative_abundance_melt$OTU %in% shared_DNA_DVM_VM])

shared_DNA_DPS_PS <- Reduce(intersect, list(unique(DPS_DVM_env_relative_abundance_melt$OTU[DPS_DVM_env_relative_abundance_melt$type=="DNA" & DPS_DVM_env_relative_abundance_melt$treatment=="PS"]), unique(DPS_DVM_env_relative_abundance_melt$OTU[DPS_DVM_env_relative_abundance_melt$type=="DNA" & DPS_DVM_env_relative_abundance_melt$treatment=="DPS" & DPS_DVM_env_relative_abundance_melt$Phylum != "Others"])))
mean(DPS_DVM_env_relative_abundance_melt$Abundance[DPS_DVM_env_relative_abundance_melt$type=="DNA" & DPS_DVM_env_relative_abundance_melt$treatment=="DPS" & DPS_DVM_env_relative_abundance_melt$Phylum != "Others" & DPS_DVM_env_relative_abundance_melt$OTU %in% shared_DNA_DPS_PS])
unique(DPS_DVM_env_relative_abundance_melt$Genus[DPS_DVM_env_relative_abundance_melt$type=="DNA" & DPS_DVM_env_relative_abundance_melt$treatment=="DPS" & DPS_DVM_env_relative_abundance_melt$Phylum != "Others" & DPS_DVM_env_relative_abundance_melt$OTU %in% shared_DNA_DPS_PS])
mean(DPS_DVM_env_relative_abundance_melt$Abundance[DPS_DVM_env_relative_abundance_melt$type=="DNA" & DPS_DVM_env_relative_abundance_melt$treatment=="PS" & DPS_DVM_env_relative_abundance_melt$OTU %in% shared_DNA_DPS_PS])



# Shared ASVs in co-cultures ------------------------------------------------------------------------------------
## Counts the number of shared ASVs between groups

# Helper function to filter unique ASVs based on type, treatment, and Phylum
filter_unique_ASVs <- function(data, type, treatment, exclude_others = TRUE) {
  if (exclude_others) {
    unique_ASVs <- unique(data$OTU[
      data$type == type &
        data$treatment == treatment &
        data$Phylum != "Others"
      ])
  } else {
    unique_ASVs <- unique(data$OTU[
      data$type == type &
        data$treatment == treatment
      ])
  }
  return(unique_ASVs)
}

# Shared ASVs between inocula VM and PS under DNA type
shared_envDNA_ASVs <- Reduce(intersect, list(
  filter_unique_ASVs(DPS_DVM_env_relative_abundance_melt, "DNA", "VM"),
  filter_unique_ASVs(DPS_DVM_env_relative_abundance_melt, "DNA", "PS")
))

# Shared ASVs between PS inoculum and DPS co-culture under DNA type
shared_envPS_DPS_ASVs <- Reduce(intersect, list(
  filter_unique_ASVs(DPS_DVM_env_relative_abundance_melt, "DNA", "PS", FALSE),
  filter_unique_ASVs(DPS_DVM_env_relative_abundance_melt, "DNA", "DPS", FALSE)
))

round(mean(DPS_DVM_env_relative_abundance_melt$Abundance[DPS_DVM_env_relative_abundance_melt$type == "DNA" &
                                                  DPS_DVM_env_relative_abundance_melt$treatment == "DPS" &
                                                  DPS_DVM_env_relative_abundance_melt$OTU %in% shared_envPS_DPS_ASVs]), 2)

round(mean(DPS_DVM_env_relative_abundance_melt$Abundance[DPS_DVM_env_relative_abundance_melt$type == "DNA" &
                                                           DPS_DVM_env_relative_abundance_melt$treatment == "PS" &
                                                           DPS_DVM_env_relative_abundance_melt$OTU %in% shared_envPS_DPS_ASVs]), 2)



# Shared ASVs between VM inoculum and DVM co-culture under DNA type
shared_envVM_DVM_ASVs <- Reduce(intersect, list(
  filter_unique_ASVs(DPS_DVM_env_relative_abundance_melt, "DNA", "VM", FALSE),
  filter_unique_ASVs(DPS_DVM_env_relative_abundance_melt, "DNA", "DVM")
))

round(mean(DPS_DVM_env_relative_abundance_melt$Abundance[DPS_DVM_env_relative_abundance_melt$type == "DNA" &
                                                           DPS_DVM_env_relative_abundance_melt$treatment == "DVM" &
                                                           DPS_DVM_env_relative_abundance_melt$OTU %in% shared_envVM_DVM_ASVs]), 2)

round(mean(DPS_DVM_env_relative_abundance_melt$Abundance[DPS_DVM_env_relative_abundance_melt$type == "DNA" &
                                                           DPS_DVM_env_relative_abundance_melt$treatment == "VM" &
                                                           DPS_DVM_env_relative_abundance_melt$OTU %in% shared_envVM_DVM_ASVs]), 2)


# Shared ASVs between VM inoculum and DPS co-culture under DNA type
shared_envVM_DPS_ASVs <- Reduce(intersect, list(
  filter_unique_ASVs(DPS_DVM_env_relative_abundance_melt, "DNA", "VM", FALSE),
  filter_unique_ASVs(DPS_DVM_env_relative_abundance_melt, "DNA", "DPS")
))

round(mean(DPS_DVM_env_relative_abundance_melt$Abundance[DPS_DVM_env_relative_abundance_melt$type == "DNA" &
                                                           DPS_DVM_env_relative_abundance_melt$treatment == "DPS" &
                                                           DPS_DVM_env_relative_abundance_melt$OTU %in% shared_envVM_DPS_ASVs]), 2)

round(mean(DPS_DVM_env_relative_abundance_melt$Abundance[DPS_DVM_env_relative_abundance_melt$type == "DNA" &
                                                           DPS_DVM_env_relative_abundance_melt$treatment == "VM" &
                                                           DPS_DVM_env_relative_abundance_melt$OTU %in% shared_envVM_DPS_ASVs]), 2)


# Shared ASVs between DVM and DPS co-cultures under DNA type
shared_DNA_ASVs <- Reduce(intersect, list(
  filter_unique_ASVs(DPS_DVM_env_relative_abundance_melt, "DNA", "DVM"),
  filter_unique_ASVs(DPS_DVM_env_relative_abundance_melt, "DNA", "DPS")
))

round(mean(DPS_DVM_env_relative_abundance_melt$Abundance[DPS_DVM_env_relative_abundance_melt$type == "DNA" &
                                                           DPS_DVM_env_relative_abundance_melt$treatment == "DPS" &
                                                           DPS_DVM_env_relative_abundance_melt$OTU %in% shared_DNA_ASVs]), 2)

round(mean(DPS_DVM_env_relative_abundance_melt$Abundance[DPS_DVM_env_relative_abundance_melt$type == "DNA" &
                                                           DPS_DVM_env_relative_abundance_melt$treatment == "DVM" &
                                                           DPS_DVM_env_relative_abundance_melt$OTU %in% shared_DNA_ASVs]), 2)



# Shared ASVs between DVM and DPS co-cultures under RNA type
shared_RNA_ASVs <- Reduce(intersect, list(
  filter_unique_ASVs(DPS_DVM_env_relative_abundance_melt, "RNA", "DVM"),
  filter_unique_ASVs(DPS_DVM_env_relative_abundance_melt, "RNA", "DPS")
))

round(mean(DPS_DVM_env_relative_abundance_melt$Abundance[DPS_DVM_env_relative_abundance_melt$type == "RNA" &
                                                           DPS_DVM_env_relative_abundance_melt$treatment == "DPS" &
                                                           DPS_DVM_env_relative_abundance_melt$OTU %in% shared_RNA_ASVs]), 2)

round(mean(DPS_DVM_env_relative_abundance_melt$Abundance[DPS_DVM_env_relative_abundance_melt$type == "RNA" &
                                                           DPS_DVM_env_relative_abundance_melt$treatment == "DVM" &
                                                           DPS_DVM_env_relative_abundance_melt$OTU %in% shared_RNA_ASVs]), 2)




# top 5 most abundant ASVs per co-culture sample -----------------------------------------------------------------------------------------------
# Initialize a DataFrame to hold RNA-specific data
experiment_amplicon_RNA <- subset(DPS_DVM_env_relative_abundance_melt, type == "RNA")

# Add a column to classify each sample's ASVs as 'top5' or 'none'
experiment_amplicon_RNA$subtype <- "none"

# Loop through each unique sample to identify its top 5 most abundant ASVs
unique_samples <- unique(experiment_amplicon_RNA$Sample)
for(i in 1:length(unique_samples)) {
  
  # Get the current unique sample name
  current_sample <- unique_samples[i]
  
  # Identify the top 5 most abundant ASVs for this sample
  top5_ASVs <- experiment_amplicon_RNA[experiment_amplicon_RNA$Sample == current_sample & experiment_amplicon_RNA$Family != "Others", ] %>% 
    group_by(Sample) %>% 
    arrange(desc(Abundance)) %>% 
    slice_head(n = 5) %>% 
    ungroup()
  
  top5_ASVs <- unique(top5_ASVs$OTU)
  
  # Update 'subtype' to 'top5' for these ASVs
  experiment_amplicon_RNA$subtype[experiment_amplicon_RNA$Sample == current_sample & experiment_amplicon_RNA$OTU %in% top5_ASVs] <- "top5"
}

# Modify OTU labels based on subtype and Family conditions
experiment_amplicon_RNA$OTU <- ifelse(
  experiment_amplicon_RNA$subtype == "top5",  experiment_amplicon_RNA$OTU, "other ASVs")

# Convert OTU to factor based on pre-defined family_colors
experiment_amplicon_RNA$OTU <- as.character(experiment_amplicon_RNA$OTU)

# Aggregate data
experiment_amplicon_RNA_aggregated <- aggregate(
  Abundance ~ OTU + Sample + treatment,
  data = experiment_amplicon_RNA[, c(1, 2, 5, 3)],
  FUN = sum
) %>%
  # Rename the column names to match the original DataFrame
  setNames(colnames(experiment_amplicon_RNA[, c(1, 2, 5, 3)]))



# Community visualization ------------------------------------------------------------------------------------

## Aggregates the counts for the second time to group all Others
DPS_DVM_env_rel_abund_melt_family <- aggregate(
  Abundance ~ OTU + Class + Family + sample + type + treatment,
  data = DPS_DVM_env_relative_abundance_melt[,c(3,4,5,6,9,11,1)], # ,5
  sum
) %>%
  # Rename column names
  setNames(colnames(DPS_DVM_env_relative_abundance_melt[,c(1,9,11,6,4,5,3)]))

colnames(DPS_DVM_env_rel_abund_melt_family)[4] <- "Sample"

DPS_DVM_env_rel_abund_melt_family$seq <- ifelse(grepl("DNA", DPS_DVM_env_rel_abund_melt_family$Sample), "16S rRNA", "16S cDNA")
DPS_DVM_env_rel_abund_melt_family$Sample <- gsub("[RD]NA", "", DPS_DVM_env_rel_abund_melt_family$Sample)
DPS_DVM_env_rel_abund_melt_family$Sample <- gsub("_", "", DPS_DVM_env_rel_abund_melt_family$Sample)


## Add the bacterial control panel
RNA_bac_contols_rel_abund_melt_vis_subset <- RNA_bac_controls_rel_abund_melt[,c(13, 8,10, 2, 4, 5, 3, 4)]
RNA_bac_contols_rel_abund_melt_vis_subset$type <- "RNA"
colnames(RNA_bac_contols_rel_abund_melt_vis_subset)[1] <- "OTU"
colnames(RNA_bac_contols_rel_abund_melt_vis_subset)[8] <- "seq"

## merge the bacterial controls with the data
DPS_DVM_env_rel_abund_melt_family <- rbind(DPS_DVM_env_rel_abund_melt_family, RNA_bac_contols_rel_abund_melt_vis_subset)
## rename the taxa names to the latest change
DPS_DVM_env_rel_abund_melt_family$Class[DPS_DVM_env_rel_abund_melt_family$Class == "Actinobacteria"] <- "Actinobacteriota"
DPS_DVM_env_rel_abund_melt_family$Class[DPS_DVM_env_rel_abund_melt_family$Class == "Bacteroidetes"] <- "Bacteroidota"
DPS_DVM_env_rel_abund_melt_family$Class[DPS_DVM_env_rel_abund_melt_family$Class == "Planctomycetes"]  <- "Planctomycetota"


## For easier visualization, aggregate all the "Others" ASVs into one panel
DPS_DVM_env_rel_abund_melt_family_vis <- DPS_DVM_env_rel_abund_melt_family[DPS_DVM_env_rel_abund_melt_family$Family != "Others",]
DPS_DVM_env_rel_abund_melt_family_others <- DPS_DVM_env_rel_abund_melt_family[DPS_DVM_env_rel_abund_melt_family$Family == "Others",]
DPS_DVM_env_rel_abund_melt_family_others_aggregate <- aggregate(DPS_DVM_env_rel_abund_melt_family_others$Abundance, by = list(DPS_DVM_env_rel_abund_melt_family_others$Class, DPS_DVM_env_rel_abund_melt_family_others$Family, DPS_DVM_env_rel_abund_melt_family_others$Sample, DPS_DVM_env_rel_abund_melt_family_others$type, DPS_DVM_env_rel_abund_melt_family_others$treatment, DPS_DVM_env_rel_abund_melt_family_others$seq), sum)
DPS_DVM_env_rel_abund_melt_family_others_aggregate$OTU <- "Others"

## Unite again the aggregated table with the ASV-level table
DPS_DVM_env_rel_abund_melt_family_others_aggregate <- DPS_DVM_env_rel_abund_melt_family_others_aggregate[,c(8, 1:5,7,6)]
colnames(DPS_DVM_env_rel_abund_melt_family_others_aggregate) <- colnames(DPS_DVM_env_rel_abund_melt_family)
DPS_DVM_env_rel_abund_melt_family_vis <- rbind(DPS_DVM_env_rel_abund_melt_family_vis, DPS_DVM_env_rel_abund_melt_family_others_aggregate)
## Substitute the "Unknown" with the Class name
for(i in unique(DPS_DVM_env_rel_abund_melt_family_vis$Class)){
  if(any(grepl("Unknown", DPS_DVM_env_rel_abund_melt_family_vis$Family[DPS_DVM_env_rel_abund_melt_family_vis$Class == i]))){
    str <- paste0("Unknown ", i)
    DPS_DVM_env_rel_abund_melt_family_vis$Family[grepl("Unknown", DPS_DVM_env_rel_abund_melt_family_vis$Family) & DPS_DVM_env_rel_abund_melt_family_vis$Class == i] <- str
  }
}

## Prepare for plotting
DPS_DVM_env_rel_abund_melt_family_vis <- DPS_DVM_env_rel_abund_melt_family_vis[order(DPS_DVM_env_rel_abund_melt_family_vis$Class, DPS_DVM_env_rel_abund_melt_family_vis$Family), ]
DPS_DVM_env_rel_abund_melt_family_vis$Family <- factor(DPS_DVM_env_rel_abund_melt_family_vis$Family, levels=unique(DPS_DVM_env_rel_abund_melt_family_vis$Family))
DPS_DVM_env_rel_abund_melt_family_vis$OTU[! DPS_DVM_env_rel_abund_melt_family_vis$OTU %in% unique(experiment_amplicon_RNA_aggregated$OTU) | DPS_DVM_env_rel_abund_melt_family_vis$type == "DNA"] <- ""
DPS_DVM_env_rel_abund_melt_family_vis$OTU <- gsub("ASV", "", DPS_DVM_env_rel_abund_melt_family_vis$OTU)
DPS_DVM_env_rel_abund_melt_family_vis$type2 <- ifelse(DPS_DVM_env_rel_abund_melt_family_vis$treatment %in% c("DPS", "DVM"), "Co-cultures", "Inocula")
DPS_DVM_env_rel_abund_melt_family_vis$type2[DPS_DVM_env_rel_abund_melt_family_vis$Sample %in% c("PS4", "PS5", "PS6", "VM1", "VM2")] <- "Bacterial Controls"
DPS_DVM_env_rel_abund_melt_family_vis$type2 <- factor(DPS_DVM_env_rel_abund_melt_family_vis$type2, levels=c("Inocula", "Co-cultures", "Bacterial Controls"))
DPS_DVM_env_rel_abund_melt_family_vis$seq <- factor(DPS_DVM_env_rel_abund_melt_family_vis$seq, levels=c("16S rRNA", "16S cDNA"))

## Generates a common color palette as gradient per Class
class_colors <- c("orange2", "hotpink", "mediumpurple3", "green1", "tomato1", "brown","lemonchiffon2", "maroon3", "palegreen", "mintcream", "plum3")
names(class_colors) <- unique(DPS_DVM_env_rel_abund_melt_family_vis$Class)

family_gradient_colors <- rep(NA, length(unique(DPS_DVM_env_rel_abund_melt_family_vis$Family)))
names(family_gradient_colors) <- unique(DPS_DVM_env_rel_abund_melt_family_vis$Family)
family_gradient_colors <- family_gradient_colors[levels(DPS_DVM_env_rel_abund_melt_family_vis$Family)]

for(i in unique(DPS_DVM_env_rel_abund_melt_family_vis$Class)){
  vec1 <- unique(DPS_DVM_env_rel_abund_melt_family_vis$Family[DPS_DVM_env_rel_abund_melt_family_vis$Class == i])
  family_gradient_colors[names(family_gradient_colors) %in% vec1] <- class_colors[names(class_colors) == i]
  if(i == "Bacteroidia"){
      col1 <- class_colors[names(class_colors) == i]
      color_vec2 <-  colorRampPalette(c(col1, "darkolivegreen"))(length(vec1))
      family_gradient_colors[names(family_gradient_colors) %in% vec1] <- color_vec2
    }
    if(i == "Alphaproteobacteria"){
      col1 <- class_colors[names(class_colors) == i]
      color_vec2 <-  colorRampPalette(c(col1, "lightskyblue", "slategray3", "royalblue", "turquoise", "navy"))(length(vec1))
      family_gradient_colors[names(family_gradient_colors) %in% vec1] <- color_vec2
    }
    if(i == "Bdellovibrionia"){
      col1 <- class_colors[names(class_colors) == i]
      color_vec2 <-  colorRampPalette(c(col1, "red"))(length(vec1))
      family_gradient_colors[names(family_gradient_colors) %in% vec1] <- color_vec2
    }
    if(i == "Gammaproteobacteria"){
      col1 <- class_colors[names(class_colors) == i]
      color_vec2 <-  colorRampPalette(c(col1, "gold"))(length(vec1))
      family_gradient_colors[names(family_gradient_colors) %in% vec1] <- color_vec2
    }
}
## Renames Rhodobacteraceae as distinct color
family_gradient_colors[8] <- "grey70"

DPS_DVM_env_rel_abund_melt_family_vis_aggregated <- aggregate(DPS_DVM_env_rel_abund_melt_family_vis$Abundance, by = list(DPS_DVM_env_rel_abund_melt_family_vis$OTU, DPS_DVM_env_rel_abund_melt_family_vis$Class, DPS_DVM_env_rel_abund_melt_family_vis$Family, DPS_DVM_env_rel_abund_melt_family_vis$Sample, DPS_DVM_env_rel_abund_melt_family_vis$type, DPS_DVM_env_rel_abund_melt_family_vis$treatment, DPS_DVM_env_rel_abund_melt_family_vis$seq, DPS_DVM_env_rel_abund_melt_family_vis$type2), sum)
colnames(DPS_DVM_env_rel_abund_melt_family_vis_aggregated) <- colnames(DPS_DVM_env_rel_abund_melt_family_vis)[c(1:5,6,8,9,7)]
DPS_DVM_env_rel_abund_melt_family_vis_aggregated$treatment2 <- ifelse(grepl("PS", DPS_DVM_env_rel_abund_melt_family_vis_aggregated$treatment), "Foreign", "Native")

experiment_fig <- ggplot(DPS_DVM_env_rel_abund_melt_family_vis_aggregated, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", color = "grey20") + scale_fill_manual(values=family_gradient_colors) + facet_nested(.~type2+treatment2+seq, scale = 'free_x', space = "free") + 
  labs(fill = "Family") + theme_bw() + 
  ylab("Abundance (%)") + xlab("") + 
  guides(fill=guide_legend(ncol=1)) +   
  geom_text(aes(label = OTU), size=1.4, colour = "black", position = position_stack(vjust = 0.5)) + 
  theme(strip.text = element_text(size = 6), legend.position="right", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid = element_line(color = "white"), strip.background =element_rect(fill="white"))


## Save the figure in different formats, the legend separately for easier visualization
#topptx(experiment_fig, "figure.pptx", width=8, height=6)
#topptx(experiment_fig, "legend.pptx", width=6, height=8)


