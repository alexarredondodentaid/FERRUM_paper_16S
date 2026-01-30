# List of packages that will be needed:

library(dada2) # BiocManager
library(phyloseq) # BiocManager
library(ANCOMBC) # BiocManager
library(Biostrings) # BiocManager
library(biomformat) # BiocManager
library(microbiome) # BiocManager
library(vegan) # CRAN
library(dplyr) # CRAN
library(ggplot2) # CRAN
library(emmeans) # CRAN
library(ggpubr) # CRAN
library(lmerTest) # CRAN
library(DT) # CRAN
library(indicspecies) # CRAN
library(tidyr) # CRAN
library(igraph) # CRAN
library(Matrix) # CRAN
library(devtools) # CRAN
library(ggraph) # CRAN
library(compositions) # CRAN
library(ggrepel) # CRAN
library(circlize) # CRAN
library(microViz) # github
library(pairwiseAdonis) # github
library(SpiecEasi) # github

# Load the .fastq files using the path they're in:
path <- "/home/alexarredondo/Documents/2025_Nitrats_v6/fastq/"
list.files(path) # make sure the files in the path are the ones you need

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names=TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names=TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), "[", 1)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Reads have been precut at the 5' end using fastqc. Here we cut the 3' end (use the values in trimLeft that apply to your data).
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(17,21),
                     maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Denoise.
dadaFs <- dada(filtFs, err=errF, multithread = TRUE)
dadaRs <- dada(filtRs, err=errR, multithread = TRUE)

# Merge the forward and reverse reads.
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

# Remove reads shorter than 400, since they are probably not merged.
seqtab <- makeSequenceTable(mergers)

llargada <- nchar(getSequences(seqtab))>400
seqtab<- seqtab[,llargada]

# Prepare a table to assess how many reads you've lost in the process.
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", verbose=TRUE, multithread = TRUE)
getN <- function(x) sum(getUniques(x))
track_v1 <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track_v1) <- c("input", "filtrades", "denoisedF", "denoisedR", "ajuntades", "sense quimeres")
rownames(track_v1) <- sample.names
track_v1 <- track_v1[-199,]
View(track_v1)

# This histogram let you see the distribution of reads per sample in the last and first step of the process, respectively.
hist(track_v1[,1], main = "Reads per mostra sense filtrar", xlab = "n reads")
hist(track_v1[,6], main = "Reads per mostra filtrats", xlab = "n reads")

# Apply your taxonomy database to the asv_table. In this case we will be using the eHOMD db v15.22.
taxa <- assignTaxonomy(seqtab.nochim, "/home/ordinador1/Documents/Alex/Software/ehomd_dada2/possible_nou/eHOMD_RefSeq_dada2_V15.22.fasta.gz", tryRC = TRUE)
taxa_spp <- dada2::addSpecies(taxa, "/home/ordinador1/Documents/Alex/Software/ehomd_dada2/eHOMD_RefSeq_dada2_assign_species_V15.22.fasta.gz", verbose = TRUE)

# Now we prepare the phyloseq object. We begin by loading the metadata object.
metadata_df <- read.csv("metadata.csv", sep = ";")
rownames(metadata_df) <- rownames(otu_table(seqtab.nochim, taxa_are_rows = FALSE)) # make sure that the rownames of your metadata file and otu table are the same

# Here we add a new variable that merges timepoint and treatment.
metadata_df$Timepoint_Treatment <- c(paste(metadata_df$Timepoint, metadata_df$Treatment, sep = "_"))

# We add the number of reads for sample, obtained from the table that we did previously.
track_v1_df <- data.frame(track_v1)
metadata_df$Reads <- c(track_v1_df$sense.quimeres, 0) # we add this 0 to account for the unintegrated samples.

# Factorize all the variables that need to.
metadata_df <- metadata_df %>%
  mutate(across(c(Timepoint, Treatment, ID, Data_extraccio, lot_kit_extraccio,
                  extractor, tanda_pcr_16s, tanda_rentat_1, tanda_pcr_index,
                  tanda_rentat_2, Timepoint_Treatment), as.factor))

# Now we create the object for the genera and spp levels.
ps_genera_v1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data(metadata_df), tax_table(taxa))
ps_spp_v1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data(metadata_df), tax_table(taxa_spp))

# Let's remove the unintegrated samples (this could have been done previously).
ps_genera_v1 <- subset_samples(ps_genera_v1, Sample != "Unintegrated")
ps_spp_v1 <- subset_samples(ps_spp_v1, Sample != "Unintegrated")

# Use the taxglom function for those analyses that need it. 
ps_glom_genus_v1 <- tax_glom(ps_genera_v1, taxrank = "Genus")
ps_glom_spp_v1 <- tax_glom(ps_spp_v1, taxrank = "Species")

# We need to check how the sequencing coverage affected the diversity of our samples.
rarecurve(as.matrix(data.frame(otu_table(ps_genera_v1))), step = 50, cex = 0.5)

# There is a sample that took to many reads and does not allow us to assess the rest of them. Let's remove it.
rarecurve(as.matrix(data.frame(otu_table(ps_genera_v1))[-78,]), step = 50, cex = 0.5)

# We performe the rarefaction curves for all the groups of the study:
rarecurve(as.matrix(data.frame(otu_table(ps_genera_v1))[c(1:22),]), step = 50, cex = 0.5)
rarecurve(as.matrix(data.frame(otu_table(ps_genera_v1))[c(23:44),]), step = 50, cex = 0.5)
rarecurve(as.matrix(data.frame(otu_table(ps_genera_v1))[c(45:66),]), step = 50, cex = 0.5)
rarecurve(as.matrix(data.frame(otu_table(ps_genera_v1))[c(67:77,79:88),]), step = 50, cex = 0.5)
rarecurve(as.matrix(data.frame(otu_table(ps_genera_v1))[c(89:110),]), step = 50, cex = 0.5)
rarecurve(as.matrix(data.frame(otu_table(ps_genera_v1))[c(111:132),]), step = 50, cex = 0.5)
rarecurve(as.matrix(data.frame(otu_table(ps_genera_v1))[c(133:154),]), step = 50, cex = 0.5)
rarecurve(as.matrix(data.frame(otu_table(ps_genera_v1))[c(155:176),]), step = 50, cex = 0.5)
rarecurve(as.matrix(data.frame(otu_table(ps_genera_v1))[c(177:198),]), step = 50, cex = 0.5)

# According to the rarefaction curves, we need to remove at least 9 samples due to the lack of stabilization of richness: 7d-CDQ-18, 7d-KNO3-18, Fi-CDQ-5, Fi-CDQ-3, Fi-CDQ-22, Fi-KNO3-14, Fi-KNO3-21, Fi-MQ-18, Fi-MQ-22.
mostres_a_borrar <- c("7d-CDQ-18", "7d-KNO3-18", "Fi-CDQ-5", "Fi-CDQ-3", "Fi-CDQ-22", 
                      "Fi-KNO3-14", "Fi-KNO3-21", "Fi-MQ-18", "Fi-MQ-22")
ps_genera_v1_pruned <- subset_samples(ps_genera_v1, !(Sample %in% mostres_a_borrar))
ps_spp_v1_pruned <- subset_samples(ps_spp_v1, !(Sample %in% mostres_a_borrar))

ps_glom_genus_v1_pruned <- tax_glom(ps_genera_v1_pruned, taxrank = "Genus")
ps_glom_spp_v1_pruned <- tax_glom(ps_spp_v1_pruned, taxrank = "Species")
```

####################################### Alpha diversity #######################################

alpha_diversity <- estimate_richness(ps_genera_v1_pruned, split = TRUE, measures = NULL)
table_alpha <- data.frame(alpha_diversity, data.frame(sample_data(ps_genera_v1_pruned)))

table_alpha$Treatment <- factor(table_alpha$Treatment, levels = c("Basal", "MQ", "KNO3", "Base", "CDQ"))
table_alpha$Timepoint <- factor(table_alpha$Timepoint, levels = c("0d", "7d", "14d"))

# Let's visualize it with some boxplots.
box_chao1 <- ggplot(table_alpha, aes(x = Treatment , y = Chao1, fill = Timepoint)) +
  geom_boxplot() + 
  theme_minimal() +
  scale_fill_manual(values = c("chartreuse1", "chocolate", "cadetblue4"))

box_shannon <- ggplot(table_alpha, aes(x = Treatment , y = Shannon, fill = Timepoint)) +
  geom_boxplot() + 
  theme_minimal() +
  scale_fill_manual(values = c("chartreuse1", "chocolate", "cadetblue4"))

box_simpson <- ggplot(table_alpha, aes(x = Treatment , y = Simpson, fill = Timepoint)) +
  geom_boxplot() + 
  theme_minimal() +
  scale_fill_manual(values = c("chartreuse1", "chocolate", "cadetblue4"))

ggarrange(box_chao1, box_shannon, box_simpson, ncol = 3, common.legend = TRUE) # this is just to plot them all together.


# Let's check for statistical significant differences.
table_alpha$Reads_scaled <- scale(table_alpha$Reads) # Since the number of reads is quite diverse, we need to scale them. 

chao1 <- lm(Chao1 ~ Treatment*Timepoint + ID + Data_extraccio + 
              lot_kit_extraccio + extractor + tanda_pcr_16s + tanda_rentat_1 +
              tanda_pcr_index + tanda_rentat_2 + Reads_scaled, data = table_alpha)
shapiro.test(residuals(chao1)) # Check for the normalitiy of the residuals. Since they are not normal, we cannot use an lm, we have to use a glm.

chao1_glm <- glmer(Chao1 ~ Treatment + Timepoint + (1|ID) + 
                     Reads_scaled + lot_kit_extraccio,
                   family = Gamma(link = "log"), data = table_alpha) # Variables without significance have been deleted from the model.
summary(chao1_glm)

# We check the pairwise comparisons using emmeans.
chao1_Timepoint_emmeans <- data.frame(emmeans::emmeans(chao1_glm, specs = pairwise ~ Timepoint, 
                                                       data = table_alpha, rg.limit = 3e5)$emmeans)
chao1_Timepoint_contrasts <- data.frame(emmeans::emmeans(chao1_glm, specs = pairwise ~ Timepoint, 
                                                         data = table_alpha, rg.limit = 3e5, nesting = NULL)$contrasts)

chao1_Treatment_emmeans <- data.frame(emmeans::emmeans(chao1_glm, specs = pairwise ~ Treatment, 
                                                       data = table_alpha, rg.limit = 3e5)$emmeans)
chao1_Treatment_contrasts <- data.frame(emmeans::emmeans(chao1_glm, specs = pairwise ~ Treatment, 
                                                         data = table_alpha, rg.limit = 3e5)$contrasts)

chao1_Treat_Time_emmeans <- data.frame(emmeans::emmeans(chao1_glm, specs = pairwise ~ Treatment + Timepoint, 
                                                        data = table_alpha, rg.limit = 3e5)$emmeans)
chao1_Treat_Time_contrasts <- data.frame(emmeans::emmeans(chao1_glm, specs = pairwise ~ Treatment + Timepoint, 
                                                          data = table_alpha, rg.limit = 3e5)$contrasts)

# Export the pairwise comparisons to a .csv file to work with them outside R.
write.table(chao1_Treat_Time_emmeans, file = "/home/alexarredondo/Documents/2025_Nitrats_v6/chao1_emmeans.csv", sep = ";")
write.table(chao1_Treat_Time_contrasts, file = "/home/alexarredondo/Documents/2025_Nitrats_v6/chao1_contrasts.csv", sep = ";")

# Repeat the process for Shannon and Simpson indices.
shannon <- lm(Shannon ~ Treatment*Timepoint + ID + Data_extraccio + 
                lot_kit_extraccio + extractor + tanda_pcr_16s + tanda_rentat_1 +
                tanda_pcr_index + tanda_rentat_2 + Reads_scaled, data = table_alpha)
shapiro.test(residuals(shannon)) 

shannon_glm <- glmer(Shannon ~ Treatment + Timepoint + (1|ID) + 
                       Reads_scaled + lot_kit_extraccio,
                     family = Gamma(link = "log"), data = table_alpha)

summary(shannon_glm)

shannon_Timepoint_emmeans <- data.frame(emmeans::emmeans(shannon_glm, specs = pairwise ~ Timepoint, 
                                                         data = table_alpha, rg.limit = 3e5)$emmeans)
shannon_Timepoint_contrasts <- data.frame(emmeans::emmeans(shannon_glm, specs = pairwise ~ Timepoint, 
                                                           data = table_alpha, rg.limit = 3e5, nesting = NULL)$contrasts)

shannon_Treatment_emmeans <- data.frame(emmeans::emmeans(shannon_glm, specs = pairwise ~ Treatment, 
                                                         data = table_alpha, rg.limit = 3e5)$emmeans)
shannon_Treatment_contrasts <- data.frame(emmeans::emmeans(shannon_glm, specs = pairwise ~ Treatment, 
                                                           data = table_alpha, rg.limit = 3e5)$contrasts)

shannon_Treat_Time_emmeans <- data.frame(emmeans::emmeans(shannon_glm, specs = pairwise ~ Treatment + Timepoint, 
                                                          data = table_alpha, rg.limit = 3e5)$emmeans)
shannon_Treat_Time_contrasts <- data.frame(emmeans::emmeans(shannon_glm, specs = pairwise ~ Treatment + Timepoint, 
                                                            data = table_alpha, rg.limit = 3e5)$contrasts)

write.table(shannon_Treat_Time_emmeans, file = "/home/alexarredondo/Documents/2025_Nitrats_v6/shannon_emmeans.csv", sep = ";")
write.table(shannon_Treat_Time_contrasts, file = "/home/alexarredondo/Documents/2025_Nitrats_v6/shannon_contrasts.csv", sep = ";")

simpson <- lm(Simpson ~ Treatment*Timepoint + ID + Data_extraccio + 
                lot_kit_extraccio + extractor + tanda_pcr_16s + tanda_rentat_1 +
                tanda_pcr_index + tanda_rentat_2 + Reads_scaled, data = table_alpha)
shapiro.test(residuals(simpson)) 

simpson_glm <- glmer(Simpson ~ Treatment + Timepoint + (1|ID) + 
                       Reads_scaled + lot_kit_extraccio,
                     family = Gamma(link = "log"), data = table_alpha) 

summary(simpson_glm)

simpson_Timepoint_emmeans <- data.frame(emmeans::emmeans(simpson_glm, specs = pairwise ~ Timepoint, 
                                                         data = table_alpha, rg.limit = 3e5)$emmeans)
simpson_Timepoint_contrasts <- data.frame(emmeans::emmeans(simpson_glm, specs = pairwise ~ Timepoint, 
                                                           data = table_alpha, rg.limit = 3e5, nesting = NULL)$contrasts)

simpson_Treatment_emmeans <- data.frame(emmeans::emmeans(simpson_glm, specs = pairwise ~ Treatment, 
                                                         data = table_alpha, rg.limit = 3e5)$emmeans)
simpson_Treatment_contrasts <- data.frame(emmeans::emmeans(simpson_glm, specs = pairwise ~ Treatment, 
                                                           data = table_alpha, rg.limit = 3e5)$contrasts)

simpson_Treat_Time_emmeans <- data.frame(emmeans::emmeans(simpson_glm, specs = pairwise ~ Treatment + Timepoint, 
                                                          data = table_alpha, rg.limit = 3e5)$emmeans)
simpson_Treat_Time_contrasts <- data.frame(emmeans::emmeans(simpson_glm, specs = pairwise ~ Treatment + Timepoint, 
                                                            data = table_alpha, rg.limit = 3e5)$contrasts)

write.table(simpson_Treat_Time_emmeans, file = "/home/alexarredondo/Documents/2025_Nitrats_v6/simpson_emmeans.csv", sep = ";")
write.table(simpson_Treat_Time_contrasts, file = "/home/alexarredondo/Documents/2025_Nitrats_v6/simpson_contrasts.csv", sep = ";")
```

####################################### Beta diversity #######################################

# We are going to use aitchison distances and start by visualizing the plots.

# Treatment and timepoint
ps_genera_v1_pruned %>%
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "Treatment", shape = "Timepoint", size = 3, auto_caption = NA) +
  scale_colour_brewer(palette = "Dark2") +
  stat_ellipse(aes(linetype = Timepoint, colour = Treatment), linewidth = 0.3) +
  ggtitle("Treatment - Timepoint") +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

# Extraction date
ps_genera_v1_pruned %>%
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "Data_extraccio", shape = "Data_extraccio", size = 3, auto_caption = NA) +
  scale_colour_brewer(palette = "Dark2") +
  stat_ellipse(aes(linetype = Data_extraccio, colour = Data_extraccio), linewidth = 0.3) +
  ggtitle("Data extracció") +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

# Extraction kit batch
ps_genera_v1_pruned %>%
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "lot_kit_extraccio", shape = "lot_kit_extraccio", size = 3, auto_caption = NA) +
  scale_colour_brewer(palette = "Dark2") +
  stat_ellipse(aes(linetype = lot_kit_extraccio, colour = lot_kit_extraccio), linewidth = 0.3) +
  ggtitle("Lot kit d'extracció") +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

# N reads
set.seed(123) 
quartiles <- quantile(sample_data(ps_genera_v1_pruned)$Reads, probs = seq(0, 1, 0.25)) # Transform continous data to categoric using quartiles.
sample_data(ps_genera_v1_pruned)$reads_quartiles <- cut(sample_data(ps_genera_v1_pruned)$Reads, breaks = quartiles, include.lowest = TRUE, 
                                                        labels = c("Q1","Q2","Q3","Q4"))

ps_genera_v1_pruned %>%
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "reads_quartiles", shape = "reads_quartiles", size = 3, auto_caption = NA) +
  scale_colour_brewer(palette = "Dark2") +
  stat_ellipse(aes(linetype = reads_quartiles, colour = reads_quartiles), linewidth = 0.3) +
  ggtitle("Quartils nº de reads") +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

# Now let's check if we can find statistical differences. We are going to use the same distances (aitchison) and use the pairwise.adonis function to assess pairwise comparisons for each variable.
ps_dist_matrix <- dist_calc(ps_genera_v1_pruned, "aitchison")
dist_get(ps_dist_matrix)
vegan::adonis2(dist_get(ps_dist_matrix) ~ paste(sample_data(ps_genera_v1_pruned)$Treatment, sample_data(ps_genera_v1_pruned)$Timepoint, sep = " "))
pairwise.adonis(dist_get(ps_dist_matrix), paste(sample_data(ps_genera_v1_pruned)$Treatment, sample_data(ps_genera_v1_pruned)$Timepoint, sep = " "))
pairwise.adonis(dist_get(ps_dist_matrix), sample_data(ps_genera_v1_pruned)$Data_extraccio)
pairwise.adonis(dist_get(ps_dist_matrix), sample_data(ps_genera_v1_pruned)$lot_kit_extraccio)
pairwise.adonis(dist_get(ps_dist_matrix), sample_data(ps_genera_v1_pruned)$ID)
pairwise.adonis(dist_get(ps_dist_matrix), sample_data(ps_genera_v1_pruned)$reads_quartiles)

####################################### Relative abundance #######################################

# This two functions will transform the otu table, which is in absolute numbers, to relative numbers and then we create a table with a few descriptive statistics. 
funcio_taula_relabund <- function(ps_object){
  taula_otus <- data.frame(otu_table(ps_object))
  taula_otus_rel <- taula_otus
  taula_taxa <- data.frame(tax_table(ps_object))
  colnames(taula_otus_rel) <- taula_taxa$Genus
  for(i in 1:ncol(taula_otus)){
    for(k in 1:nrow(taula_otus)){
      taula_otus_rel[k,i] <- (taula_otus[k,i]/sum(taula_otus[k,]))*100
    }
  }
  return(taula_otus_rel)
}

funcio_top10_relabund <- function(taula_rel){
  vector_sum <- c()
  vector_sd <- c()
  for(i in 1:ncol(taula_rel)){
    vector_sum[i] <- mean(taula_rel[,i])
    vector_sd[i] <- sd(taula_rel[,i])
  }
  taula <- data.frame(Taxa = colnames(taula_rel),
                      Mean = vector_sum,
                      SD = vector_sd)
  taula <- taula[order(taula$Mean, decreasing = TRUE),]
  return(taula)
}

# First we apply the function to transform to relative abundance.
ps_rel <- funcio_taula_relabund(ps_glom_genus_v1_pruned)

# Subset phyloseq objects based on the name of the samples, which include the variables that we want to study.
ps_rel_Basal <- ps_rel[grep("Saliva", rownames(ps_rel)),]
ps_rel_7d_Base <- ps_rel[grep("7d-Base", rownames(ps_rel)),]
ps_rel_14d_Base <- ps_rel[grep("Fi-Base", rownames(ps_rel)),]
ps_rel_7d_CDQ <- ps_rel[grep("7d-CDQ", rownames(ps_rel)),]
ps_rel_14d_CDQ <- ps_rel[grep("Fi-CDQ", rownames(ps_rel)),]
ps_rel_7d_KNO3 <- ps_rel[grep("7d-KNO3", rownames(ps_rel)),]
ps_rel_14d_KNO3 <- ps_rel[grep("Fi-KNO3", rownames(ps_rel)),]
ps_rel_7d_MQ <- ps_rel[grep("7d-MQ", rownames(ps_rel)),]
ps_rel_14d_MQ <- ps_rel[grep("Fi-MQ", rownames(ps_rel)),]

# Extract the top10 with the higher mean from the general population.
ps_rel_top10 <- funcio_top10_relabund(ps_rel)[1:10,]

# Now we extract the ones with the higher mean from each population that match the ones from the general population. This is an optional step, and data can be retrieved from each population separately. This has been done to compare the same genera in each population.
top10_rel_Basal <- funcio_top10_relabund(select(ps_rel_Basal, ps_rel_top10$Taxa))
top10_rel_7d_Base <- funcio_top10_relabund(select(ps_rel_7d_Base, ps_rel_top10$Taxa))
top10_rel_14d_Base <- funcio_top10_relabund(select(ps_rel_14d_Base, ps_rel_top10$Taxa))
top10_rel_7d_CDQ <- funcio_top10_relabund(select(ps_rel_7d_CDQ, ps_rel_top10$Taxa))
top10_rel_14d_CDQ <- funcio_top10_relabund(select(ps_rel_14d_CDQ, ps_rel_top10$Taxa)) 
top10_rel_7d_KNO3 <- funcio_top10_relabund(select(ps_rel_7d_KNO3, ps_rel_top10$Taxa))  
top10_rel_14d_KNO3 <- funcio_top10_relabund(select(ps_rel_14d_KNO3, ps_rel_top10$Taxa))
top10_rel_7d_MQ <- funcio_top10_relabund(select(ps_rel_7d_MQ, ps_rel_top10$Taxa))
top10_rel_14d_MQ <- funcio_top10_relabund(select(ps_rel_14d_MQ, ps_rel_top10$Taxa))

# And we add the element "other", wich contains all the genera not described in the top10 of the general population. 
top10_rel_Basal <- rbind(top10_rel_Basal, c("Other", 100-sum(as.numeric(top10_rel_Basal$Mean)), NA))
top10_rel_7d_Base <- rbind(top10_rel_7d_Base, c("Other", 100-sum(as.numeric(top10_rel_7d_Base$Mean)), NA))
top10_rel_14d_Base <- rbind(top10_rel_14d_Base, c("Other", 100-sum(as.numeric(top10_rel_14d_Base$Mean)), NA))
top10_rel_7d_CDQ <- rbind(top10_rel_7d_CDQ, c("Other", 100-sum(as.numeric(top10_rel_7d_CDQ$Mean)), NA))
top10_rel_14d_CDQ <- rbind(top10_rel_14d_CDQ, c("Other", 100-sum(as.numeric(top10_rel_14d_CDQ$Mean)), NA))
top10_rel_7d_KNO3 <- rbind(top10_rel_7d_KNO3, c("Other", 100-sum(as.numeric(top10_rel_7d_KNO3$Mean)), NA)) 
top10_rel_14d_KNO3 <- rbind(top10_rel_14d_KNO3, c("Other", 100-sum(as.numeric(top10_rel_14d_KNO3$Mean)), NA))
top10_rel_7d_MQ <- rbind(top10_rel_7d_MQ, c("Other", 100-sum(as.numeric(top10_rel_7d_MQ$Mean)), NA))
top10_rel_14d_MQ <- rbind(top10_rel_14d_MQ, c("Other", 100-sum(as.numeric(top10_rel_14d_MQ$Mean)), NA))

# Merge all the data to prepare a single plot.
ps_global <- rbind(top10_rel_Basal,
                   top10_rel_7d_Base,
                   top10_rel_14d_Base,
                   top10_rel_7d_CDQ,
                   top10_rel_14d_CDQ,
                   top10_rel_7d_KNO3,
                   top10_rel_14d_KNO3,
                   top10_rel_7d_MQ,
                   top10_rel_14d_MQ)

ps_global <- cbind(ps_global, Comparison = c(rep("Basal", 11), 
                                             rep("7d - Base", 11), 
                                             rep("14d - Base", 11), 
                                             rep("7d - CDQ", 11), 
                                             rep("14d - CDQ", 11), 
                                             rep("7d - KNO3", 11),
                                             rep("14d - KNO3", 11),
                                             rep("7d - MQ", 11),
                                             rep("14d - MQ", 11))) # add a new variable for ggplot.

ps_global$Comparison <- factor(ps_global$Comparison, levels = c("Basal", "7d - MQ", "14d - MQ", "7d - KNO3",
                                                                "14d - KNO3", "7d - Base", "14d - Base", "7d - CDQ",
                                                                "14d - CDQ")) # factorize it to order the x-axis

global_graph <- ggplot(ps_global, aes(x = Comparison, y = as.numeric(Mean), fill = Taxa)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_brewer(palette = "Spectral") + 
  ggtitle("Abundància relativa top10 gèneres") +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) +
  labs(y = "% Abundància relativa") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

# Relative differential abundance. We prepare a function to compare between timepoints within the same treatment and between treatments within the same timepoint. 
set.seed(123)

# First we start comparing timepoints within the same treatment.  
funcio_timepoint <- function(ps_object){
  output <- ancombc2(data = ps_object, meta_data = sample_data(ps_object),
                     fix_formula = "Timepoint", rand_formula = "(1|ID)",
                     p_adj_method = "holm", pseudo_sens = TRUE,
                     prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                     group = "Timepoint", struc_zero = FALSE, neg_lb = FALSE,
                     alpha = 0.05, n_cl = 2, verbose = TRUE,
                     global = FALSE, pairwise = TRUE, 
                     dunnet = FALSE, trend = FALSE,
                     iter_control = list(tol = 1e-5, max_iter = 20, 
                                         verbose = FALSE),
                     em_control = list(tol = 1e-5, max_iter = 100),
                     mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                     trend_control = NULL)
  resultats <- output$res
  return(resultats)
}

taxa_names(ps_glom_genus_v1_pruned) <- c(tax_table(ps_glom_genus_v1_pruned)[,6])
diff_MQ <- subset_samples(ps_glom_genus_v1_pruned, Treatment == "MQ")
res_MQ <- funcio_timepoint(diff_MQ)

diff_Base <- subset_samples(ps_glom_genus_v1_pruned, Treatment == "Base")
res_Base <- funcio_timepoint(diff_Base)

diff_KNO3 <- subset_samples(ps_glom_genus_v1_pruned, Treatment == "KNO3")
res_KNO3 <- funcio_timepoint(diff_KNO3)

diff_CDQ <- subset_samples(ps_glom_genus_v1_pruned, Treatment == "CDQ")
res_CDQ <- funcio_timepoint(diff_CDQ)

# Prepare simplified tables with the results. 
taula_res_MQ <- data.frame(Taxa = res_MQ$taxon,
                           LFC= res_MQ$lfc_Timepoint7d,
                           p_value = res_MQ$p_Timepoint7d,
                           Comparison = rep("7d vs 14d - MQ", nrow(res_MQ)))

taula_res_Base <- data.frame(Taxa = res_Base$taxon,
                             LFC= res_Base$lfc_Timepoint7d,
                             p_value = res_Base$p_Timepoint7d,
                             Comparison = rep("7d vs 14d - Base", nrow(res_Base)))

taula_res_KNO3 <- data.frame(Taxa = res_KNO3$taxon,
                             LFC= res_KNO3$lfc_Timepoint7d,
                             p_value = res_KNO3$p_Timepoint7d,
                             Comparison = rep("7d vs 14d - KNO3", nrow(res_KNO3)))

taula_res_CDQ <- data.frame(Taxa = res_CDQ$taxon,
                            LFC= res_CDQ$lfc_Timepoint7d,
                            p_value = res_CDQ$p_Timepoint7d,
                            Comparison = rep("7d vs 14d - CDQ", nrow(res_CDQ)))

# Filter results by p value.
taula_res_MQ_p <- subset(taula_res_MQ, taula_res_MQ$p_value < 0.05)
taula_res_Base_p <- subset(taula_res_Base, taula_res_Base$p_value < 0.05)
taula_res_KNO3_p <- subset(taula_res_KNO3, taula_res_KNO3$p_value < 0.05)
taula_res_CDQ_p <- subset(taula_res_CDQ, taula_res_CDQ$p_value < 0.05)

taula_timepoints <- rbind(taula_res_MQ_p, taula_res_Base_p, taula_res_KNO3_p, taula_res_CDQ_p)
taula_timepoints$Comparison <- factor(taula_timepoints$Comparison, levels = c("7d vs 14d - MQ", "7d vs 14d - KNO3",
                                                                              "7d vs 14d - Base", "7d vs 14d - CDQ"))

# Prepare some heatmaps to visualize results.
lo_intra = floor(min(taula_timepoints$LFC))
up_intra = ceiling(max(taula_timepoints$LFC))
mid = 0
taula_timepoints %>%
  ggplot(aes(x = Comparison, y = Taxa, fill = LFC)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo_inter, up_inter),
                       name = NULL) +
  geom_text(aes(Comparison, Taxa, label = round(LFC,2)), size = 4) +
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = "Log fold change dels diferents tractaments \nen els seus timepoints - p-value < 0.05") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_discrete(limits=rev)

# Now we do the same but comparing treatments within the same timepoints.  
funcio_treatment <- function(ps_object){
  output <- ancombc2(data = ps_object, meta_data = sample_data(ps_object),
                     fix_formula = "Treatment", rand_formula = "(1|ID)",
                     p_adj_method = "holm", pseudo_sens = TRUE,
                     prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                     group = "Treatment", struc_zero = FALSE, neg_lb = FALSE,
                     alpha = 0.05, n_cl = 2, verbose = TRUE,
                     global = FALSE, pairwise = TRUE, 
                     dunnet = FALSE, trend = FALSE,
                     iter_control = list(tol = 1e-5, max_iter = 20, 
                                         verbose = FALSE),
                     em_control = list(tol = 1e-5, max_iter = 100),
                     mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                     trend_control = NULL)
  resultats <- output$res
  return(resultats)
}

diff_MQvsKNO3 <- subset_samples(ps_glom_genus_v1_pruned, Treatment == "MQ" | Treatment == "KNO3")
diff_MQvsKNO3_7d <- subset_samples(diff_MQvsKNO3, Timepoint == "7d")
diff_MQvsKNO3_14d <- subset_samples(diff_MQvsKNO3, Timepoint == "14d")

res_MQvsKNO3_7d <- funcio_treatment(diff_MQvsKNO3_7d)
res_MQvsKNO3_14d <- funcio_treatment(diff_MQvsKNO3_14d)

diff_BasevsCDQ <- subset_samples(ps_glom_genus_v1_pruned, Treatment == "Base" | Treatment == "CDQ")
diff_BasevsCDQ_7d <- subset_samples(diff_BasevsCDQ, Timepoint == "7d")
diff_BasevsCDQ_14d <- subset_samples(diff_BasevsCDQ, Timepoint == "14d")

res_BasevsCDQ_7d <- funcio_treatment(diff_BasevsCDQ_7d)
res_BasevsCDQ_14d <- funcio_treatment(diff_BasevsCDQ_14d)

diff_BasevsMQ <- subset_samples(ps_glom_genus_v1_pruned, Treatment == "Base" | Treatment == "MQ")
diff_BasevsMQ_7d <- subset_samples(diff_BasevsMQ, Timepoint == "7d")
diff_BasevsMQ_14d <- subset_samples(diff_BasevsMQ, Timepoint == "14d")

res_BasevsMQ_7d <- funcio_treatment(diff_BasevsMQ_7d)
res_BasevsMQ_14d <- funcio_treatment(diff_BasevsMQ_14d)

diff_KNO3vsCDQ <- subset_samples(ps_glom_genus_v1_pruned, Treatment == "KNO3" | Treatment == "CDQ")
diff_KNO3vsCDQ_7d <- subset_samples(diff_KNO3vsCDQ, Timepoint == "7d")
diff_KNO3vsCDQ_14d <- subset_samples(diff_KNO3vsCDQ, Timepoint == "14d")

res_KNO3vsCDQ_7d <- funcio_treatment(diff_KNO3vsCDQ_7d)
res_KNO3vsCDQ_14d <- funcio_treatment(diff_KNO3vsCDQ_14d)

diff_MQvsCDQ <- subset_samples(ps_glom_genus_v1_pruned, Treatment == "MQ" | Treatment == "CDQ")
diff_MQvsCDQ_7d <- subset_samples(diff_MQvsCDQ, Timepoint == "7d")
diff_MQvsCDQ_14d <- subset_samples(diff_MQvsCDQ, Timepoint == "14d")

res_MQvsCDQ_7d <- funcio_treatment(diff_MQvsCDQ_7d)
res_MQvsCDQ_14d <- funcio_treatment(diff_MQvsCDQ_14d)

taula_res_MQvsKNO3_7d <- data.frame(Taxa = res_MQvsKNO3_7d$taxon,
                                    LFC= res_MQvsKNO3_7d$lfc_TreatmentMQ,
                                    p_value = res_MQvsKNO3_7d$p_TreatmentMQ,
                                    Comparison = rep("MQ vs KNO3 - 7d", nrow(res_MQvsKNO3_7d)))

taula_res_MQvsKNO3_14d <- data.frame(Taxa = res_MQvsKNO3_14d$taxon,
                                     LFC= res_MQvsKNO3_14d$lfc_TreatmentMQ,
                                     p_value = res_MQvsKNO3_14d$p_TreatmentMQ,
                                     Comparison = rep("MQ vs KNO3 - 14d", nrow(res_MQvsKNO3_14d)))

taula_res_BasevsCDQ_7d <- data.frame(Taxa = res_BasevsCDQ_7d$taxon,
                                     LFC= res_BasevsCDQ_7d$lfc_TreatmentCDQ,
                                     p_value = res_BasevsCDQ_7d$p_TreatmentCDQ,
                                     Comparison = rep("CDQ vs Base - 7d", nrow(res_BasevsCDQ_7d)))

taula_res_BasevsCDQ_14d <- data.frame(Taxa = res_BasevsCDQ_14d$taxon,
                                      LFC= res_BasevsCDQ_14d$lfc_TreatmentCDQ,
                                      p_value = res_BasevsCDQ_14d$p_TreatmentCDQ,
                                      Comparison = rep("CDQ vs Base - 14d", nrow(res_BasevsCDQ_14d)))

taula_res_BasevsMQ_7d <- data.frame(Taxa = res_BasevsMQ_7d$taxon,
                                    LFC= res_BasevsMQ_7d$lfc_TreatmentMQ,
                                    p_value = res_BasevsMQ_7d$p_TreatmentMQ,
                                    Comparison = rep("MQ vs Base - 7d", nrow(res_BasevsMQ_7d)))

taula_res_BasevsMQ_14d <- data.frame(Taxa = res_BasevsMQ_14d$taxon,
                                     LFC= res_BasevsMQ_14d$lfc_TreatmentMQ,
                                     p_value = res_BasevsMQ_14d$p_TreatmentMQ,
                                     Comparison = rep("MQ vs Base - 14d", nrow(res_BasevsMQ_14d)))

taula_res_KNO3vsCDQ_7d <- data.frame(Taxa = res_KNO3vsCDQ_7d$taxon,
                                     LFC= res_KNO3vsCDQ_7d$lfc_TreatmentKNO3,
                                     p_value = res_KNO3vsCDQ_7d$p_TreatmentKNO3,
                                     Comparison = rep("KNO3 vs CDQ - 7d", nrow(res_KNO3vsCDQ_7d)))

taula_res_KNO3vsCDQ_14d <- data.frame(Taxa = res_KNO3vsCDQ_14d$taxon,
                                      LFC= res_KNO3vsCDQ_14d$lfc_TreatmentKNO3,
                                      p_value = res_KNO3vsCDQ_14d$p_TreatmentKNO3,
                                      Comparison = rep("KNO3 vs CDQ - 14d", nrow(res_KNO3vsCDQ_14d)))

taula_res_MQvsCDQ_7d <- data.frame(Taxa = res_MQvsCDQ_7d$taxon,
                                   LFC= res_MQvsCDQ_7d$lfc_TreatmentMQ,
                                   p_value = res_MQvsCDQ_7d$p_TreatmentMQ,
                                   Comparison = rep("MQ vs CDQ - 7d", nrow(res_MQvsCDQ_7d)))

taula_res_MQvsCDQ_14d <- data.frame(Taxa = res_MQvsCDQ_14d$taxon,
                                    LFC= res_MQvsCDQ_14d$lfc_TreatmentMQ,
                                    p_value = res_MQvsCDQ_14d$p_TreatmentMQ,
                                    Comparison = rep("MQ vs CDQ - 14d", nrow(res_MQvsCDQ_14d)))

taula_res_MQvsKNO3_7d_p <- subset(taula_res_MQvsKNO3_7d, taula_res_MQvsKNO3_7d$p_value < 0.05)
taula_res_MQvsKNO3_14d_p <- subset(taula_res_MQvsKNO3_14d, taula_res_MQvsKNO3_14d$p_value < 0.05)
taula_res_BasevsCDQ_7d_p <- subset(taula_res_BasevsCDQ_7d, taula_res_BasevsCDQ_7d$p_value < 0.05)
taula_res_BasevsCDQ_14d_p <- subset(taula_res_BasevsCDQ_14d, taula_res_BasevsCDQ_14d$p_value < 0.05)
taula_res_BasevsMQ_7d_p <- subset(taula_res_BasevsMQ_7d, taula_res_BasevsMQ_7d$p_value < 0.05)
taula_res_BasevsMQ_14d_p <- subset(taula_res_BasevsMQ_14d, taula_res_BasevsMQ_14d$p_value < 0.05)
taula_res_KNO3vsCDQ_7d_p <- subset(taula_res_KNO3vsCDQ_7d, taula_res_KNO3vsCDQ_7d$p_value < 0.05)
taula_res_KNO3vsCDQ_14d_p <- subset(taula_res_KNO3vsCDQ_14d, taula_res_KNO3vsCDQ_14d$p_value < 0.05)
taula_res_MQvsCDQ_7d_p <- subset(taula_res_MQvsCDQ_7d, taula_res_MQvsCDQ_7d$p_value < 0.05)
taula_res_MQvsCDQ_14d_p <- subset(taula_res_MQvsCDQ_14d, taula_res_MQvsCDQ_14d$p_value < 0.05)

taula_treatments <- rbind(taula_res_MQvsKNO3_7d_p, taula_res_MQvsKNO3_14d_p, 
                          taula_res_BasevsCDQ_7d_p, taula_res_BasevsCDQ_14d_p)

taula_treatments$Comparison <- factor(taula_treatments$Comparison, levels = c("MQ vs KNO3 - 7d", "MQ vs KNO3 - 14d",
                                                                              "CDQ vs Base - 7d", "CDQ vs Base - 14d"))

lo_inter = floor(min(taula_treatments$LFC))
up_inter = ceiling(max(taula_treatments$LFC))
mid = 0
taula_treatments %>%
  ggplot(aes(x = Comparison, y = Taxa, fill = LFC)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo_inter, up_inter),
                       name = NULL) +
  geom_text(aes(Comparison, Taxa, label = round(LFC,2)), size = 4) +
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = "Log fold change dels diferents tractaments \nseparats per timepoints - p-value < 0.05") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_discrete(limits=rev)


# We repeat the ancombc2 analysis with the object with spp level and withou aglommeration of taxa. The code is largely the same as before.
funcio_treatment <- function(ps_object){
  output <- ancombc2(data = ps_object, meta_data = sample_data(ps_object),
                     fix_formula = "Treatment", rand_formula = "(1|ID)",
                     p_adj_method = "holm", pseudo_sens = TRUE,
                     prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                     group = "Treatment", struc_zero = FALSE, neg_lb = FALSE,
                     alpha = 0.05, n_cl = 2, verbose = TRUE,
                     global = FALSE, pairwise = TRUE, 
                     dunnet = FALSE, trend = FALSE,
                     iter_control = list(tol = 1e-5, max_iter = 20, 
                                         verbose = FALSE),
                     em_control = list(tol = 1e-5, max_iter = 100),
                     mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                     trend_control = NULL)
  resultats <- output$res
  return(resultats)
}

vector_names <- paste(data.frame(tax_table(ps_spp_v1_pruned))$Genus, data.frame(tax_table(ps_spp_v1_pruned))$Species, sep = "_")
vector_names <- make.unique(vector_names, sep = "_")

taxa_names(ps_spp_v1_pruned) <- vector_names

spp_diff_MQvsKNO3 <- subset_samples(ps_spp_v1_pruned, Treatment == "MQ" | Treatment == "KNO3")
spp_diff_MQvsKNO3_7d <- subset_samples(spp_diff_MQvsKNO3, Timepoint == "7d")
spp_diff_MQvsKNO3_14d <- subset_samples(spp_diff_MQvsKNO3, Timepoint == "14d")

spp_res_MQvsKNO3_7d <- funcio_treatment(spp_diff_MQvsKNO3_7d)
spp_res_MQvsKNO3_14d <- funcio_treatment(spp_diff_MQvsKNO3_14d)

spp_diff_BasevsCDQ <- subset_samples(ps_spp_v1_pruned, Treatment == "Base" | Treatment == "CDQ")
spp_diff_BasevsCDQ_7d <- subset_samples(spp_diff_BasevsCDQ, Timepoint == "7d")
spp_diff_BasevsCDQ_14d <- subset_samples(spp_diff_BasevsCDQ, Timepoint == "14d")

spp_res_BasevsCDQ_7d <- funcio_treatment(spp_diff_BasevsCDQ_7d)
spp_res_BasevsCDQ_14d <- funcio_treatment(spp_diff_BasevsCDQ_14d)

spp_taula_res_MQvsKNO3_7d <- data.frame(Taxa = spp_res_MQvsKNO3_7d$taxon,
                                        LFC= spp_res_MQvsKNO3_7d$lfc_TreatmentMQ,
                                        p_value = spp_res_MQvsKNO3_7d$p_TreatmentMQ,
                                        Comparison = rep("MQ vs KNO3 - 7d", nrow(spp_res_MQvsKNO3_7d)))

spp_taula_res_MQvsKNO3_14d <- data.frame(Taxa = spp_res_MQvsKNO3_14d$taxon,
                                         LFC= spp_res_MQvsKNO3_14d$lfc_TreatmentMQ,
                                         p_value = spp_res_MQvsKNO3_14d$p_TreatmentMQ,
                                         Comparison = rep("MQ vs KNO3 - 14d", nrow(spp_res_MQvsKNO3_14d)))

spp_taula_res_BasevsCDQ_7d <- data.frame(Taxa = spp_res_BasevsCDQ_7d$taxon,
                                         LFC= spp_res_BasevsCDQ_7d$lfc_TreatmentCDQ,
                                         p_value = spp_res_BasevsCDQ_7d$p_TreatmentCDQ,
                                         Comparison = rep("CDQ vs Base - 7d", nrow(spp_res_BasevsCDQ_7d)))

spp_taula_res_BasevsCDQ_14d <- data.frame(Taxa = spp_res_BasevsCDQ_14d$taxon,
                                          LFC= spp_res_BasevsCDQ_14d$lfc_TreatmentCDQ,
                                          p_value = spp_res_BasevsCDQ_14d$p_TreatmentCDQ,
                                          Comparison = rep("CDQ vs Base - 14d", nrow(spp_res_BasevsCDQ_14d)))

spp_taula_res_MQvsKNO3_7d_p <- subset(spp_taula_res_MQvsKNO3_7d, spp_taula_res_MQvsKNO3_7d$p_value < 0.05)
spp_taula_res_MQvsKNO3_14d_p <- subset(spp_taula_res_MQvsKNO3_14d, spp_taula_res_MQvsKNO3_14d$p_value < 0.05)
spp_taula_res_BasevsCDQ_7d_p <- subset(spp_taula_res_BasevsCDQ_7d, spp_taula_res_BasevsCDQ_7d$p_value < 0.05)
spp_taula_res_BasevsCDQ_14d_p <- subset(spp_taula_res_BasevsCDQ_14d, spp_taula_res_BasevsCDQ_14d$p_value < 0.05)

spp_taula_treatments <- rbind(spp_taula_res_MQvsKNO3_7d_p, spp_taula_res_MQvsKNO3_14d_p, 
                              spp_taula_res_BasevsCDQ_7d_p, spp_taula_res_BasevsCDQ_14d_p)

spp_taula_treatments$Comparison <- factor(spp_taula_treatments$Comparison, levels = c("MQ vs KNO3 - 7d", "MQ vs KNO3 - 14d",
                                                                                      "CDQ vs Base - 7d", "CDQ vs Base - 14d"))

spp_lo_inter = floor(min(spp_taula_treatments$LFC))
spp_up_inter = ceiling(max(spp_taula_treatments$LFC))
mid = 0
grafic <- spp_taula_treatments %>%
  ggplot(aes(x = Comparison, y = Taxa, fill = LFC)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(spp_lo_inter, spp_up_inter),
                       name = NULL) +
  geom_text(aes(Comparison, Taxa, label = round(LFC,2)), size = 4) +
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = "Log fold change dels diferents tractaments \nseparats per timepoints - p-value < 0.05") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_discrete(limits=rev) 

set.seed(123)

funcio_timepoint <- function(ps_object){
  output <- ancombc2(data = ps_object, meta_data = sample_data(ps_object),
                     fix_formula = "Timepoint", rand_formula = "(1|ID)",
                     p_adj_method = "holm", pseudo_sens = TRUE,
                     prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                     group = "Timepoint", struc_zero = FALSE, neg_lb = FALSE,
                     alpha = 0.05, n_cl = 2, verbose = TRUE,
                     global = FALSE, pairwise = TRUE, 
                     dunnet = FALSE, trend = FALSE,
                     iter_control = list(tol = 1e-5, max_iter = 20, 
                                         verbose = FALSE),
                     em_control = list(tol = 1e-5, max_iter = 100),
                     mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                     trend_control = NULL)
  resultats <- output$res
  return(resultats)
}

spp_diff_MQ <- subset_samples(ps_spp_v1_pruned, Treatment == "MQ")
spp_res_MQ <- funcio_timepoint(spp_diff_MQ)

spp_diff_Base <- subset_samples(ps_spp_v1_pruned, Treatment == "Base")
spp_res_Base <- funcio_timepoint(spp_diff_Base)

spp_diff_KNO3 <- subset_samples(ps_spp_v1_pruned, Treatment == "KNO3")
spp_res_KNO3 <- funcio_timepoint(spp_diff_KNO3)

spp_diff_CDQ <- subset_samples(ps_spp_v1_pruned, Treatment == "CDQ")
spp_res_CDQ <- funcio_timepoint(spp_diff_CDQ)

spp_taula_res_MQ <- data.frame(Taxa = spp_res_MQ$taxon,
                               LFC= spp_res_MQ$lfc_Timepoint7d,
                               p_value = spp_res_MQ$p_Timepoint7d,
                               Comparison = rep("7d vs 14d - MQ", nrow(spp_res_MQ)))

spp_taula_res_Base <- data.frame(Taxa = spp_res_Base$taxon,
                                 LFC= spp_res_Base$lfc_Timepoint7d,
                                 p_value = spp_res_Base$p_Timepoint7d,
                                 Comparison = rep("7d vs 14d - Base", nrow(spp_res_Base)))

spp_taula_res_KNO3 <- data.frame(Taxa = spp_res_KNO3$taxon,
                                 LFC= spp_res_KNO3$lfc_Timepoint7d,
                                 p_value = spp_res_KNO3$p_Timepoint7d,
                                 Comparison = rep("7d vs 14d - KNO3", nrow(spp_res_KNO3)))

spp_taula_res_CDQ <- data.frame(Taxa = spp_res_CDQ$taxon,
                                LFC= spp_res_CDQ$lfc_Timepoint7d,
                                p_value = spp_res_CDQ$p_Timepoint7d,
                                Comparison = rep("7d vs 14d - CDQ", nrow(spp_res_CDQ)))

spp_taula_res_MQ_p <- subset(spp_taula_res_MQ, spp_taula_res_MQ$p_value < 0.05)
spp_taula_res_Base_p <- subset(spp_taula_res_Base, spp_taula_res_Base$p_value < 0.05)
spp_taula_res_KNO3_p <- subset(spp_taula_res_KNO3, spp_taula_res_KNO3$p_value < 0.05)
spp_taula_res_CDQ_p <- subset(spp_taula_res_CDQ, spp_taula_res_CDQ$p_value < 0.05)

spp_taula_timepoints <- rbind(spp_taula_res_MQ_p, spp_taula_res_Base_p, spp_taula_res_KNO3_p, spp_taula_res_CDQ_p)
spp_taula_timepoints$Comparison <- factor(spp_taula_timepoints$Comparison, levels = c("7d vs 14d - MQ", "7d vs 14d - KNO3",
                                                                                      "7d vs 14d - Base", "7d vs 14d - CDQ"))

spp_lo_intra = floor(min(spp_taula_timepoints$LFC))
spp_up_intra = ceiling(max(spp_taula_timepoints$LFC))
mid = 0
grafic2 <- spp_taula_timepoints %>%
  ggplot(aes(x = Comparison, y = Taxa, fill = LFC)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo_inter, up_inter),
                       name = NULL) +
  geom_text(aes(Comparison, Taxa, label = round(LFC,2)), size = 4) +
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = "Log fold change dels diferents tractaments \nen els seus timepoints - p-value < 0.05") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_discrete(limits=rev)

####################################### Functional profiles #######################################

# We use Picrust2 and the tutorial obtained from https://forum.qiime2.org/t/importing-dada2-and-phyloseq-objects-to-qiime-2/4683 that allows to use data from a phyloseq object. 

# Export taxonomy table.
tax <- as(tax_table(ps_genera_v1_pruned),"matrix")
tax_cols <- colnames(tax)
tax <- as.data.frame(tax)
tax$taxonomy <- do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co] <- NULL
write.table(tax, "/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/tax.txt", quote=FALSE, col.names=FALSE, sep="\t")

# Export the ASV table as a biom file.
otu <- t(as(otu_table(ps_genera_v1_pruned),"matrix")) # 't' to transform if taxa_are_rows=FALSE
otu_biom <- make_biom(data=otu)
write_biom(otu_biom,"/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/otu_biom.biom")

# Export metadata
write.table(sample_data(ps_genera_v1_pruned),
            "/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/sample-metadata.txt", 
            sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# Extract the .fasta file.
dna <- DNAStringSet(colnames(seqtab.nochim))
names(dna) <- paste0("ASV", seq_len(length(dna)))
writeXStringSet(dna, "/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/study_seqs.fna")

# The order in picrust2 has been this one:
./picrust2_pipeline.py -s /home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/study_seqs.fna -i /home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/otu_biom.biom -o /home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/output -p 100

# Import the file back to R and prepare two jerarchy groupings. First the top level and then the second level.

# Top level.
abun <- read.table("/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/output/pathways_out/path_abun_unstrat.tsv", 
                   header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
meta_top <- read.table("/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/metacyc_pathways_info_prokaryotes_top_level.tsv", 
                       header = TRUE, sep = "\t", quote = "", comment.char = "")
colnames(meta_top) <- c("pathway_id", "top_level_category")
merged_top <- merge(meta_top, abun, by.x = "pathway_id", by.y = "row.names", all.y = TRUE)
category_abun_top <- merged_top %>%
  group_by(top_level_category) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

# Second level.
meta_2 <- read.table("/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/metacyc_pathways_info_prokaryotes_sec_level.tsv", 
                     header = TRUE, sep = "\t", quote = "", comment.char = "")
colnames(meta_2) <- c("pathway_id", "level2_category")
merged_2 <- merge(meta_2, abun, by.x = "pathway_id", by.y = "row.names", all.y = TRUE)
category_abun_2 <- merged_2 %>%
  group_by(level2_category) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

# Now perform som statistical analysis.

# First at level 2.
keep <- rowMeans(category_abun_2 > 0) > 0.1 # Filter those pathways with less than 10%.
category_abun_2_pruned <- category_abun_2[keep, ]
meta <- read.csv("/home/alexarredondo/Documents/2025_Nitrats_v6/metadata.csv", sep = ";")

# Make sure rownames are the same. 
category_abun_2_pruned <- category_abun_2_pruned[, rownames(meta)]

noms_files <- category_abun_2_pruned$level2_category
category_abun_2_pruned <- category_abun_2_pruned[,-1]
rownames(category_abun_2_pruned) <- noms_files
rownames(meta) <- meta$Sample

# Create phyloseq object were the "ASV" are now the pathways and subset by groups of interest.
otu_nivell2 <- otu_table(as.matrix(category_abun_2_pruned), taxa_are_rows = TRUE)
sam_nivell2 <- sample_data(meta)
ps_nivell2 <- phyloseq(otu_nivell2, sam_nivell2)
ps_nivell2_7d <- subset_samples(ps_nivell2, Timepoint == "7d")
ps_nivell2_14d <- subset_samples(ps_nivell2, Timepoint == "14d")
ps_nivell2_MQ <- subset_samples(ps_nivell2, Treatment == "MQ")
ps_nivell2_CDQ <- subset_samples(ps_nivell2, Treatment == "CDQ")
ps_nivell2_KNO3 <- subset_samples(ps_nivell2, Treatment == "KNO3")
ps_nivell2_Base <- subset_samples(ps_nivell2, Treatment == "Base")

# Use ancombc2 to detect differences.
out_7d <- ancombc2(
  data = ps_nivell2_7d,
  assay_name = "category_abun_2_pruned",      
  fix_formula = "Treatment",  
  p_adj_method = "fdr",
  group = "Treatment",
  prv_cut = 0.1,              
  lib_cut = 0,                
  neg_lb = TRUE,
  struc_zero = TRUE,
  pairwise = TRUE
)

out_14d <- ancombc2(
  data = ps_nivell2_14d,
  assay_name = "category_abun_2_pruned",      
  fix_formula = "Treatment",  
  p_adj_method = "fdr",
  group = "Treatment",
  prv_cut = 0.1,              
  lib_cut = 0,                
  neg_lb = TRUE,
  struc_zero = TRUE,
  pairwise = TRUE
)

out_Base <- ancombc2(
  data = ps_nivell2_Base,
  assay_name = "category_abun_2_pruned",      
  fix_formula = "Timepoint",  
  p_adj_method = "fdr",
  group = "Timepoint",
  prv_cut = 0.1,              
  lib_cut = 0,                
  neg_lb = TRUE,
  struc_zero = TRUE,
  pairwise = TRUE
)

out_MQ <- ancombc2(
  data = ps_nivell2_MQ,
  assay_name = "category_abun_2_pruned",      
  fix_formula = "Timepoint",  
  p_adj_method = "fdr",
  group = "Timepoint",
  prv_cut = 0.1,              
  lib_cut = 0,                
  neg_lb = TRUE,
  struc_zero = TRUE,
  pairwise = TRUE
)

out_KNO3 <- ancombc2(
  data = ps_nivell2_KNO3,
  assay_name = "category_abun_2_pruned",      
  fix_formula = "Timepoint",  
  p_adj_method = "fdr",
  group = "Timepoint",
  prv_cut = 0.1,              
  lib_cut = 0,                
  neg_lb = TRUE,
  struc_zero = TRUE,
  pairwise = TRUE
)

out_CDQ <- ancombc2(
  data = ps_nivell2_CDQ,
  assay_name = "category_abun_2_pruned",      
  fix_formula = "Timepoint",  
  p_adj_method = "fdr",
  group = "Timepoint",
  prv_cut = 0.1,              
  lib_cut = 0,                
  neg_lb = TRUE,
  struc_zero = TRUE,
  pairwise = TRUE
)

# Export results to tables to work outside R. 
write.csv(out_7d$res_pair, "/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/resultats_7d.csv", sep = ";")
write.csv(out_14d$res_pair, "/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/resultats_14d.csv", sep = ";")
write.csv(out_KNO3$res, "/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/resultats_KNO3.csv", sep = ";")
write.csv(out_Base$res, "/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/resultats_Base.csv", sep = ";")
write.csv(out_MQ$res, "/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/resultats_MQ.csv", sep = ";")
write.csv(out_CDQ$res, "/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/resultats_CDQ.csv", sep = ";")

# Let's do the same at the top level.
keep <- rowMeans(category_abun_top > 0) > 0.1 
category_abun_top_pruned <- category_abun_top[keep, ]
meta <- read.csv("/home/alexarredondo/Documents/2025_Nitrats_v6/metadata.csv", sep = ";")

noms_files <- category_abun_top_pruned$top_level_category
category_abun_top_pruned <- category_abun_top_pruned[,-1]
rownames(category_abun_top_pruned) <- noms_files
rownames(meta) <- meta$Sample

otu_nivell_top <- otu_table(as.matrix(category_abun_top_pruned), taxa_are_rows = TRUE)
sam_nivell_top <- sample_data(meta)
ps_nivell_top <- phyloseq(otu_nivell_top, sam_nivell_top)
ps_nivell_top_7d <- subset_samples(ps_nivell_top, Timepoint == "7d")
ps_nivell_top_14d <- subset_samples(ps_nivell_top, Timepoint == "14d")
ps_nivell_top_MQ <- subset_samples(ps_nivell_top, Treatment == "MQ")
ps_nivell_top_CDQ <- subset_samples(ps_nivell_top, Treatment == "CDQ")
ps_nivell_top_KNO3 <- subset_samples(ps_nivell_top, Treatment == "KNO3")
ps_nivell_top_Base <- subset_samples(ps_nivell_top, Treatment == "Base")

out_7d_top <- ancombc2(
  data = ps_nivell_top_7d,
  assay_name = "category_abun_2_pruned",      
  fix_formula = "Treatment",  
  p_adj_method = "fdr",
  group = "Treatment",
  prv_cut = 0.1,              
  lib_cut = 0,                
  neg_lb = TRUE,
  struc_zero = TRUE,
  pairwise = TRUE
)

out_14d_top <- ancombc2(
  data = ps_nivell_top_14d,
  assay_name = "category_abun_2_pruned",      
  fix_formula = "Treatment",  
  p_adj_method = "fdr",
  group = "Treatment",
  prv_cut = 0.1,              
  lib_cut = 0,                
  neg_lb = TRUE,
  struc_zero = TRUE,
  pairwise = TRUE
)

out_Base_top <- ancombc2(
  data = ps_nivell_top_Base,
  assay_name = "category_abun_2_pruned",      
  fix_formula = "Timepoint",  
  p_adj_method = "fdr",
  group = "Timepoint",
  prv_cut = 0.1,              
  lib_cut = 0,                
  neg_lb = TRUE,
  struc_zero = TRUE,
  pairwise = TRUE
)

out_MQ_top <- ancombc2(
  data = ps_nivell_top_MQ,
  assay_name = "category_abun_2_pruned",      
  fix_formula = "Timepoint",  
  p_adj_method = "fdr",
  group = "Timepoint",
  prv_cut = 0.1,              
  lib_cut = 0,                
  neg_lb = TRUE,
  struc_zero = TRUE,
  pairwise = TRUE
)

out_KNO3_top <- ancombc2(
  data = ps_nivell_top_KNO3,
  assay_name = "category_abun_2_pruned",      
  fix_formula = "Timepoint",  
  p_adj_method = "fdr",
  group = "Timepoint",
  prv_cut = 0.1,              
  lib_cut = 0,                
  neg_lb = TRUE,
  struc_zero = TRUE,
  pairwise = TRUE
)

out_CDQ_top <- ancombc2(
  data = ps_nivell_top_CDQ,
  assay_name = "category_abun_2_pruned",      
  fix_formula = "Timepoint",  
  p_adj_method = "fdr",
  group = "Timepoint",
  prv_cut = 0.1,              
  lib_cut = 0,                
  neg_lb = TRUE,
  struc_zero = TRUE,
  pairwise = TRUE
)

write.csv(out_7d_top$res_pair, "/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/resultats_7d_top.csv", sep = ";")
write.csv(out_14d_top$res_pair, "/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/resultats_14d_top.csv", sep = ";")
write.csv(out_KNO3_top$res, "/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/resultats_KNO3_top.csv", sep = ";")
write.csv(out_Base_top$res, "/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/resultats_Base_top.csv", sep = ";")
write.csv(out_MQ_top$res, "/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/resultats_MQ_top.csv", sep = ";")
write.csv(out_CDQ_top$res, "/home/alexarredondo/Documents/2025_Nitrats_v6/picrust2/resultats_CDQ_top.csv", sep = ";")


# Let's visualize the data using the same type of heatmaps that we used when visualizing differential relative abundance.
taula_grafic_picrust <- read.csv("/Users/alex.arredondo/OneDrive - DENTAID/Projectes/Interns/Ferrum/Paper/taula_picrust2_figura.csv", sep = ";")

lo_picrust = floor(min(taula_grafic_picrust$LFC))
up_picrust = ceiling(max(taula_grafic_picrust$LFC))
mid = 0
taula_grafic_picrust %>%
  ggplot(aes(x = Comparison, y = Pathway, fill = LFC)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo_picrust, up_picrust),
                       name = NULL) +
  geom_text(aes(Comparison, Pathway, label = round(LFC,2)), size = 2) +
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = "Log fold change dels diferents tractaments \nseparats per timepoints - p-value < 0.05") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_discrete(limits=rev)

####################################### Core microbiome #######################################

# We are going to assess the core microbiota of each group using indicspecies. We start by filtering taxa with a prevalence lower than 10%.
ps_glom_genus_v1_pruned_filtered <- filter_taxa(ps_glom_genus_v1_pruned, function(x) sum(x > 0) > (0.1 * length(x)), prune = TRUE)

# Normalize to relative abundance. 
seq_rel <- transform_sample_counts(ps_glom_genus_v1_pruned_filtered, function(x) x / sum(x))

# Extract the asv table and the grouping by which indicspecies will manage the data.
otu_rel <- as.data.frame(otu_table(seq_rel))
if(taxa_are_rows(seq_rel)) otu_rel <- t(otu_rel)
groups <- as.factor(sample_data(ps_glom_genus_v1_pruned_filtered)$Timepoint_Treatment
                    
# Normalize to hellinger (reported as the best normalization for indicspecies).
otu_hel <- decostand(otu_rel, method = "hellinger")
                    
# Run indicspecies.
set.seed(123)
res_indval <- multipatt(otu_hel, groups, func = "r.g", control = how(nperm = 999))
summary(res_indval)

####################################### Co-occurence networks #######################################
# The best way to analyze any type of networks in microbial comunities is through spiec-easi, which take into account the compositional nature of the data. However, it needs a large dataset which was not the case in this study. Therefore, we have opted for a correlation analysis with clr-normalized data and the Spearman test.

# Prepare a function to filter by prevalence and counts (or not, but be ready to have unreadable plots)
filtra_taxa <- function(ps_object, min_prevalence = 0.2, min_count = 10) {
  otu <- as(otu_table(ps_object), "matrix")
  if(taxa_are_rows(ps_object)) otu <- t(otu)
  
  preval <- apply(otu, 2, function(x) mean(x > 0))
  taxa_keep <- colnames(otu)[preval >= min_prevalence]
  ps_object <- prune_taxa(taxa_keep, ps_object)
  
  taxa_total <- taxa_sums(ps_object)
  ps_object <- prune_taxa(taxa_total >= min_count, ps_object)
  return(ps_object)
}

# Prepare a function to normalize to CLR
clr_transform <- function(ps_object, pseudo_count = 1) {
  otu <- as(otu_table(ps_object), "matrix")
  if(taxa_are_rows(ps_object)) otu <- t(otu) # make sure the table is well oriented.
  otu <- otu + pseudo_count # add pseudocounts.
  clr_mat <- t(apply(otu, 1, compositions::clr)) # transform to clr.
  clr_mat <- as.matrix(clr_mat)
  mode(clr_mat) <- "numeric" # make sure all data is numeric.
  clr_mat[!is.finite(clr_mat)] <- 0 # add 0's to NA or -inf if something happened.
  clr_mat <- clr_mat[, apply(clr_mat, 2, sd, na.rm = TRUE) > 0, drop = FALSE] # delete columns (taxa) without variation, which are most likely 0's.
  return(clr_mat)
}

# Prepare a function for the Spearman test.
cor_spearman <- function(clr_mat, p_adj_method = "BH", p_thresh = 0.05) {
  n_taxa <- ncol(clr_mat)
  cor_mat <- matrix(NA, n_taxa, n_taxa)
  p_mat   <- matrix(NA, n_taxa, n_taxa)
  colnames(cor_mat) <- colnames(clr_mat)
  rownames(cor_mat) <- colnames(clr_mat)
  colnames(p_mat) <- colnames(clr_mat)
  rownames(p_mat) <- colnames(clr_mat)
  
  for (i in 1:(n_taxa - 1)) {
    for (j in (i + 1):n_taxa) {
      test <- try(cor.test(as.numeric(clr_mat[, i]),
                           as.numeric(clr_mat[, j]),
                           method = "spearman"),
                  silent = TRUE)
      if (!inherits(test, "try-error")) {
        cor_mat[i, j] <- cor_mat[j, i] <- test$estimate
        p_mat[i, j]   <- p_mat[j, i]   <- test$p.value
      }
    }
  }
  
  p_mat_adj <- matrix(p.adjust(p_mat, method = p_adj_method), n_taxa, n_taxa)
  colnames(p_mat_adj) <- colnames(clr_mat)
  rownames(p_mat_adj) <- colnames(clr_mat)
  
  cor_mat_filtered <- cor_mat
  cor_mat_filtered[p_mat_adj > p_thresh] <- 0
  
  return(list(cor = cor_mat_filtered, p_adj = p_mat_adj))
}

# Now merge all functions into one.
pipeline_spearman <- function(ps_object,
                              min_prevalence = 0.2,
                              min_count = 10,
                              pseudo_count = 1,
                              p_adj_method = "BH",
                              p_thresh = 0.05,
                              abs_threshold = 0.3) {
  ps_filt <- filtra_taxa(ps_object, min_prevalence, min_count)
  clr_mat <- clr_transform(ps_filt, pseudo_count)
  cor_res <- cor_spearman(clr_mat, p_adj_method, p_thresh)
  ig <- cor_to_igraph(cor_res$cor, abs_threshold)
  
  return(list(cor_matrix = cor_res$cor,
              p_adj_matrix = cor_res$p_adj,
              ps_filt = ps_filt))
}


# Finally, a function to extract correlations with p values < 0.05
obte_correlacions_significatives <- function(res,
                                             p_thresh = 0.05,
                                             abs_r_thresh = 0.3) {
  cor_df <- as.data.frame(as.table(res$cor_matrix))
  p_df   <- as.data.frame(as.table(res$p_adj_matrix))
  df <- merge(cor_df, p_df, by = c("Var1", "Var2"))
  colnames(df) <- c("Taxa1", "Taxa2", "Spearman_rho", "p_adj")
  # Delete symetric data.
  df <- df[df$Taxa1 != df$Taxa2, ]
  df <- df[!duplicated(t(apply(df[,1:2], 1, sort))), ]
  df_sig <- df %>%
    dplyr::filter(!is.na(Spearman_rho),
                  !is.na(p_adj),
                  abs(Spearman_rho) >= abs_r_thresh,
                  p_adj < p_thresh) %>%
    dplyr::arrange(p_adj)
  return(df_sig)
}

# Subset the data that will be analyzed.
spearman_basal <- subset_samples(ps_glom_genus_v1_pruned, Timepoint_Treatment == "0d_Basal")
spearman_7d_MQ <- subset_samples(ps_glom_genus_v1_pruned, Timepoint_Treatment == "7d_MQ")
spearman_14d_MQ <- subset_samples(ps_glom_genus_v1_pruned, Timepoint_Treatment == "14d_MQ")
spearman_7d_KNO3 <- subset_samples(ps_glom_genus_v1_pruned, Timepoint_Treatment == "7d_KNO3")
spearman_14d_KNO3 <- subset_samples(ps_glom_genus_v1_pruned, Timepoint_Treatment == "14d_KNO3")
spearman_7d_Base <- subset_samples(ps_glom_genus_v1_pruned, Timepoint_Treatment == "7d_Base")
spearman_14d_Base <- subset_samples(ps_glom_genus_v1_pruned, Timepoint_Treatment == "14d_Base")
spearman_7d_CDQ <- subset_samples(ps_glom_genus_v1_pruned, Timepoint_Treatment == "7d_CDQ")
spearman_14d_CDQ <- subset_samples(ps_glom_genus_v1_pruned, Timepoint_Treatment == "14d_CDQ")

# Run the pipeline and export the table to work with it outside R.
res_basal <- pipeline_spearman(spearman_basal)
df_res_basal <- obte_correlacions_significatives(res_basal)
write.table(df_res_basal, "/home/alexarredondo/Documents/2025_Nitrats_v6/resultats_spearman/res_basal.csv", sep = ";")

res_7d_MQ <- pipeline_spearman(spearman_7d_MQ)
df_res_7d_MQ <- obte_correlacions_significatives(res_7d_MQ)
write.table(df_res_7d_MQ, "/home/alexarredondo/Documents/2025_Nitrats_v6/resultats_spearman/res_7d_MQ.csv", sep = ";")

res_14d_MQ <- pipeline_spearman(spearman_14d_MQ)
df_res_14d_MQ <- obte_correlacions_significatives(res_14d_MQ)
write.table(df_res_14d_MQ, "/home/alexarredondo/Documents/2025_Nitrats_v6/resultats_spearman/res_14d_MQ.csv", sep = ";")

res_7d_KNO3 <- pipeline_spearman(spearman_7d_KNO3)
df_res_7d_KNO3 <- obte_correlacions_significatives(res_7d_KNO3)
write.table(df_res_7d_KNO3, "/home/alexarredondo/Documents/2025_Nitrats_v6/resultats_spearman/res_7d_KNO3.csv", sep = ";")

res_14d_KNO3 <- pipeline_spearman(spearman_14d_KNO3)
df_res_14d_KNO3 <- obte_correlacions_significatives(res_14d_KNO3)
write.table(df_res_14d_KNO3, "/home/alexarredondo/Documents/2025_Nitrats_v6/resultats_spearman/res_14d_KNO3.csv", sep = ";")

res_7d_Base <- pipeline_spearman(spearman_7d_Base)
df_res_7d_Base <- obte_correlacions_significatives(res_7d_Base)
write.table(df_res_7d_Base, "/home/alexarredondo/Documents/2025_Nitrats_v6/resultats_spearman/res_7d_Base.csv", sep = ";")

res_14d_Base <- pipeline_spearman(spearman_14d_Base)
df_res_14d_Base <- obte_correlacions_significatives(res_14d_Base)
write.table(df_res_14d_Base, "/home/alexarredondo/Documents/2025_Nitrats_v6/resultats_spearman/res_14d_Base.csv", sep = ";")

res_7d_CDQ <- pipeline_spearman(spearman_7d_CDQ)
df_res_7d_CDQ <- obte_correlacions_significatives(res_7d_CDQ)
write.table(df_res_7d_CDQ, "/home/alexarredondo/Documents/2025_Nitrats_v6/resultats_spearman/res_7d_CDQ.csv", sep = ";")

res_14d_CDQ <- pipeline_spearman(spearman_14d_CDQ)
df_res_14d_CDQ <- obte_correlacions_significatives(res_14d_CDQ)
write.table(df_res_14d_CDQ, "/home/alexarredondo/Documents/2025_Nitrats_v6/resultats_spearman/res_14d_CDQ.csv", sep = ";")


# Now let's plot it.
taula_corplots <- read.csv("/Users/alex.arredondo/OneDrive - DENTAID/Projectes/Interns/Ferrum/Paper/taula_corplot_figura.csv", sep = ";")

# Add a new variable to set the colour for positive/negative correlations. 
taula_corplots$color <- ifelse(taula_corplots$Spearman.Rho > 0, "forestgreen", "firebrick")

taula_corplots_baseline <- taula_corplots[taula_corplots$Comparison == "Baseline", ]
taula_corplots_mq7d <- taula_corplots[taula_corplots$Comparison == "MQ 7d", ]
taula_corplots_mq14d <- taula_corplots[taula_corplots$Comparison == "MQ 14d", ]
taula_corplots_base7d <- taula_corplots[taula_corplots$Comparison == "Base 7d", ]
taula_corplots_kno314d <- taula_corplots[taula_corplots$Comparison == "KNO3 14d", ]
taula_corplots_cdq7d <- taula_corplots[taula_corplots$Comparison == "CDQ 7d", ]
taula_corplots_cdq14d <- taula_corplots[taula_corplots$Comparison == "CDQ 14d", ]


# Circos is a bit tricky, so I suggest running each plot in a separate chunk if you want to reproduce the same output. 

# Chunk1
circos.clear()

circos.par(canvas.xlim = c(-1.6, 1.6),
           canvas.ylim = c(-1., 1.6))

chordDiagram(
  x = taula_corplots_baseline[, c("Taxa.1", "Taxa.2", "Spearman.Rho")],
  grid.col = NULL,
  col = taula_corplots_baseline$color,
  transparency = 0.3,
  directional = 0,
  annotationTrack = "grid",
  preAllocateTracks = 1
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    circos.text(
      x = mean(get.cell.meta.data("xlim")),
      y = 0,
      labels = sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.6
    )
  },
  bg.border = NA
)

title("Correlation plot - Baseline")

# Chunk2
circos.clear()

circos.par(canvas.xlim = c(-1.6, 1.6),
           canvas.ylim = c(-1., 1.6))

chordDiagram(
  x = taula_corplots_mq7d[, c("Taxa.1", "Taxa.2", "Spearman.Rho")],
  grid.col = NULL,
  col = taula_corplots_mq7d$color,
  transparency = 0.3,
  directional = 0,
  annotationTrack = "grid",
  preAllocateTracks = 1
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    circos.text(
      x = mean(get.cell.meta.data("xlim")),
      y = 0,
      labels = sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.6
    )
  },
  bg.border = NA
)

title("Correlation plot - 7d - MQ")

# Chunk3
circos.clear()

circos.par(canvas.xlim = c(-1.6, 1.6),
           canvas.ylim = c(-1., 1.6))

chordDiagram(
  x = taula_corplots_mq14d[, c("Taxa.1", "Taxa.2", "Spearman.Rho")],
  grid.col = NULL,
  col = taula_corplots_mq14d$color,
  transparency = 0.3,
  directional = 0,
  annotationTrack = "grid",
  preAllocateTracks = 1
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    circos.text(
      x = mean(get.cell.meta.data("xlim")),
      y = 0,
      labels = sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.6
    )
  },
  bg.border = NA
)

title("Correlation plot - 14d - MQ")

# Chunk4
circos.clear()

circos.par(canvas.xlim = c(-1.6, 1.6),
           canvas.ylim = c(-1., 1.6))

chordDiagram(
  x = taula_corplots_base7d[, c("Taxa.1", "Taxa.2", "Spearman.Rho")],
  grid.col = NULL,
  col = taula_corplots_base7d$color,
  transparency = 0.3,
  directional = 0,
  annotationTrack = "grid",
  preAllocateTracks = 1
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    circos.text(
      x = mean(get.cell.meta.data("xlim")),
      y = 0,
      labels = sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.6
    )
  },
  bg.border = NA
)

title("Correlation plot - 7d - Base")

# Chunk5
circos.clear()

circos.par(canvas.xlim = c(-1.6, 1.6),
           canvas.ylim = c(-1., 1.6))

chordDiagram(
  x = taula_corplots_kno314d[, c("Taxa.1", "Taxa.2", "Spearman.Rho")],
  grid.col = NULL,
  col = taula_corplots_kno314d$color,
  transparency = 0.3,
  directional = 0,
  annotationTrack = "grid",
  preAllocateTracks = 1
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    circos.text(
      x = mean(get.cell.meta.data("xlim")),
      y = 0,
      labels = sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.6
    )
  },
  bg.border = NA
)

title("Correlation plot - 14d - KNO3")

# Chunk6
circos.clear()

circos.par(canvas.xlim = c(-1.6, 1.6),
           canvas.ylim = c(-1., 1.6))

chordDiagram(
  x = taula_corplots_cdq7d[, c("Taxa.1", "Taxa.2", "Spearman.Rho")],
  grid.col = NULL,
  col = taula_corplots_cdq7d$color,
  transparency = 0.3,
  directional = 0,
  annotationTrack = "grid",
  preAllocateTracks = 1
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    circos.text(
      x = mean(get.cell.meta.data("xlim")),
      y = 0,
      labels = sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.6
    )
  },
  bg.border = NA
)

title("Correlation plot - 7d - CDQ")

# Chunk7
circos.clear()

circos.par(canvas.xlim = c(-1.6, 1.6),
           canvas.ylim = c(-1., 1.6))

chordDiagram(
  x = taula_corplots_cdq14d[, c("Taxa.1", "Taxa.2", "Spearman.Rho")],
  grid.col = NULL,
  col = taula_corplots_cdq14d$color,
  transparency = 0.3,
  directional = 0,
  annotationTrack = "grid",
  preAllocateTracks = 1
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    circos.text(
      x = mean(get.cell.meta.data("xlim")),
      y = 0,
      labels = sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.6
    )
  },
  bg.border = NA
)

title("Correlation plot - 14d - CDQ")

# It is always helpful to have this function at the end of your .rmd file to save and/or load your ambient faster.
save.image("/home/alexarredondo/Documents/2025_Nitrats_v6/ambient_nitrats_v6_1.Rdata")
Septembre2025_
Genomica2021
load("/home/alexarredondo/Documents/2025_Nitrats_v6/ambient_nitrats_v6_1.Rdata")
session_info() # Very important to report the package's version you are using.

