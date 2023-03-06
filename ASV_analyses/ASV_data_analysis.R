# PHYLOSEQ ANALYSIS

# Author: Vanessa Brisson, PhD
# Email: brisson2@llnl.gov
# Date Created: 2023-01-11

# Copyright (c) Vanessa Brisson, 2023

# ------------------------------------------------------------------------------
# LOAD PACKAGES

library(dplyr)
library(phyloseq)
library(ggplot2)
library(vegan)
library(metagenomeSeq)
library(ANCOMBC)
library(microbiome)
library(tibble)
library(ape)
library(DECIPHER)
library(phangorn)
library(ggtree)

#-------------------------------------------------------------------------------
# LOAD DATA

# change this path variable to wherever the ASV data are stored
path <- 'C:Users/username/Documents'

# The ASVtable in and taxonomy tables are the output from 
# "raw_sequence_data_analysis_and_ASV_table_generation.R" and are provided in 
# the repository for convenience.

seqtab <- readRDS(file.path(path, 'ASVtable.rds'))    # sequence table from 
                                                      #  DADA2 analysis

tax <- readRDS(file.path(path, 'taxtable.rds'))       # taxonomy table from 
                                                      #  DADA2 analysis - 
                                                      #  Silva database

tax2 <- readRDS(file.path(path, 'taxtable_RDP.rds')) # alternate taxonomy table 
                                                      #  from DADA2 analaysis - 
                                                      #  RDP database

tax3 <- cbind(tax, tax2[,"Genus"])
colnames(tax3)[7] <- "RDPGenus"

# metadata file
metadata <- read.csv(file.path(path, 'metadata.csv'), header=TRUE, row.names=1)
metadata$time2 <- as.factor(metadata$time)

# create phyloseq object from loaded data
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(tax3))

# add dna seqs to phyloseq object and rename ASVs
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
tax4 <- cbind(tax_table(ps), taxa_names(ps))
colnames(tax4)[8] <- "Species"
tax_table(ps) <- tax_table(tax4)


# remove Microchloropsis data, this is other data not included in this study
ps <- ps %>%
  subset_samples(algae != 'Microchloropsis')

# remove inoculum sequences for this analysis
ps <- ps %>%
  subset_samples(algae != 'inoculum')

# add sample data columns indicating whether alga and/or exudates are present
sample_data(ps)$algae_present <- sample_data(ps)$algae == "Phaeodactylum"
sample_data(ps)$exudates_present <- sample_data(ps)$algae != "none"

#-------------------------------------------------------------------------------
# FILTER DATA

# filter to include only Bacterial and Acchaeal sequences and remove any 
#  Chloroplast of Mitochondrial reads
ps <- ps %>% 
  subset_taxa((Kingdom == "Bacteria" | Kingdom == "Archaea") &
                (is.na(Order) |  Order != "Chloroplast" ) &
                (is.na(Family) | Family != "Mitochondria"))

total_reads = sample_sums(ps)
hist(total_reads, breaks=100)
print(total_reads[sample_sums(ps) < 5e3]) # check to see if any samples have 
                                          # fewer than 5000 reads
rarecurve((otu_table(ps)),   # rarefaction curve
          step=50, cex=0.5)

# keep taxa detected in at least 4 samples, this can be changed and was 
# determined based on the experimental design and the data
ps <- filter_taxa(ps, function(x) sum(x >= 10) >= 4, TRUE) 
ps

# Generate a phylogenetic tree of the ASV sequences
# This tree construction is based on tree on Callahan et. al. 2016
alignment <- AlignSeqs(refseq(ps), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)                           # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
phy_tree(ps) <- phy_tree(fitGTR$tree)

#-------------------------------------------------------------------------------
# CORRECT FOR COMPOSITIONALLITY USING ANCOMBC

out <- suppressWarnings(ancombc(phyloseq = ps, tax_level = 'Species',
                                formula = "algae + metabolite + time", 
                                p_adj_method = "fdr",alpha = 0.05,
                                prv_cut = 0.10, lib_cut = 1000,
                                group = "algae", struc_zero = FALSE))
samp_frac = out$samp_frac
samp_frac[is.na(samp_frac)] = 0 
log_obs_abn = log(abundances(ps) + 1) 
log_obs_abn_adj = t(t(log_obs_abn) - samp_frac)

# This is the compositionality corrected phyloseq object
ps.ABC <- phyloseq(otu_table(t(exp(log_obs_abn)), taxa_are_rows=FALSE), 
                           sample_data(sample_data(ps)), 
                           tax_table(tax_table(ps)))

# this is the compositionality corrected phyloseq object converted to relative
# abundance
ps.ABC.RA <- transform_sample_counts(ps.ABC, function(x) x / sum(x))

# Save the compositionality corrected ASV table for analysis in other scripts
write.csv(as.data.frame(otu_table(ps.ABC.RA)), file.path(path, 'OTUtable_ABC_RA.csv'))
write.csv(as.data.frame(sample_data(ps.ABC.RA)), file.path(path, 'SampleData.csv'))
write.csv(as.data.frame(tax_table(ps.ABC.RA)), file.path(path, 'TaxTable.csv'))

#-------------------------------------------------------------------------------
# SUBSET DATA

# Subset data to look only at the no metabolite added experiment
ps.ctrl <- ps %>% subset_samples(metabolite == 'none')
ps.ABC.ctrl <- ps.ABC %>% subset_samples(metabolite == 'none') 

#plot the phylogenetic tree of ASVs
plot_tree(merge_samples(ps.ctrl, 'metabolite'), ladderize="left", 
          color="Class", label.tips="taxa_names", nodelabf=nodeplotblank,
          title='Maximum Likelihood Tree of ASV Sequences') +
  scale_size_continuous(range = c(0.1, 0.1))

#-------------------------------------------------------------------------------
# ALPHA DIVERSITY

# ANOVA
sample_data(ps.ABC.ctrl)$Shannon <- estimate_richness(ps.ABC.ctrl,measures='Shannon')$Shannon
sample_data_frame_ps <- data.frame(sample_data(ps.ABC.ctrl))
options(contrasts = c('contr.sum','contr.poly'))
model <- lm(Shannon ~ algae*time2, sample_data_frame_ps)
drop1(model, .~., test='F')

# Tukey's test
tukey <- TukeyHSD(aov(Shannon ~ algae*time2, data=sample_data_frame_ps))

# Save data to make a prettier plot
write.csv(sample_data(ps.ABC.ctrl)[,'Shannon'], file.path(path, 'ctrShannon_ABC.csv'))

#-------------------------------------------------------------------------------
# BETA DIVERSITY

# NO METABOLITE ADDED EXPERIMENT


# PERMANOVA
set.seed(1)
ps.ABC.ctrl.dist <- phyloseq::distance(ps.ABC.ctrl, method='jsd')
sampledf <- data.frame(sample_data(ps.ABC.ctrl))
adonis2(ps.ABC.ctrl.dist ~ algae*time2, data=sampledf, permutations = 9999)

# PCoA with negative eigenvalue correction
ps.ABC.ctrl.PCoA <- pcoa(ps.ABC.ctrl.dist, correction='cailliez')
plot(ps.ABC.ctrl.PCoA$vectors.cor[,'Axis.1'], 
     ps.ABC.ctrl.PCoA$vectors.cor[,'Axis.2'], type="p")

write.csv(as.data.frame(ps.ABC.ctrl.PCoA$vectors.cor) %>% 
            select(c('Axis.1','Axis.2')), 
          file.path(path, 'ctrlPCoA_ABC_cor.csv'))
print('Percent Variance Explained')
print(ps.ABC.ctrl.PCoA$values$Corr_eig[1]/sum(ps.ABC.ctrl.PCoA$values$Corr_eig))
print(ps.ABC.ctrl.PCoA$values$Corr_eig[2]/sum(ps.ABC.ctrl.PCoA$values$Corr_eig))

# METABOLITE ADDITION EXPERIMENTs

# PERMANOVA
set.seed(1)
ps.ABC.dist <- phyloseq::distance(ps.ABC, method='jsd')
sampledf <- data.frame(sample_data(ps.ABC))
adonis2(ps.ABC.dist ~ algae*metabolite*time2, data=sampledf, permutations = 9999)

# PCoA with negative eigenvalue correction
ps.ABC.PCoA <- pcoa(ps.ABC.dist, correction='cailliez')
plot(ps.ABC.PCoA$vectors.cor[,'Axis.1'], 
     ps.ABC.PCoA$vectors.cor[,'Axis.2'], type="p")

write.csv(as.data.frame(ps.ABC.PCoA$vectors.cor) %>% 
            select(c('Axis.1','Axis.2', 'Axis.3')), 
          file.path(path, 'allPCoA_ABC_cor.csv'))

print('Percent Variance Explained')
print(ps.ABC.PCoA$values$Corr_eig[1]/sum(ps.ABC.PCoA$values$Corr_eig))
print(ps.ABC.PCoA$values$Corr_eig[2]/sum(ps.ABC.PCoA$values$Corr_eig))
print(ps.ABC.PCoA$values$Corr_eig[3]/sum(ps.ABC.PCoA$values$Corr_eig))

#-------------------------------------------------------------------------------
# DIFFERENTIAL ABUNDANCE

# IMPORTANT: Note that the input phyloseq object that is not corrected for 
# compositionality. This is because the ANCOMBC differential abundance analysis 
# includes that correction, and we don't want to repeat the correction on the 
# already corrected data.

# NO METABOLITE ADDED EXPERIMENT

# Conduct pairwise differential abundance comparisons between the three conditions: 
#     P = P. tricornutum present
#     PS = P. tricornutum spent medium
#     C = no alga or exudates
out_PStoP <- suppressWarnings(ancombc(phyloseq = ps.ctrl %>% 
                                       subset_samples(algae != "none"), 
                                      tax_level = 'Species',
                                      formula = "time + algae", 
                                      p_adj_method = "fdr", alpha = 0.05,
                                      prv_cut = 0.10, lib_cut = 1000, 
                                      group = 'algae', struc_zero = TRUE))
out_PStoC <- suppressWarnings(ancombc(phyloseq = ps.ctrl%>% 
                                        subset_samples(algae != "Phaeodactylum"), 
                                      tax_level = 'Species',
                                      formula = "time + algae", 
                                      p_adj_method = "fdr",alpha = 0.05,
                                      prv_cut = 0.10, lib_cut = 1000, 
                                      group = 'algae', struc_zero = TRUE))
out_PtoC <- suppressWarnings(ancombc(phyloseq = ps.ctrl %>% 
                                       subset_samples(algae != "Phaeodactylum spent medium"), 
                                     tax_level = 'Species',
                                     formula = "time + algae", 
                                     p_adj_method = "fdr",alpha = 0.05,
                                     prv_cut = 0.10, lib_cut = 1000, 
                                     group = 'algae', struc_zero = TRUE))

# put the reults together for output
DA.PStoP <- merge(out_PStoP$res$diff_abn['algae1'], 
                  out_PStoP$res$lfc['algae1'], by=0) %>%
  rename('algae1.x' = 'PStoP') %>%
  rename('algae1.y' = 'LFC.PStoP')
rownames(DA.PStoP) <- DA.PStoP$Row.names

DA.PtoC <- merge(out_PtoC$res$diff_abn['algae1'], 
                  out_PtoC$res$lfc['algae1'], by=0) %>%
  rename('algae1.x' = 'PtoC') %>%
  rename('algae1.y' = 'LFC.PtoC')
rownames(DA.PtoC) <- DA.PtoC$Row.names

DA.PStoC <- merge(out_PStoC$res$diff_abn['algae1'], 
                 out_PStoC$res$lfc['algae1'], by=0) %>%
  rename('algae1.x' = 'PStoC') %>%
  rename('algae1.y' = 'LFC.PStoC')
rownames(DA.PStoC) <- DA.PStoC$Row.names

means <- colMeans(otu_table(ps.ABC.RA))   # mean relative abundance (after
                                          # adjustment for compositionallity)

# merge differential abundance analysis results
DA.ctrl <- merge(DA.PStoP, DA.PtoC, by=0, all=TRUE) %>%
  subset(select = c('Row.names','PStoP','LFC.PStoP',
                    'PtoC','LFC.PtoC')) %>%
  merge(DA.PStoC, by.x=1, by.y=0, all=TRUE) %>%
  subset(select = c('Row.names', 'PStoP','LFC.PStoP',
                    'PtoC','LFC.PtoC', 
                    'PStoC','LFC.PStoC')) %>%
  merge(as.data.frame(means), by.x=1,by.y=0, all=TRUE) %>%
  merge(as.data.frame(tax_table(ps.ABC.RA)), all=TRUE, by.x=1, by.y=0)
rownames(DA.ctrl) <- DA.ctrl$Row.names

# save the results for plotting
write.csv(DA.ctrl, file.path(path, 'DA_Ctrl.csv'))

# METABOLITE ADDITION EXPERIMENTS

# 4-HYDROXYBENZOIC ACID

# Conduct pairwise differential abundance comparisons
out.4hba.all <- suppressWarnings(ancombc(phyloseq = ps %>% 
                                           subset_samples(metabolite != 'lumichrome'), 
                                     tax_level = 'Species',
                                     formula = "algae + time + metabolite", 
                                     p_adj_method = "fdr",alpha = 0.05,
                                     prv_cut = 0.10, lib_cut = 1000, 
                                     group = 'metabolite', struc_zero = TRUE))

out.4hba.C <- suppressWarnings(ancombc(phyloseq = ps %>% 
                                           subset_samples(metabolite != 'lumichrome') %>%
                                           subset_samples(algae == 'none'), 
                                         tax_level = 'Species',
                                         formula = "time + metabolite", 
                                         p_adj_method = "fdr",alpha = 0.05,
                                         prv_cut = 0.10, lib_cut = 1000, 
                                         group = 'metabolite', struc_zero = TRUE))

out.4hba.PS <- suppressWarnings(ancombc(phyloseq = ps %>% 
                                           subset_samples(metabolite != 'lumichrome') %>%
                                           subset_samples(algae == 'Phaeodactylum spent medium'),
                                         tax_level = 'Species',
                                         formula = "time + metabolite", 
                                         p_adj_method = "fdr",alpha = 0.05,
                                         prv_cut = 0.10, lib_cut = 1000, 
                                         group = 'metabolite', struc_zero = TRUE))

out.4hba.P <- suppressWarnings(ancombc(phyloseq = ps %>% 
                                           subset_samples(metabolite != 'lumichrome') %>%
                                           subset_samples(algae == 'Phaeodactylum'), 
                                         tax_level = 'Species',
                                         formula = "time + metabolite", 
                                         p_adj_method = "fdr",alpha = 0.05,
                                         prv_cut = 0.10, lib_cut = 1000, 
                                         group = 'metabolite', struc_zero = TRUE))

# merge differential abundance analysis results
DA.4hba.all <- merge(out.4hba.all$res$diff_abn['metabolite1'], 
                  out.4hba.all$res$q_val['metabolite1'], by=0) %>%
  rename('metabolite1.x' = 'all') %>%
  rename('metabolite1.y' = 'q_val.all') %>%
  column_to_rownames(var='Row.names') %>%
  merge(out.4hba.all$res$lfc['metabolite1'], by=0) %>%
  rename('metabolite1' = 'LFC.all') %>%
  column_to_rownames(var='Row.names')

DA.4hba.C <- merge(out.4hba.C$res$diff_abn['metabolite1'], 
                     out.4hba.C$res$q_val['metabolite1'], by=0) %>%
  rename('metabolite1.x' = 'C') %>%
  rename('metabolite1.y' = 'q_val.C') %>%
  column_to_rownames(var='Row.names') %>%
  merge(out.4hba.C$res$lfc['metabolite1'], by=0) %>%
  rename('metabolite1' = 'LFC.C') %>%
  column_to_rownames(var='Row.names')

DA.4hba.PS <- merge(out.4hba.PS$res$diff_abn['metabolite1'], 
                     out.4hba.PS$res$q_val['metabolite1'], by=0) %>%
  rename('metabolite1.x' = 'PS') %>%
  rename('metabolite1.y' = 'q_val.PS') %>%
  column_to_rownames(var='Row.names') %>%
  merge(out.4hba.PS$res$lfc['metabolite1'], by=0) %>%
  rename('metabolite1' = 'LFC.PS') %>%
  column_to_rownames(var='Row.names')

DA.4hba.P <- merge(out.4hba.P$res$diff_abn['metabolite1'], 
                     out.4hba.P$res$q_val['metabolite1'], by=0) %>%
  rename('metabolite1.x' = 'P') %>%
  rename('metabolite1.y' = 'q_val.P') %>%
  column_to_rownames(var='Row.names') %>%
  merge(out.4hba.P$res$lfc['metabolite1'], by=0) %>%
  rename('metabolite1' = 'LFC.P') %>%
  column_to_rownames(var='Row.names')

means <- colMeans(otu_table(ps.ABC.RA))  # mean relative abundance (after
                                         # adjustment for compositionallity)

DA.4hba <- merge(DA.4hba.all, DA.4hba.C, by=0, all=TRUE) %>%
  subset(select = c('Row.names','all','q_val.all','LFC.all',
                    'C','q_val.C','LFC.C')) %>%
  merge(DA.4hba.PS, by.x=1, by.y=0, all=TRUE) %>%
  subset(select = c('Row.names', 'all','q_val.all','LFC.all',
                    'C','q_val.C','LFC.C', 
                    'PS','q_val.PS','LFC.PS')) %>%
  merge(DA.4hba.P, by.x=1, by.y=0, all=TRUE) %>%
  subset(select = c('Row.names', 'all','q_val.all','LFC.all',
                    'C','q_val.C','LFC.C', 
                    'PS','q_val.PS','LFC.PS',
                    'P','q_val.P','LFC.P')) %>%
  merge(as.data.frame(means), by.x=1,by.y=0, all=TRUE) %>%
  merge(as.data.frame(tax_table(ps.ABC.RA)), all=TRUE, by.x=1, by.y=0) %>%
  column_to_rownames(var='Row.names')

# save the results for plotting
write.csv(DA.4hba, file.path(path, 'DA_4hba.csv'))

# LUMICHROME

# Conduct pairwise differential abundance comparisons
out.lum.all <- suppressWarnings(ancombc(phyloseq = ps %>% 
                                           subset_samples(metabolite != '4hydroxybenzoic_acid'), 
                                         tax_level = 'Species',
                                         formula = "algae + time + metabolite", 
                                         p_adj_method = "fdr",alpha = 0.05,
                                         prv_cut = 0.10, lib_cut = 1000, 
                                         group = 'metabolite', struc_zero = TRUE))

out.lum.C <- suppressWarnings(ancombc(phyloseq = ps %>% 
                                         subset_samples(metabolite != '4hydroxybenzoic_acid') %>%
                                         subset_samples(algae == 'none'), 
                                       tax_level = 'Species',
                                       formula = "time + metabolite", 
                                       p_adj_method = "fdr",alpha = 0.05,
                                       prv_cut = 0.10, lib_cut = 1000, 
                                       group = 'metabolite', struc_zero = TRUE))

out.lum.PS <- suppressWarnings(ancombc(phyloseq = ps %>% 
                                          subset_samples(metabolite != '4hydroxybenzoic_acid') %>%
                                          subset_samples(algae == 'Phaeodactylum spent medium'),
                                        tax_level = 'Species',
                                        formula = "time + metabolite", 
                                        p_adj_method = "fdr",alpha = 0.05,
                                        prv_cut = 0.10, lib_cut = 1000, 
                                        group = 'metabolite', struc_zero = TRUE))

out.lum.P <- suppressWarnings(ancombc(phyloseq = ps %>% 
                                         subset_samples(metabolite != '4hydroxybenzoic_acid') %>%
                                         subset_samples(algae == 'Phaeodactylum'), 
                                       tax_level = 'Species',
                                       formula = "time + metabolite", 
                                       p_adj_method = "fdr",alpha = 0.05,
                                       prv_cut = 0.10, lib_cut = 1000, 
                                       group = 'metabolite', struc_zero = TRUE))

# merge differential abundance analysis results
DA.lum.all <- merge(out.lum.all$res$diff_abn['metabolite1'], 
                     out.lum.all$res$q_val['metabolite1'], by=0) %>%
  rename('metabolite1.x' = 'all') %>%
  rename('metabolite1.y' = 'q_val.all') %>%
  column_to_rownames(var='Row.names') %>%
  merge(out.lum.all$res$lfc['metabolite1'], by=0) %>%
  rename('metabolite1' = 'LFC.all') %>%
  column_to_rownames(var='Row.names')

DA.lum.C <- merge(out.lum.C$res$diff_abn['metabolite1'], 
                   out.lum.C$res$q_val['metabolite1'], by=0) %>%
  rename('metabolite1.x' = 'C') %>%
  rename('metabolite1.y' = 'q_val.C') %>%
  column_to_rownames(var='Row.names') %>%
  merge(out.lum.C$res$lfc['metabolite1'], by=0) %>%
  rename('metabolite1' = 'LFC.C') %>%
  column_to_rownames(var='Row.names')

DA.lum.PS <- merge(out.lum.PS$res$diff_abn['metabolite1'], 
                    out.lum.PS$res$q_val['metabolite1'], by=0) %>%
  rename('metabolite1.x' = 'PS') %>%
  rename('metabolite1.y' = 'q_val.PS') %>%
  column_to_rownames(var='Row.names') %>%
  merge(out.lum.PS$res$lfc['metabolite1'], by=0) %>%
  rename('metabolite1' = 'LFC.PS') %>%
  column_to_rownames(var='Row.names')

DA.lum.P <- merge(out.lum.P$res$diff_abn['metabolite1'], 
                   out.lum.P$res$q_val['metabolite1'], by=0) %>%
  rename('metabolite1.x' = 'P') %>%
  rename('metabolite1.y' = 'q_val.P') %>%
  column_to_rownames(var='Row.names') %>%
  merge(out.lum.P$res$lfc['metabolite1'], by=0) %>%
  rename('metabolite1' = 'LFC.P') %>%
  column_to_rownames(var='Row.names')

means <- colMeans(otu_table(ps.ABC.RA))  # mean relative abundance (after
                                         # adjustment for compositionallity)

DA.lum <- merge(DA.lum.all, DA.lum.C, by=0, all=TRUE) %>%
  subset(select = c('Row.names','all','q_val.all','LFC.all',
                    'C','q_val.C','LFC.C')) %>%
  merge(DA.lum.PS, by.x=1, by.y=0, all=TRUE) %>%
  subset(select = c('Row.names', 'all','q_val.all','LFC.all',
                    'C','q_val.C','LFC.C', 
                    'PS','q_val.PS','LFC.PS')) %>%
  merge(DA.lum.P, by.x=1, by.y=0, all=TRUE) %>%
  subset(select = c('Row.names', 'all','q_val.all','LFC.all',
                    'C','q_val.C','LFC.C', 
                    'PS','q_val.PS','LFC.PS',
                    'P','q_val.P','LFC.P')) %>%
  merge(as.data.frame(means), by.x=1,by.y=0, all=TRUE) %>%
  merge(as.data.frame(tax_table(ps.ABC.RA)), all=TRUE, by.x=1, by.y=0) %>%
  column_to_rownames(var='Row.names')

# save the results for plotting
write.csv(DA.lum, file.path(path, 'DA_lum.csv'))
