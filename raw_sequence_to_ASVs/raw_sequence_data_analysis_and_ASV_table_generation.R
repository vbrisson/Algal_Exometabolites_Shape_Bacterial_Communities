# DADA2 ANALYSIS

# Author: Vanessa Brisson, PhD
# Email: brisson2@llnl.gov
# Date Created: 2023-01-11

# based on the tutorial at https://benjjneb.github.io/dada2/tutorial.html

# ------------------------------------------------------------------------------
# load DADA2 package
library(dada2)

# You will need the sequencing data, which is available through the NCBI 
# Sequence Read Archive under BioProject number PRJNA803592:
#   https://www.ncbi.nlm.nih.gov/bioproject/PRJNA803592/
#
# You will also need the taxonomy reference databases, which are available at:
#   https://benjjneb.github.io/dada2/training.html
#
#   The referencec databases used for the paper were:
#      silva_nr99_v138.1_train_set.fa
#      rdp_train_set_18.fa


# change this path variable to wherever the sequencing data are stored
path <- 'C:/Users/username/Documents'

# get file names
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# assess quality profiles
plotQualityProfile(fnFs[51:55])
plotQualityProfile(fnRs[51:55])
# both forward an revers reads look good so trim both to 140

# filter based on quality and trim to 140 bases
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(140,140),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)

# learn error rates
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# sample inference
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
dadaFs[[1]]
dadaRs[[1]]

# merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:255]
dim(seqtab2)
table(nchar(getSequences(seqtab2)))

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# track reads through different steps
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, file.path(path,"tax/silva_nr99_v138.1_train_set.fa.gz"), multithread=FALSE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

taxa2 <- assignTaxonomy(seqtab.nochim, file.path(path,"tax/rdp_train_set_18.fa.gz"), multithread=FALSE)
taxa2.print <- taxa2 # Removing sequence rownames for display only
rownames(taxa2.print) <- NULL
head(taxa2.print)

saveRDS(seqtab.nochim, file = "C:Users/username/Documents/ASVtable.rds")
saveRDS(taxa, file = "C:/Users/username/Documents/taxtable.rds")
saveRDS(taxa2, file = "C:/Users/username/Documents/taxtable_RDP.rds")
