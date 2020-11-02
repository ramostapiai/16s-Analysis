library(dada2); packageVersion("dada2")
setwd("path")
path <- "path/sequences" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Perform filtering and trimming
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

#Filter the forward and reverse reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
                     maxN=0, maxEE=c(2,2), truncQ=0, rm.phix=TRUE,
                     compress=TRUE, multithread=16) # On Windows set multithread=FALSE
head(out)

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=16)
errR <- learnErrors(filtRs, multithread=16)

plotErrors(errF, nominalQ=TRUE)

#Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Sample Inference
#Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=16)
dadaRs <- dada(derepRs, err=errR, multithread=16)

#Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=16, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

saveRDS(object = seqtab.nochim, file = "path/seqtab.nochim.RDS")

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

#Save track table
write.table(track, "path/track.tsv", sep = "\t")

#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa", multithread=16)
taxa <- addSpecies(taxa, "silva_species_assignment_v132.fa")

#Let’s inspect the taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(object = taxa, file = "path/taxa.RDS")

###Construction the phylogenetic tree###

#HACER ALINEAMIENTO(DECIPHER)
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs
library(DECIPHER)
DNA <- DNAStringSet(seqs)
align<-AlignSeqs(DNA, processors = 16)
saveRDS(object = align, file = "path/align.RDS")
writeXStringSet(align, "path/align_seqs")

#Make phylogenetic tree with fasttree
system("~/miniconda3/bin/fasttree -gtr -nt align_seqs > tree_align_seqs")

library(phyloseq)
#Filogenia con FASTTREE
tree <- read_tree(treefile = "path/tree_align_seqs")
#Handoff to phyloseq
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa),
               phy_tree(tree))

saveRDS(object = ps, file = "path/phyloseq.RDS")


