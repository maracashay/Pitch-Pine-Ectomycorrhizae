
library(dada2)
#packageVersion("dada2")

##### 1. Set working directory ####

setwd("~/Downloads/Berts_ITS_Data/demultiplexed_fastq")

path<- setwd("~/Downloads/Berts_ITS_Data/demultiplexed_fastq")

#### 2. Sort the forward and reverse reads #####
fnFs <- sort(list.files(path, pattern=".fastq", full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), ".fastq"), `[`, 1)

#### 3. Examine the quality profiles ####

## forward reads ##
quartz()
plotQualityProfile(fnFs[1:6])

#### 4. Assign Filtering samples ####
filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))

out.2<- filterAndTrim(fnFs, filtFs, truncLen=400, 
                      minLen = 200, maxN=0, maxEE=5, truncQ = 2, 
                      multithread=TRUE, compress=TRUE)


#### 5. Learn the Error Rates ####

errf<- learnErrors(filtFs, multithread=TRUE) #set multithread=FALSE if on windows

plotErrors(errf, nominalQ = TRUE)

#### 6. Dereplicate the sequences ####

derepFs <- derepFastq(filtFs, verbose=TRUE)

names(derepFs) <- sample.names

##### 7. Call ESVs ####

dadaFs<-dada(derepFs, err=errf, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32, USE_QUALS=TRUE, multithread=TRUE)

dadaFs[[2]]

seqtab<-makeSequenceTable(dadaFs)
dim(seqtab)


#### 8. Remove Chimeric Sequences ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))

track <- cbind(out.2, sapply(dadaFs, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)

##### 9. Assign Taxonomy with Unite ####
taxa <- assignTaxonomy(seqtab.nochim, "/Users/admin/Downloads/Unite/sh_general_release_dynamic_01.12.2017.fasta", multithread=TRUE)

#### 10. Format data for analysis using Phyloseq object ####

library(phyloseq)
library(tibble)

Map<-import_qiime_sample_data("~/Box Sync/Berts_ITS_Data/BertsMetadata-1.txt")

esv.table<-otu_table(seqtab.nochim.mod, taxa_are_rows=FALSE)

# build phyloseq object
ps <- phyloseq(esv.table, tax_table(taxa),Map)

#Now, let's make sure that the correct information is included in our phyloseq object

ps@otu_table #Should include ESV info
ps@tax_table #Should include taxonomic info, but we get a NULL output
ps@sam_data

ps.1 <- ps
# exchange sequences with ASVs
dna <- Biostrings::DNAStringSet(taxa_names(ps.1))
names(dna) <- taxa_names(ps.1)
ps.1 <- merge_phyloseq(ps.1, dna)
taxa_names(ps.1) <- paste0("ASV", seq(ntaxa(ps.1)))
ps.1


#### 11. General processing and analysis using silva phyloseq object####

sample_sums(ps.1)
# the sample with the smallest sequence count is E28 with 6, so let's remove samples with fewer than 2,000 seqs
ps.pruned <- prune_samples(sample_sums(ps.1)>=1400, ps)
sample_sums(ps.pruned)


full.phyloseq = subset_taxa(ps.pruned, Kingdom %in% "k__Fungi")
