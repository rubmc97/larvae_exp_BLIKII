library(dada2)

path = "/project/genomics/Illumina/rubenmc_larvae/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

forwplot<-ggplotly(plotQualityProfile(fnFs[1:length(fnRs)], aggregate = TRUE) + geom_hline(yintercept=c(15,25,35), color=c("red","blue","green"), size=0.5), width =750)
forwplot

revwplot<-ggplotly(plotQualityProfile(fnRs[1:length(fnRs)], aggregate = TRUE) + geom_hline(yintercept=c(15,25,35), color=c("red","blue","green"), size=0.5), width =750)
revwplot

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(20, 16), truncLen=c(255,190),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE (only needed for filterAndTrim)
#If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse reads (eg. maxEE=c(2,5)), and reducing the truncLen to remove low quality tails.
out

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
#Forward
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#Reverse

dadaFs[[1]]
#Inspect object

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 329:351]
#Depending on the values of table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, file = "track_table_larvae.csv")
#save the table resuming the reads along the filtering steps

taxa = assignTaxonomy(seqtab.nochim, "silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)
taxa = addSpecies(taxa, "silva_v138.2_assignSpecies.fa.gz")

samdf = read.csv("meta_larvae.csv")
samdf = samdf %>% column_to_rownames("SampleID")

phy = phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(taxa))

saveRDS(object = phy, file = "phyloseq_larvae_raw.rds")