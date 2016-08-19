#!/usr/bin/Rscript

library(Biostrings)
library(splitstackshape)
library(dplyr)

# All lincRNA's
fastaFile <- readDNAStringSet("All.lincRNAs.fa")
lincRNA_ID = names(fastaFile)
sequence = paste(fastaFile)
size = width(fastaFile)
df <- data.frame(lincRNA_ID, size)

# Overlapping known lincRNA's
fastaFile_over <- readDNAStringSet("lincRNAs.overlapping.known.lincs.fa")
lincRNA_ID_over <- names(fastaFile_over)
len <- length(lincRNA_ID_over)
newmat <- matrix(lincRNA_ID_over, ncol =1 , nrow = len, byrow = T)
overlap <- cSplit(as.data.table(newmat), "V1", "_o")

# CAGE supported lincRNA's
fastaFile_Cage <- readDNAStringSet("lincRNAs.with.CAGE.support.annotated.fa")
lincRNA_ID_Cage <- names(fastaFile_Cage)
len1 <- length(lincRNA_ID_Cage)
newmat1 <- matrix(lincRNA_ID_Cage, ncol =1 , nrow = len1, byrow = T)
cage <- cSplit(as.data.table(newmat1), "V1", "_C")

# Merging All with overlapping known lincRNA's
merge1 <- merge(x = df, y = overlap, by = 1, all = TRUE)

# Merging above with CAGE supported lincRNA's
merge2 <- merge(x = merge1, y = cage, by = 1, all = TRUE)

merge2$V1_2.x <- as.character(merge2$V1_2.x)
merge2$V1_2.y <- as.character(merge2$V1_2.y)
merge2$V1_2.x[!is.na(merge2$V1_2.x)] <- "Yes"
merge2$V1_2.y[!is.na(merge2$V1_2.y)] <- "Yes"
merge2[is.na(merge2)] <- "No"

colnames(merge2)[2] <- "Size(bp)"
colnames(merge2)[3] <- "Overlapping_known_lincRNA"
colnames(merge2)[4] <- "Has_TSS_data"

# Bed file
bed_File <- read.table("lincRNA.bed")

t <- paste0(bed_File$V14 , ".gene=", bed_File$V11)
bind <- as.data.frame(cbind(t,bed_File$V17))
bind$t <- as.character(bind$t)
bind$V2 <- as.character(bind$V2)

bedfinal <- bind %>% group_by(t) %>% summarise(V2 = max(V2))

# Final merge
merge3 <- merge(x = merge2, y = bedfinal, by = 1, all = TRUE)
colnames(merge3)[5] <- "Number_of_exons"

# Writing data
write.table(merge3, file = "Final_Summary_table_evolinc-I.tsv", row.names = F, col.names = T, quote = F, sep = "\t")

