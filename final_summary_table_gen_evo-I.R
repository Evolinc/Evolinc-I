#!/usr/bin/Rscript

# Install dependencies
library(Biostrings)
library(splitstackshape)
library(dplyr)
library(getopt)


args<-commandArgs(TRUE)

#############################################################
## Command-line preparations and getopt parsing ##
#############################################################

options<-matrix(c('lincRNA',	'l',	1,	"character",
		  		  'lincRNAbed',	'b',	1,	"character",
		  		  'overlap',	'v',	2,	"character",
		  		  'tss',		't',	2,	"character",
		  		  'help',		'h', 	0,  "logical"),
		  		   ncol=4,byrow=TRUE)

ret.opts<-getopt(options,args)

if ( !is.null(ret.opts$help) ) {
  cat(getopt(options, usage=TRUE));
  q(status=1);
}

# No TSS file here
if(is.null(ret.opts$tss)) {
    # All lincRNA's
    fastaFile <- readDNAStringSet(ret.opts$lincRNA)
    lincRNA_ID = names(fastaFile)
    sequence = paste(fastaFile)
    size = width(fastaFile)
    df <- data.frame(lincRNA_ID, size)
    
    # Overlapping lincRNA's
    fastaFile_over <- readDNAStringSet(ret.opts$overlap)
    lincRNA_ID_over <- names(fastaFile_over)
    len <- length(lincRNA_ID_over)
    newmat <- matrix(lincRNA_ID_over, ncol =1 , nrow = len, byrow = T)
    overlap <- cSplit(as.data.table(newmat), "V1", "_o")
    
    # Finding out the lincRNA gene id:
    info = file.info("./intersect_output.txt")
    if (info$size == 0) {
        data4 <- data.frame(gene_id = rep("NA",nrow(df)))}
    else {    
        data1 <- read.table("./intersect_output.txt", sep="\t")
        names(data1) <- c("ID","gene")
        data1$gene <- as.character(data1$gene)
        gene2 <- unlist(strsplit(data1$gene, "_"))
        data2 <- matrix(gene2, ncol = 4, byrow = T)
        data3 <- as.data.frame(cbind(data1[1], data2[,2]))
        names(data3) <- c("id", "gene")
        data4 <- as.data.frame(data3 %>% group_by(id) %>% summarise(gene[1]))
        names(data4)[2] <- "gene_id" 
    }

    # Bed file
    bed_File <- read.table(ret.opts$lincRNAbed)
    t <- paste0(bed_File$V14 , ".gene=", bed_File$V11)
    bind <- as.data.frame(cbind(t,bed_File$V17))
    bind$t <- as.character(bind$t)
    bind$V2 <- as.character(bind$V2)
    bedfinal <- bind %>% group_by(t) %>% summarise(V2 = max(V2))
    ##Non present file
    cols <- ncol(bedfinal)
    cols <- cols + 1 
    bedfinal[,cols] <- NA
    
    # Merging All three files
    merge1 <- merge(x = df, y = overlap, by = 1, all = TRUE)
    if (info$size == 0) {
        merge2 <- cbind(merge1, data4) }
    else {
        merge2 <- merge(x = merge1, y = data4, by = 1, all = TRUE)
    }
    merge2$V1_2 <- as.character(merge2$V1_2)
    merge2$gene_id <- as.character(merge2$gene_id)
    merge2$V1_2[!is.na(merge2$V1_2)] <- "Yes"
    merge2$V1_2[is.na(merge2$V1_2)] <- "No"            
    colnames(merge2)[2] <- "Size(bp)"
    colnames(merge2)[3] <- "Overlapping_known_lincRNA"
    merge2[is.na(merge2)] <- "NA"
    merge3 <- merge(x = merge2, y = bedfinal, by = 1, all = TRUE)
    colnames(merge3)[5] <- "Number_of_exons"
    colnames(merge3)[6] <- "Has_TSS_data"
    merge3 <- merge3[,c(1:4,6,5)]
    # Writing data
    write.table(merge3, file = "final_Summary_table_evolinc-I.tsv", row.names = F, col.names = T, quote = F, sep = "\t")

# No Overlap file
} else if(is.null(ret.opts$overlap)) {
	# All lincRNA's
    fastaFile <- readDNAStringSet(ret.opts$lincRNA)
    lincRNA_ID = names(fastaFile)
    sequence = paste(fastaFile)
    size = width(fastaFile)
    df <- data.frame(lincRNA_ID, size)
    
    # CAGE supported lincRNA's
    fastaFile_Cage <- readDNAStringSet("lincRNAs.with.CAGE.support.annotated.fa")
    lincRNA_ID_Cage <- names(fastaFile_Cage)
    len1 <- length(lincRNA_ID_Cage)
    newmat1 <- matrix(lincRNA_ID_Cage, ncol =1 , nrow = len1, byrow = T)
    cage <- cSplit(as.data.table(newmat1), "V1", "_C")
    
    # Bed file
    bed_File <- read.table(ret.opts$lincRNAbed)
    t <- paste0(bed_File$V14 , ".gene=", bed_File$V11)
    bind <- as.data.frame(cbind(t,bed_File$V17))
    bind$t <- as.character(bind$t)
    bind$V2 <- as.character(bind$V2)
    bedfinal <- bind %>% group_by(t) %>% summarise(V2 = max(V2))
    #non present columns
    cols <- ncol(bedfinal)
    cols <- cols + 1 
    bedfinal[,cols] <- NA
    cols <- cols + 1 
    bedfinal[,cols] <- NA
    
    # Merging All three files
    merge1 <- merge(x = df, y = cage, by = 1, all = TRUE)
    merge2 <- merge(x = merge1, y = bedfinal, by = 1, all = TRUE)
    merge2$V1_2 <- as.character(merge2$V1_2)
    merge2$V1_2[is.na(merge2$V1_2)] <- "No"
    merge2$V1_2 <- sub("AGE_PLUS", "Yes", merge2$V1_2)
    merge2$V1_2 <- as.factor(merge2$V1_2)
    colnames(merge2)[2] <- "Size(bp)"
    colnames(merge2)[3] <- "Has_TSS_data"
    colnames(merge2)[4] <- "Number_of_exons"
    colnames(merge2)[5] <- "Overlapping_known_lincRNA"
    colnames(merge2)[6] <- "gene_id"
    merge2[,order(colnames(merge2))]
    merge2 <- merge2[,c(1:2,5:6,3:4)]
    # Writing data
    write.table(merge2, file = "final_Summary_table_evolinc-I.tsv", row.names = F, col.names = T, quote = F, sep = "\t")

# All present
} else {
    # All lincRNA's
    fastaFile <- readDNAStringSet(ret.opts$lincRNA)
    lincRNA_ID = names(fastaFile)
    sequence = paste(fastaFile)
    size = width(fastaFile)
    df <- data.frame(lincRNA_ID, size)
    
    # Overlapping lincRNA's
    fastaFile_over <- readDNAStringSet(ret.opts$overlap)
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
  
    # Bed file
    bed_File <- read.table(ret.opts$lincRNAbed)
    t <- paste0(bed_File$V14 , ".gene=", bed_File$V11)
    bind <- as.data.frame(cbind(t,bed_File$V17))
    bind$t <- as.character(bind$t)
    bind$V2 <- as.character(bind$V2)
    bedfinal <- bind %>% group_by(t) %>% summarise(V2 = max(V2))

    # Finding out the lincRNA gene id:
    info = file.info("./intersect_output.txt")
    if (info$size == 0) {
        data4 <- data.frame(gene_id = rep("NA",nrow(df)))}
    else {    
        data1 <- read.table("./intersect_output.txt", sep="\t")
        names(data1) <- c("ID","gene")
        data1$gene <- as.character(data1$gene)
        gene2 <- unlist(strsplit(data1$gene, "_"))
        data2 <- matrix(gene2, ncol = 4, byrow = T)
        data3 <- as.data.frame(cbind(data1[1], data2[,2]))
        names(data3) <- c("id", "gene")
        data4 <- as.data.frame(data3 %>% group_by(id) %>% summarise(gene[1]))
        names(data4)[2] <- "gene_id" 
    } 
    
    # Merging All three files
    merge1 <- merge(x = df, y = overlap, by = 1, all = TRUE)
    if (info$size == 0) {
        merge2 <- cbind(merge1, data4) }
    else {
        merge2 <- merge(x = merge1, y = data4, by = 1, all = TRUE)
    }
    merge2$V1_2.x <- as.character(merge2$V1_2.x)
    merge2$V1_2.y <- as.character(merge2$V1_2.y)
    merge2$V1_2.x[!is.na(merge2$V1_2.x)] <- "Yes"
    merge2$V1_2.y[!is.na(merge2$V1_2.y)] <- "Yes"
    merge2[is.na(merge2)] <- "No"
    colnames(merge2)[2] <- "Size(bp)"
    colnames(merge2)[3] <- "Overlapping_known_lincRNA"
    colnames(merge2)[4] <- "Has_TSS_data"
    merge3 <- merge(x = merge2, y = bedfinal, by = 1, all = TRUE)
    colnames(merge3)[5] <- "Number_of_exons"
    merge4 <- merge(merge3, data4, by=1, all = TRUE)
    merge4 <- merge4[,c(1:3,6,4,5)]
    merge4$gene <- as.character(merge4$gene)
    merge4$gene[is.na(merge4$gene)] <- "NA"
    merge4[,order(colnames(merge4))]
    
    # Writing data
    write.table(merge4, file = "final_Summary_table_evolinc-I.tsv", row.names = F, col.names = T, quote = F, sep = "\t")
}