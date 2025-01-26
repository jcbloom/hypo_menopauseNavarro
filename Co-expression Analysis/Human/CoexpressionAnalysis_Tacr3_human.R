library(dplyr)
library(tibble)
library(matrixStats)
library(stringr) 
library(reshape2)
library(readr)
library(Hmisc)

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut],
    stringsAsFactors = FALSE
  )
}

flattenCorrMatrix_mutRank <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    mut_Rank = (cormat)[ut],
    stringsAsFactors = FALSE
  )
}

loess_data <- read.delim("LOESS_output_HUMAN.txt", header=TRUE, check.names=FALSE)

loess_data2 <- as.data.frame(t(loess_data), header=FALSE)
colnames(loess_data2) <- loess_data2[1,]
# Give a name to the row names column
loess_data2 <- rownames_to_column(loess_data2, "gene")

loess_norm_matrix <- t(as.matrix(loess_data2))
col_names <- as.character(loess_norm_matrix[1, ])
# Set the column names of the data
colnames(loess_norm_matrix) <- col_names
# Remove the first row (which is now the column names)
loess_norm_matrix <- loess_norm_matrix[-1, ]
# Set row names to sequential numbers
rownames(loess_norm_matrix) <- 1:nrow(loess_norm_matrix)

loess_corr <- rcorr(as.matrix(loess_norm_matrix))
loess_corr_flat <- flattenCorrMatrix(loess_corr$r, loess_corr$P)
loess_corr_flat$padj <- p.adjust(p = loess_corr_flat$p, method="BH")
loess_corr_flat$gene_pair <- paste(loess_corr_flat$row, loess_corr_flat$column, sep="_")

#To get the mutual rank - 
corr_ranks <- loess_corr$r

##Gene 1: Get order of correlation with each other gene --> go across rows. 
for(gene in colnames(loess_corr$r)){
  gene <- "TACR3"
  #Get the flattened table for the gene of interest
  loess_corr_flat_gene <- loess_corr_flat[loess_corr_flat$row == gene | loess_corr_flat$column == gene,]
  #Add row for the gene in question
  geneData <- data.frame(gene,gene,0,1,1,paste(gene,gene,sep="_"), stringsAsFactors = FALSE)
  colnames(geneData) <- colnames(loess_corr_flat_gene)
  loess_corr_flat_gene <- rbind(loess_corr_flat_gene, geneData)
  #head(loess_corr_flat_gene, 15)
  loess_corr_flat_gene$rank <- 15883 - rank(abs(loess_corr_flat_gene$cor))
  
  #get comparison gene
  loess_corr_flat_gene$compGene <- "NA"
  loess_corr_flat_gene[loess_corr_flat_gene$row == gene, "compGene"] <- loess_corr_flat_gene[loess_corr_flat_gene$row == gene, "column"]
  loess_corr_flat_gene[loess_corr_flat_gene$column == gene, "compGene"] <- loess_corr_flat_gene[loess_corr_flat_gene$column == gene, "row"]
  rownames(loess_corr_flat_gene) <- loess_corr_flat_gene$compGene
  
  corr_ranks[,gene] <- loess_corr_flat_gene[rownames(corr_ranks),"rank"]
}

mutual_ranks <- corr_ranks

write.table(x=loess_corr_flat_gene, file="Tacr3_coexp_human.txt",quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
