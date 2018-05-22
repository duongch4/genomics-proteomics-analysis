#### Load libraries and Source other files if needed ####
sourceFunctions <- "C:/Users/DBC/OneDrive/R_myFunctions"
source(paste(sourceFunctions, "Graphs.R", sep = "/"))
source(paste(sourceFunctions, "gplots_Dendrogram.R", sep = "/"))

# Load libraries
library(edgeR)
library(data.table)
library(tidyr)
#library(HTSFilter) #when there are technical replicates for each condition
#library(gridExtra) # Arrange ggplot2
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")

#### Setup ####
setwd("C:/Users/DBC/OneDrive/Courses/Summer_2017/RNASeq")

# # Import HTSeq count data for each cell line and treatment ####
# M_G561_ctrl <- read.table("G561_ctrlAligned_reverse_HTseq_count.txt",stringsAsFactors=F, row.names=1, header=F)
# dim(M_G561_ctrl)
# class(M_G561_ctrl[,1])
# 
# M_G561_drug=read.table("G561_drugAligned_reverse_HTseq_count.txt",stringsAsFactors=F, row.names=1, header=F)
# dim(M_G561_drug)
# class(M_G561_drug[,1])
# 
# 
# M_G564_ctrl=read.table("G564_ctrlAligned_reverse_HTseq_count.txt",stringsAsFactors=F, row.names=1, header=F)
# dim(M_G564_ctrl)
# class(M_G564_ctrl[,1])
# 
# M_G564_drug=read.table("G564_drugAligned_reverse_HTseq_count.txt",stringsAsFactors=F, row.names=1, header=F)
# dim(M_G564_drug)
# class(M_G564_drug[,1])
# 
# 
# M_G583_ctrl=read.table("G583_ctrlAligned_reverse_HTseq_count.txt",stringsAsFactors=F, row.names=1, header=F)
# dim(M_G583_ctrl)
# class(M_G583_ctrl[,1])
# 
# M_G583_drug=read.table("G583_drugAligned_reverse_HTseq_count.txt",stringsAsFactors=F, row.names=1, header=F)
# dim(M_G583_drug)
# class(M_G583_drug[,1])
# 
# # Check that each dataset has the same rows (Quality checking) ####
# table(rownames(M_G561_ctrl)==rownames(M_G561_drug))
# table(rownames(M_G561_ctrl)==rownames(M_G564_ctrl))
# table(rownames(M_G561_ctrl)==rownames(M_G564_drug))
# table(rownames(M_G561_ctrl)==rownames(M_G583_ctrl))
# table(rownames(M_G561_ctrl)==rownames(M_G583_drug))
# 
# # The last 5 rows give the summary for each data set ####
# tail(rownames(M_G561_ctrl))
# 
# ## > tail(rownames(M_G561_ctrl))
# ## [1] "ZZZ3"                   "__no_feature"           "__ambiguous"           
# ## [4] "__too_low_aQual"        "__not_aligned"          "__alignment_not_unique"
# 
# # Extract the summary info into one matrix
# M_info_count=cbind("G561_ctrl"=M_G561_ctrl[26486:26490,1],"G564_ctrl"=M_G564_ctrl[26486:26490,1],"G583_ctrl"=M_G583_ctrl[26486:26490,1],"G561_drug"=M_G561_drug[26486:26490,1],"G564_drug"=M_G564_drug[26486:26490,1],"G583_drug"=M_G583_drug[26486:26490,1])
# rownames(M_info_count)=rownames(M_G561_ctrl)[26486:26490]
# head(M_info_count)
# 
# # Extract the actual counts data into one matrix ####
# M_count=cbind("G561_ctrl"=M_G561_ctrl[1:26485,1],"G564_ctrl"=M_G564_ctrl[1:26485,1],"G583_ctrl"=M_G583_ctrl[1:26485,1],"G561_drug"=M_G561_drug[1:26485,1],"G564_drug"=M_G564_drug[1:26485,1],"G583_drug"=M_G583_drug[1:26485,1])
# rownames(M_count) <- rownames(M_G561_ctrl)[1:26485]
load("M_count.RData") ####
class(M_count)
dim(M_count)
names(M_count)

# Create the DGEList (dge) object ####
dge <- DGEList(counts = M_count[, .SD, .SDcols = grep("G", names(M_count))], genes = M_count[["rn"]])

# Filtering
# Keep genes with more than 1 CPM in at least 3 samples
keep <- rowSums(cpm(dge) > 1) >= 3
addmargins(table(keep))

dge <- dge[keep, , keep.lib.sizes = FALSE] # Change the lib.sizes accordingly

# TMM (trimmed mean of M-values, where M = log2(ratio)) normalisation
dge <- calcNormFactors(dge, method = "TMM")
dge$samples

fwrite(x = as.data.table(dge$genes), file = "library_rnaseq.csv", row.names = FALSE, col.names = FALSE)

# Normalised cpm
tabCount <- as.data.table(cpm(dge[["counts"]], prior.count = 2, log = TRUE))
tabCount[, genes := dge[["genes"]]]

# Design matrix
Line <- sapply(strsplit(names(tabCount)[grep(pattern = paste("G561", "G564", "G583", sep = "|"), x = names(tabCount))], split = "_"), "[[", 1)
Condition <- sapply(strsplit(names(tabCount)[grep(pattern = paste("G561", "G564", "G583", sep = "|"), x = names(tabCount))], split = "_"), "[[", 2)

Line <- factor(Line)
Condition <- factor(Condition, levels = c("SGC2096", "GSK591"))

design <- model.matrix(~ Line + Condition)
rownames(design) <- colnames(dge)

design # Check for the coefficients meaning

# Estimate Dispersion (all of common, trend, and tagwise)
dge <- estimateDisp(dge, design, robust = TRUE)

dge$common.dispersion
sqrt(dge$common.dispersion)

# Mean-Var plot 
plotMeanVar(dge, show.tagwise.vars = TRUE, NBline = TRUE, main = "Variance vs Mean", nbin = 100)
legend.txt <- c("Red cross: average raw variance", 
                "Blue dot: genewise variance", "Blue line: expected NegBin M-V with the common dispersion",
                "Black line: expected Poisson M-V")
legend("topleft", legend = legend.txt, text.col = c("red", "blue", "blue", "black"), cex = 0.7)

# BCV plot
plotBCV(dge, main = "BCV vs Abundance", cex = 0.4)

# ANALYSIS using Neg Bin model ####

# LRT
fitLRT <- glmFit(dge, design)
gof(fitLRT, plot = TRUE)
colnames(design)

# # QL
# fitQLF <- glmQLFit(dge, design, robust = TRUE)
# gof(fitQLF, plot = TRUE)
# plotQLDisp(fitQLF)

# Check if there was a need to adjust for the cell lines, use both LRT and QL approaches
# Test the DE btw 3 cell lines
# => a lot DE => yes there is a batch/block effect
LRT <- glmLRT(fitLRT, coef = 2:3)
topTags(LRT)
FDR <- p.adjust(LRT$table$PValue, method = "BH")
sum(FDR < 0.05)

# QLF <- glmQLFTest(fitQLF, coef = 2:3)
# topTags(QLF)
# FDR <- p.adjust(QLF$table$PValue, method = "BH")
# sum(FDR < 0.05)

# Treatment effect, should use independent filtering (HTSFilter)
# when there are technical replicates in each condition => After LRT <- glmLRT(fitLRT, coef = 4)
# LRT <- HTSFilter(LRT, DGEGLM = fitLRT, s.len = 25, plot = FALSE)$filteredData
LRT <- glmLRT(fitLRT, coef = 4)
result <- setDT(topTags(LRT, n = nrow(dge))$table, key = "genes")

# Join result with tabCount
result <- result[tabCount, on = "genes", nomatch = 0L]
setkeyv(result, cols = "FDR")

# SAVE ####
saveRDS(result, file = "result_dge.rds")
saveRDS(tabCount, file = "tabCount_dge.rds")

sum(result$FDR < 0.05)

# QLF <- glmQLFTest(fitQLF, coef = 4)
# tabDE_QLF <- topTags(QLF, n = nrow(dge))$table
# sum(tabDE_QLF$FDR < 0.05)

# Summary of up/down-regulated DE genes at 5% FDR
summary(decideTests(LRT, p.value = 0.05))

#summary(decideTests(QLF, p.value = 0.05))

# # Plot all LFCs vs Average count size
# # Darkgreen lines indicate 2-fold up or down
# plotMD(LRT, cex = 0.4, main = "Drug vs Control")
# abline(h = c(-1,1), col = "darkgreen")
# 
# plotSmear(LRT, cex = 0.4, de.tags = DE, main = "Drug vs Control")
# abline(h = c(-1,1), col = "darkgreen")
# 
# plotMD(QLF, cex = 0.4, main = "Drug vs Control")
# abline(h = c(-1,1), col = "darkgreen")

# # Comparison with DESeq2 ####
# edgeR <- topTags(LRT, n = nrow(dge))$table
# DESeq2 <- as.data.frame(res05)
# 
# addmargins(table(sig.edgeR = edgeR$FDR < 0.05, sig.DESeq2 = DESeq2$padj < 0.05))
# 
# tab <- merge(x = edgeR, y = DESeq2, by = "row.names")
# 
# with(tab, 
#      plot(logFC, log2FoldChange, pch = 20, col = "black",
#           xlab = "LFC egdeR", ylab = "LFC DESeq2",
#           main = "LFC for edgeR vs DESeq2"))
# 
# with(subset(tab, FDR < 0.05), 
#      points(logFC, log2FoldChange, pch = 20, col = "red"))
# 
# with(subset(tab, padj < 0.05),
#      points(logFC, log2FoldChange, pch = 20, col = "green"))
# 
# legend("topleft", xjust = 1, yjust = 1, pch = 20,
#        legend = c("FDR < 0.05 edgeR only", "FDR < 0.05 both", "FDR > 0.05"),
#        col = c("red", "green", "black"), cex = 0.7, text.width = 0.7 )
# 
# edgeR_DE <- row.names(edgeR[which(edgeR$FDR < 0.05), ])
# DESeq2_DE <- row.names(DESeq2[which(DESeq2$padj < 0.05), ])
# 
# diff <- setdiff(edgeR_DE, DESeq2_DE)
# write.csv(diff, "DEGenes_SetDiff_edgeR_DESeq2.csv")
