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
