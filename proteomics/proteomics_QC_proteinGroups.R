#### Load libraries and Source other files if needed ####
sourceFunctions <- "C:/Users/DBC/OneDrive/R_myFunctions"
source(paste(sourceFunctions, "Graphs.R", sep = "/"))

library(data.table)
library(tidyr) # for separate()
library(FactoMineR) # for PCA
library(factoextra) # for PCA visualisation

#### STORAGE ####
setwd("C:/Users/DBC/OneDrive/Courses/Summer_2017/Proteomics")

dataType <- "proteomics"
folderName <- paste(dataType, "QualityControl", sep = "_")
currentStorage <- paste(getwd(), folderName, sep = "/")
dir.create(path = currentStorage, showWarnings = TRUE)

### Load peptides data ####
peptides_norm <- readRDS("peptides_norm.rds")

## Use melt() from reshape2
peptides_plot <- melt(data = peptides_norm,
                      id.vars = "Proteins", measure.vars = grep(pattern = paste("G561", "G564", "G583", sep = "|"), x = names(peptides_norm)),
                      variable.name = "Sample", value.name = "Intensity")

peptides_plot <- separate(data = peptides_plot,
                          col = Sample, into = c("Cell_Line", "Treatment", "Replicate"), sep = "_",
                          remove = FALSE, extra = "merge")

### Plot density distribution ####
keep <- !(is.na(peptides_plot$Intensity))
table(keep)
peptides_plot <- peptides_plot[keep, ] 

Num1_Cat3(data = peptides_plot[, .SD, .SDcols = c("Cell_Line", "Intensity", "Treatment", "Replicate")],
          PlotTitle = "Overall experiment", width_in = 5, height_in = 4,
          filename = paste(currentStorage,
                           paste("OverallPlot", dataType, sep = "_"),
                           sep = "/"))

Boxplot2(data = peptides_plot[, .SD, .SDcols = c("Sample", "Intensity", "Treatment")],
         PlotTitle = "Boxplot for Intensity",
         filename = paste(currentStorage,
                          paste("BoxplotIntensity", dataType, sep = "_"),
                          sep = "/"))

OverlaidHist(data = peptides_plot[, .SD, .SDcols = c("Intensity", "Treatment")], showDensity = TRUE,
             Binwidth = 0.2, seq_break_start = 15, seq_break_end = 25, seq_break_by = 5,
             PlotTitle = "Histogram for each treatment", width_in = 5, height_in = 4,
             filename = paste(currentStorage,
                              paste("DistributionIntensity", dataType, sep = "_"),
                              sep = "/"))


### MDS plot ####
plotMDS(data = peptides_norm[,5:ncol(peptides_norm)], distMetric = "euclidean", sepCols = c("Line", "Treatment", "Replicate"), sep = "_",
        PlotTitle = "MDS plot - Proteomics", width_in = 5, height_in = 4,
        filename = paste(currentStorage,
                         paste("MDS", dataType, sep = "_"),
                         sep = "/"))

### PCA ####
# Imputation
# nbdim <- estim_ncpPCA(peptides_norm[,5:ncol(peptides_norm)], ncp.min = 0, ncp.max = 5, method.cv = "Kfold", nbsim = 50)
# imputed <- imputePCA(peptides_norm[,5:ncol(peptides_norm)], ncp = nbdim$ncp)
# saveRDS(nbdim, file = "estim_ncpPCA.rds")
# saveRDS(imputed, file = "imputePCA.rds")
# pca_result <- PCA(res.comp$completeObs, graph = FALSE)

pca_no_na <- na.omit(peptides_norm[,5:ncol(peptides_norm)])
pca_result <- PCA(pca_no_na, graph = FALSE)

# Get eigenvalues / variance
get_eig(pca_result)
# Visualise eigenvalues / variance
fviz_screeplot(pca_result, addlabels = TRUE)
ggsave(filename = paste(paste(currentStorage,
                        paste("PCA_Scree", dataType, sep = "_"),
                        sep = "/"), "png", sep = "."),
       width = 10, height = 7, units = "in")

# Contributions of variables to PC1
fviz_contrib(pca_result, choice = "var", axes = 1, top = 20)
ggsave(filename = paste(paste(currentStorage,
                              paste("PC1_Contribution", dataType, sep = "_"),
                              sep = "/"), "png", sep = "."),
       width = 10, height = 7, units = "in")

# Contributions of variables to PC2
fviz_contrib(pca_result, choice = "var", axes = 2, top = 10)
ggsave(filename = paste(paste(currentStorage,
                              paste("PC2_Contribution", dataType, sep = "_"),
                              sep = "/"), "png", sep = "."),
       width = 10, height = 7, units = "in")

# Control variable colors using their contributions
#Length of the loading arrows approximates the standard deviation of original variables 
#(squared length approximates variance), scalar products between any two arrows approximate the covariance between them,
#and cosines of the angles between arrows approximate correlations between original variables
png(filename = paste(paste(currentStorage,
                           paste("PCA_Loadings", dataType, sep = "_"),
                           sep = "/"), "png", sep = "."),
    width = 12, height = 12, units = "in", res = 300)

grid.arrange(fviz_pca_var(pca_result, col.var = "contrib", geom = c("point", "text"), col.circle = "darkgreen",
                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                          repel = TRUE), 
             fviz_pca_var(pca_result, col.var="contrib", geom = c("arrow", "text"), col.circle = "darkgreen",
                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                          repel = TRUE),
             nrow = 2, ncol = 1)

dev.off()

# # Graph of individuals
# # 1. Use repel = TRUE to avoid overplotting
# # 2. Control automatically the color of individuals using the cos2
# # cos2 = the quality of the individuals on the factor map
# # Use points only
# # 3. Use gradient color
# fviz_pca_ind(pca_result, col.ind = "cos2", label = "none",
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE # Avoid text overlapping (slow if many points)
# )
# 
# # Biplot of individuals and variables
# fviz_pca_biplot(res.pca, repel = TRUE)