#### LIBRARIES and SOURCES ####
library(gplots) # for heatmap.2()
library(STRINGdb)
#library(biomaRt)
library(wordcloud2)
library(webshot); #webshot::install_phantomjs()
library("htmlwidgets") # for Saving wordcloud2
library(data.table)

sourceFunctions <- "C:/Users/DBC/OneDrive/R_myFunctions"
source(paste(sourceFunctions, "Graphs.R", sep = "/"))
source(paste(sourceFunctions, "gplots_Dendrogram.R", sep = "/"))

#source("http://bioconductor.org/biocLite.R")

#biocLite("biomaRt")
#biocLite("STRINGdb")

############################################### LOAD DATA ####
getwd()
setwd("C:/Users/DBC/OneDrive/Courses/Summer_2017/RNASeq")

result <- readRDS(file = "result_dge.rds")

##########################################################

################################################### SETUP ####
dataType <- "rnaseq"
annoType <- "RNA-Seq"
contrast <- "GSK591_Vs_SGC2096"
cell <- "All cell lines"
condition1 <- "GSK591"
condition2 <- "SGC2096"
condition3 <- NULL

LFC_col <- "logFC"
pval_col <- "PValue"
qval_col <- "FDR"
Genes_col <- "genes"

# Set FDR threshold for DEGs and for Enrichment
fdr_DEG <- 0.05
fdr_enrichment <- 0.05

# Set size for word clouds
DEG_wordsize <- list(0.6, 0.65, 0.75)
DEG_word_width <- list(1.05e3, 1.05e3, 1.05e3)
DEG_word_height <- list(0.7e3, 0.7e3, 0.7e3)

GKPI_wordsize <- list(0.13, 0.8, 0.2, 0.12)
GKPI_word_width <- list(1.05e3, 1.4e3, 1.2e3, 1.30e3)
GKPI_word_height <- list(0.7e3, 0.5e3, 0.7e3, 0.78e3)

# Set size for gene names on heatmap (full DEGs list), and sample names
heatmap_row_wordsize <- 0.08
heatmap_col_wordsize <- 0.8

##########################################################

#### STORAGE ####
folderName <- paste(dataType, contrast, sep = "_")
currentStorage <- paste(getwd(), folderName, sep = "/")
dir.create(path = currentStorage, showWarnings = TRUE)

dir.create(path = paste(currentStorage, "Fig", sep = "/"), showWarnings = TRUE)

myFolders <- c("DEG", "QualityControl", "Volcano", "Heatmap", "ProteinNetwork_PPI", "GO", "KEGG", "PFAM", "InterPro")
for (i in 1:length(myFolders)) {dir.create(path = paste(currentStorage, "Fig", myFolders[i], sep = "/"), showWarnings = TRUE)}

GKPI_folders <- list("GO", "KEGG", "PFAM", "InterPro")

## Are the non-adjusted up-values okay => should be U-shape or right skew ####
Hist(data = subset(result, select = c(pval_col)), 
     Binwidth = 0.05, seq_break_start = 0.25, seq_break_end = 0.75, seq_break_by = 0.25, 
     PlotTitle = "Distribution of non-adjusted p-values", width_in = 5, height_in = 4,
     filename = paste(currentStorage, "Fig", "QualityControl", paste("Hist_pvalue", dataType, contrast, sep = "_"), sep = "/"))

## MDS ####
plotMDS(data = result[, .SD, .SDcols = grep(pattern = paste("G561", "G564", "G583", sep = "|"), x = names(result))],
        distMetric = "euclidean", sepCols = c("Line", "Treatment"), sep = "_",
        PlotTitle = "MDS plot : RNA-Seq", width_in = 5, height_in = 4,
        filename = paste(currentStorage, "Fig", "QualityControl", paste("MDS", dataType, contrast, sep = "_"), sep = "/"))

## DEG ####
DEG_folders <- list("DEG_All", "DEG_DownReg", "DEG_UpReg")
for (i in 1:length(DEG_folders)) {dir.create(path = paste(currentStorage, "Fig", "DEG", DEG_folders[[i]], sep = "/"), showWarnings = TRUE)}

DEG <- list(result[result[[qval_col]] < fdr_DEG],
            result[ (result[[qval_col]] < fdr_DEG) & (result[[LFC_col]] < 0) ],
            result[ (result[[qval_col]] < fdr_DEG) & (result[[LFC_col]] > 0) ])
names(DEG) <- DEG_folders

for (i in 1:length(DEG_folders)) {
  
  # Export .csv
  fwrite(x = DEG[[i]],
         file = paste(currentStorage, "Fig", "DEG", DEG_folders[[i]],
                      paste(paste(DEG_folders[[i]], dataType, contrast, paste("DEGfdr", fdr_DEG*100, sep = ""), sep = "_"),
                            "csv", sep = "."),
                      sep = "/"))
  
  # Word cloud
  DEG[[i]]$logFDR <- -log(DEG[[i]][[qval_col]])
  
  word <- wordcloud2(DEG[[i]][, .SD, .SDcols = c(Genes_col, "logFDR")], minRotation = 0, maxRotation = 0, size = DEG_wordsize[[i]])
  
  # Save in html
  saveWidget(widget = word,
             file = paste(currentStorage, "Fig", "DEG", DEG_folders[[i]],
                          paste(
                            paste(DEG_folders[[i]], "Word", dataType, contrast,
                                  paste("DEGfdr", fdr_DEG*100, sep = ""),
                                  sep = "_"),
                            "html",
                            sep = "."),
                          sep = "/"), 
             selfcontained = F)
  
  # Save in png
  webshot::webshot(url = paste(currentStorage, "Fig", "DEG", DEG_folders[[i]],
                               paste(
                                 paste(DEG_folders[[i]], "Word", dataType, contrast,
                                       paste("DEGfdr", fdr_DEG*100, sep = ""),
                                       sep = "_"),
                                 "html",
                                 sep = "."),
                               sep = "/"),
                   file = paste(currentStorage, "Fig", "DEG", DEG_folders[[i]],
                                paste(
                                  paste(DEG_folders[[i]], "Word", dataType, contrast,
                                        paste("DEGfdr", fdr_DEG*100, sep = ""),
                                        sep = "_"),
                                  "png",
                                  sep = "."),
                                sep = "/"),
                   delay = 20, vwidth = DEG_word_width[[i]], vheight = DEG_word_height[[i]])
}



## Volcano plot ####
# GSK vs SGC
Volcano(data = subset(result, select = c(LFC_col, qval_col, Genes_col)), 
        PlotTitle = paste(paste("Volcano", annoType, sep = " : "), paste(condition1, "vs", condition2), sep = " : "), 
        filename = paste(currentStorage, "Fig", "Volcano", paste("Volcano", dataType, contrast, sep = "_"), sep = "/"))

##########################################################

## Heatmap with package "gplots" ####
# Use makeDendrogram() function to compute a dendrogram 
# based on Euclidean / Pearson distance metric, then supply this in the 
# heatmap.2 function as the "Rowv" argument

DE <- copy(result)

table(is.na(DE))
sapply(DE, class)
names(DE)

# Get indices for cols
keep_ind <- grep(pattern = paste(Genes_col, qval_col, sep = "|"), x = colnames(DE), perl = TRUE)
condition1_ind <- grep(pattern = condition1, x = colnames(DE), perl = TRUE)
condition2_ind <- grep(pattern = condition2, x = colnames(DE), perl = TRUE)
if (is.null(condition3)) {
  message("No condition 3"); condition3_ind <- NULL
} else {condition3_ind <- grep(pattern = condition3, x = colnames(DE), perl = TRUE)}

DE <- DE[, .SD, .SDcols = c(keep_ind, condition1_ind, condition2_ind, condition3_ind)]

if (is.null(condition3)) {
  condition_ind <- grep(pattern = paste(condition1, condition2, sep = "|"), x = colnames(DE), perl = TRUE)
} else {condition_ind <- grep(pattern = paste(condition1, condition2, condition3, sep = "|"), x = colnames(DE), perl = TRUE)}

# Imputation by means for proteomics
if (dataType == "proteomics") {
  impute_mean(DE,
              myCols = condition_ind,
              number_of_skipped_columns = ncol(DE) - length(condition_ind)
  )
} else {DE <- na.omit(DE)}

table(is.na(DE))

# FDR < 0.05
DE <- DE[which(DE[[qval_col]] < fdr_DEG)]
names(DE)

# Full
myDendrogram <- makeDendrogram(DE[, .SD, .SDcols = c(condition_ind)], method = "euclidean")
png(paste(currentStorage, "Fig", "Heatmap", paste("Heatmap", dataType, contrast, "euclidean_full.png", sep = "_"), sep = "/"),
    width = 15, height = 15, units = "in", res = 1e3)
par(font = 2, font.lab = 2, font.axis = 2, cex.main = 2, cex.axis = 2, lwd = 2)
heatmap.2(as.matrix(DE[, .SD, .SDcols = c(condition_ind)]), scale = "row",
          Rowv = myDendrogram, Colv = F,
          dendrogram = "row",
          density.info = "none",
          trace = "none", cexRow = heatmap_row_wordsize, cexCol = heatmap_col_wordsize,
          labRow = DE[[Genes_col]], srtCol = +45, offsetRow = 5e-10,
          lhei = c(1,8), lwid = c(4,8), margins = c(10,20))
dev.off()

# Top 100 DE Genes
myDendrogram <- makeDendrogram(DE[1:100, .SD, .SDcols = c(condition_ind)], method = "euclidean")
png(paste(currentStorage, "Fig", "Heatmap", paste("Heatmap", dataType, contrast, "euclidean_top100.png", sep = "_"), sep = "/"),
    width = 15, height = 15, units = "in", res = 3e2)
par(font = 2, font.lab = 2, font.axis = 2, cex.main = 2, cex.axis = 2, lwd = 2)
heatmap.2(as.matrix(DE[1:100, .SD, .SDcols = c(condition_ind)]), scale = "row",
          Rowv = myDendrogram, Colv = F, 
          dendrogram = "row", 
          density.info = "none", 
          trace = "none", cexRow = 0.7, cexCol = heatmap_col_wordsize,
          labRow = DE[[Genes_col]], srtCol = +45, offsetRow = 5e-10,
          lhei = c(1,8), lwid = c(4,8), margins = c(10,25))
dev.off()

# Different distance metric - Pearson
myDendrogram <- makeDendrogram(DE[1:100, .SD, .SDcols = c(condition_ind)], method = "pearson")
png(paste(currentStorage, "Fig", "Heatmap", paste("Heatmap", dataType, contrast, "pearson_top100.png", sep = "_"), sep = "/"),
    width = 15, height = 15, units = "in", res = 3e2)
par(font = 2, font.lab = 2, font.axis = 2, cex.main = 2, cex.axis = 2, lwd = 2)
heatmap.2(as.matrix(DE[1:100, .SD, .SDcols = c(condition_ind)]), scale = "row",
          Rowv = myDendrogram, Colv = F, 
          dendrogram = "row", 
          density.info = "none", 
          trace = "none", cexRow = 0.7, cexCol = heatmap_col_wordsize,
          labRow = DE[[Genes_col]], srtCol = +45, offsetRow = 5e-10,
          lhei = c(1,8), lwid = c(4,8), margins = c(10,25))
dev.off()

# # BIOMART ####
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 
# ensemblFilters <- listFilters(ensembl)
# ensemblFilters[1:5,]
# 
# ensemblAttributes <- listAttributes(ensembl)
# ensemblAttributes[1:5,]
# 
# tabDE$Ensembl_GeneStable_ID <- getBM(attributes = "ensembl_gene_id",
#                                      filters = 'hgnc_symbol',
#                                      values = tabDE$Gene.names,
#                                      mart = ensembl)

## STRING ####
# Get the taxonomy list, homo sapiens is "9606"
taxList <- get_STRING_species(version="10", species_name=NULL)

string_db <- STRINGdb$new(version="10", species = 9606, score_threshold = 0)

# List all the methods available.
STRINGdb$methods() 

# Map the gene symbols to the STRING database
geneMap <- as.data.table(string_db$map(my_data_frame = result, 
                                       my_data_frame_id_col_names = Genes_col,
                                       removeUnmappedRows = FALSE))

# Get the genes without STRING_id
unmappedGene <- geneMap[which(is.na(geneMap$STRING_id)), .SD, .SDcols = Genes_col]

# Subset geneMap without NA and order by q-value to get the top 200 and top 100 DEGs
geneMap <- geneMap[which(!(is.na(geneMap$STRING_id))),]
setorderv(geneMap, cols = qval_col)

# Get the top 100 and 200 hits and make a network ####
top200 <- geneMap$STRING_id[1:200]
top100 <- geneMap$STRING_id[1:100]

# Filter by p-adj (already did) and add a color column 
# green for down-regulated gened and red for up-regulated genes
geneMap <- string_db$add_diff_exp_color(screen = geneMap, logFcColStr = LFC_col)

# Post payload information to the STRING server
payload_id <- string_db$post_payload(geneMap$STRING_id, colors = geneMap$color)

# Display a STRING network png with the "halo"
legend.txt <- c("- Network nodes represent proteins",
                "- Edges represent protein-protein interactions",
                "- Strong upreguated to downregulated genes",
                "  change colour intensity of halo from red to green")

pdf(paste(currentStorage, "Fig", "ProteinNetwork_PPI", paste("ProteinNetwork", dataType, contrast, "top200.pdf", sep = "_"), sep = "/"), width = 12, height = 10)
string_db$plot_network(top200, payload_id = payload_id)
legend("topleft", legend = legend.txt, title = "Legend", cex = 0.7)
dev.off()

pdf(paste(currentStorage, "Fig", "ProteinNetwork_PPI", paste("ProteinNetwork", dataType, contrast, "top100.pdf", sep = "_"), sep = "/"), width = 12, height = 10)
string_db$plot_network(top100, payload_id = payload_id)
legend("topleft", legend = legend.txt, title = "Legend", cex = 0.7)
dev.off()

# # Plot the PPI Enrichment for the best 1000 genes ####
# # If the experiment was successful, 
# # and the top hits have protein-protein interactions, 
# # you should see more enrichment at the beginning of the list than at the end.
# # Use the enrichment graph to help deﬁne a threshold on the number of proteins to consider 
# # for example, if see a strong enrichment up to position 600 on your list, 
# # this means that the signal is probably sparsed to cover the best 600 genes
# png("edgeR_ProteinNetwork_Enrichment.png", width = 4000, height = 4000, res = 300)
# string_db$plot_ppi_enrichment(geneMap$STRING_id[1:1000], quiet=TRUE)
# dev.off()

## Enrichment GKPI ####

# Set the background on all DE genes => get the correct p-value
backgroundV <- geneMap$STRING_id
string_db$set_background(backgroundV)

# Get DEGs
DEG <- unlist(geneMap[which(geneMap[[qval_col]] < fdr_DEG), c("STRING_id")])

# Enrichment GO ####
GO <- setDT(string_db$get_enrichment(DEG, category = "Process", methodMT = "fdr", iea = FALSE))
GO <- GO[pvalue_fdr < fdr_enrichment]
GO$logFDR <- -log(GO$pvalue_fdr)
GO <- GO[setDT(string_db$get_term_proteins(term_ids = GO$term_id, string_ids = DEG, enableIEA = FALSE),
               key = "term_id")[ , .(Genes = paste(preferred_name, collapse = " | ")), by = "term_id"],
         nomatch = 0L, on = "term_id"]
setorder(GO, cols = pvalue_fdr)

# Enrichment KEGG ####
KEGG <- setDT(string_db$get_enrichment(DEG, category = "KEGG", methodMT = "fdr", iea = TRUE))
KEGG <- KEGG[pvalue_fdr < fdr_enrichment]
KEGG$logFDR <- -log(KEGG$pvalue_fdr)
KEGG <- KEGG[setDT(string_db$get_term_proteins(term_ids = KEGG$term_id, string_ids = DEG, enableIEA = TRUE),
                   key = "term_id")[, .(Genes = paste(preferred_name, collapse = " | ")), by = "term_id"],
             nomatch = 0L, on = "term_id"]
setorder(KEGG, cols = pvalue_fdr)

# Enrichment PFAM ####
PFAM <- setDT(string_db$get_enrichment(DEG, category = "Pfam", methodMT = "fdr"))
PFAM <- PFAM[pvalue_fdr < fdr_enrichment]
PFAM$logFDR <- -log(PFAM$pvalue_fdr)
PFAM <- PFAM[setDT(string_db$get_term_proteins(term_ids = PFAM$term_id, string_ids = DEG),
                   key = "term_id")[, .(Genes = paste(preferred_name, collapse = " | ")), by = "term_id"],
             nomatch = 0L, on = "term_id"]
setorder(PFAM, cols = pvalue_fdr)

# Enrichment InterPro ####
InterPro <- setDT(string_db$get_enrichment(DEG, category = "InterPro", methodMT = "fdr"))
InterPro <- InterPro[pvalue_fdr < fdr_enrichment]
InterPro$logFDR <- -log(InterPro$pvalue_fdr)
InterPro <- InterPro[setDT(string_db$get_term_proteins(term_ids = InterPro$term_id, string_ids = DEG),
                           key = "term_id")[, .(Genes = paste(preferred_name, collapse = " | ")), by = "term_id"],
                     nomatch = 0L, on = "term_id"]
setorder(InterPro, cols = pvalue_fdr)

# Put all enrichments in a list ####
GKPI <- list(GO, KEGG, PFAM, InterPro)
names(GKPI) <- GKPI_folders

# Print to .csv and create Word clouds ####
for (i in 1:length(GKPI_folders)) {
  
  # Save in csv
  fwrite(x = GKPI[[i]], file = paste(currentStorage, "Fig", GKPI_folders[[i]],
                                     paste(
                                       paste(GKPI_folders[[i]], dataType, contrast,
                                             paste("fdr", fdr_enrichment*100, sep = ""),
                                             paste("DEGfdr", fdr_DEG*100, sep = ""),
                                             sep = "_"),
                                       "csv",
                                       sep = "."),
                                     sep = "/"),
         row.names = FALSE)
  
  # Word cloud
  word <- wordcloud2(GKPI[[i]][, .SD, .SDcols = c("term_description", "logFDR")], minRotation = 0, maxRotation = 0, size = GKPI_wordsize[[i]])
  
  # Save in html
  saveWidget(widget = word,
             file = paste(currentStorage, "Fig", GKPI_folders[[i]],
                          paste(
                            paste(GKPI_folders[[i]], "Word", dataType, contrast,
                                  paste("fdr", fdr_enrichment*100, sep = ""),
                                  paste("DEGfdr", fdr_DEG*100, sep = ""),
                                  sep = "_"),
                            "html",
                            sep = "."),
                          sep = "/"), 
             selfcontained = F)
  
  # Save in png
  webshot::webshot(url = paste(currentStorage, "Fig", GKPI_folders[[i]],
                               paste(
                                 paste(GKPI_folders[[i]], "Word", dataType, contrast,
                                       paste("fdr", fdr_enrichment*100, sep = ""),
                                       paste("DEGfdr", fdr_DEG*100, sep = ""),
                                       sep = "_"),
                                 "html",
                                 sep = "."),
                               sep = "/"),
                   file = paste(currentStorage, "Fig", GKPI_folders[[i]],
                                paste(
                                  paste(GKPI_folders[[i]], "Word", dataType, contrast,
                                        paste("fdr", fdr_enrichment*100, sep = ""),
                                        paste("DEGfdr", fdr_DEG*100, sep = ""),
                                        sep = "_"),
                                  "png",
                                  sep = "."),
                                sep = "/"),
                   delay = 20, vwidth = GKPI_word_width[[i]], vheight = GKPI_word_height[[i]])
}

## Get Specific Protein Information ####

# Main protein/gene
tp53 <- string_db$mp("tp53")

# Get neighbors
tp53_neighbors <- data.table(STRING_id = string_db$get_neighbors(tp53), key = "STRING_id")

# Join with DEG
tp53_neighbors_join_DEG <- setkeyv(geneMap[tp53_neighbors, on = "STRING_id", nomatch = 0L], cols = "STRING_id")
tp53_neighbors_join_DEG <- tp53_neighbors_join_DEG[tp53_neighbors_join_DEG[[qval_col]] < fdr_DEG]

# Get interaction score
tp53_interactions_score <- list()
for (i in 1:length(tp53_neighbors_join_DEG$STRING_id)) {
  tp53_interactions_score[[i]] <- setDT(string_db$get_interactions(c(tp53, tp53_neighbors_join_DEG$STRING_id[i])), key = "to")
  tp53_interactions_score[[i]] <- within(tp53_interactions_score[[i]], {
    ifelse(to == tp53, to <- from, to <- to)
    from <- tp53
  })
}
tp53_interactions_score <- setDT(rbindlist(tp53_interactions_score), key = "to")

# Join with q-value
tp53_interactions_score <- tp53_interactions_score[tp53_neighbors_join_DEG, nomatch = 0L]

fwrite(setorderv(tp53_interactions_score[, .SD, .SDcols = c(Genes_col, "combined_score", LFC_col, pval_col, qval_col)],
                 cols = "combined_score", order = -1),
       file = paste(currentStorage,
                    paste(
                      paste("tp53_interaction_score", dataType, contrast,
                            paste("DEGfdr", fdr_DEG*100, sep = ""),
                            sep = "_"),
                      "csv",
                      sep = "."),
                    sep = "/"),
       row.names = FALSE)

# string_db$get_pubmed_interaction(tp53,itga3)

## # get the reciprocal best hits of the following protein in all the STRING species
## string_db$get_homologs_besthits(tp53, symbets = TRUE)

## # get the homologs of the following two proteins in the mouse (i.e. species_id=10090)
## string_db$get_homologs(c(tp53, atm), target_species_id=10090, bitscore_threshold=60 )

# LOAD if already have the GKPI files ####
GKPI_files <- list()
for (i in 1:length(GKPI_folders)) {
  GKPI_files[[i]] <- paste(paste(GKPI_folders[[i]], dataType, contrast, "fdr5_DEGfdr5", sep = "_"), "csv", sep = ".")
}

GKPI <- list()
for (i in 1:length(GKPI_folders)) {
  GKPI[[i]] <- fread(paste(currentStorage, "Fig", GKPI_folders[[i]], GKPI_files[[i]], sep = "/"))
}
names(GKPI) <- GKPI_files
