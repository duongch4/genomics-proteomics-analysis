# source("http://bioconductor.org/biocLite.R")
# #### INSTALLING ####
# #install.packages("Matrix","lme4","MASS","zoo","plyr","devtools")
# install.packages("devtools", repos = "http://cran.us.r-project.org")
# library(devtools)
# library(MSqRobData)
# #### Download the package from github directly ####
# devtools::install_github("statOmics/MSqRob")
# Download the example data sets
# devtools::install_github("statOmics/MSqRobData")
# #### Download the MSqRob folder from github and manually build the package #### 
# # Create a .zip binary file that can be installed:
# build(pkg = paste(source_dir, "MSqRob-master", sep = "/"), path = NULL, binary = TRUE, vignettes = FALSE, args = NULL, quiet = FALSE)
# # Install the .zip file:
# install.packages(paste(source_dir, "MSqRob_0.2.zip", sep = "/"), repos=NULL)

#### Load libraries and Source other files if needed ####
sourceFunctions <- "C:/Users/chiba/OneDrive/R_myFunctions"
source(paste(sourceFunctions, "Graphs.R", sep = "/"))

# Load libraries
library(MSqRob)
library(MSnbase)
library(Biobase)
library(data.table)
library(splitstackshape)

### 1. Load peptides.txt dataset ####
setwd("C:/Users/chiba/OneDrive/Courses/Summer_2017/Proteomics")
file_peptides_txt <- paste(getwd(), "data", "peptides.txt", sep = "/")
options(datatable.integer64 = "numeric")
peptides <- fread(file = file_peptides_txt, sep = "\t", header = TRUE, quote = "", check.names = TRUE)

colnames(peptides) <- sub(pattern = "G581", replacement = "G583", fixed = TRUE, x = colnames(peptides))
colnames(peptides) <- sub(pattern = "G564_LLY283_R1_CID", replacement = "G564_LLY283_R4", fixed = TRUE, x = colnames(peptides))
colnames(peptides) <- sub(pattern = "(_HCD)|(_CID)", replacement = "", perl = TRUE, x = colnames(peptides))
colnames(peptides)

### 2. Pre-processing with proteinGroups.txt file ####
file_proteinGroups <- paste(getwd(), "data", "proteinGroups.txt", sep = "/")
proteinGroups <- fread(file = file_proteinGroups, sep = "\t", header = TRUE, quote = "", check.names = TRUE)

# Extract the column that indicates which proteins are only identified by a modification site
only_site <- proteinGroups[["Only.identified.by.site"]]

# If there are no such proteins (as is the case here), the column will be completely empty and R will import this by default as "NA"
# However, we want an empty value instead of "NA" to keep the col consistent with the "Contaminant|Potential.contaminant" and "Reverse" cols
only_site[is.na(only_site)] <- ""
# Select the protein accessions that are only identified by one or more modified peptides
filter_symbol <- "+"
removed_proteins <- proteinGroups[["Protein.IDs"]][only_site == filter_symbol]

# Create a logical variable "removed" that holds "TRUE" if a row in the peptide data frame should be removed and FALSE otherwise.
removed <- unlist(peptides[["Proteins"]]) %in% removed_proteins

# Create a new column "only_site" in "peptide" that indicates with a "+" which rows should be removed.
peptides$only_site <- ""
peptides$only_site[removed] <- "+"

## Create an experimental annotation data frame ####
run <- colnames(peptides)[grepl(pattern = "Intensity.", x = colnames(peptides))]
run <- sub(pattern = "Intensity.", replacement = "", x = run)
run

# Cell lines
cell <- factor(sapply(strsplit(run, "_"), "[[", 1))
cell

# Treatment
treatment <- factor(sapply(strsplit(run, "_"), "[[", 2))
treatment

# cell_treatment as one factor (IMPORTANT!!!)
cell_treatment <- factor(paste(cell, treatment, sep = "_"))
cell_treatment

# # "Biological" replicates (replicates per cell_treatment)
# biorep <- list()
# for (i in 1:length(levels(cell_treatment))) {
#   biorep[[i]] <- rep(paste("b_", i, sep =""), length(grep(pattern = levels(cell_treatment)[i], x = cell_treatment)))
# }
# biorep <- factor(unlist(biorep))
# biorep

# Combine
exp_annotation <- data.frame(run = run, cell_treatment = cell_treatment)
exp_annotation

system.time(peptides2 <- preprocess_wide(peptides, accession = "Proteins", split = ";",
                                         exp_annotation = exp_annotation,
                                         quant_cols = "Intensity.",
                                         aggr_by = "Sequence",
                                         aggr_function = "sum",
                                         logtransform = TRUE, base = 2,
                                         normalisation = "quantiles",
                                         smallestUniqueGroups = TRUE,
                                         useful_properties = c("Gene.names", "Protein.names", "Sequence"),
                                         filter = colnames(peptides)[which(grepl(pattern = "(Reverse)|([Cc]ontaminant)|(only_site)",
                                                                                 x = colnames(peptides), perl = TRUE))],
                                         filter_symbol="+", minIdentified=2)
)

attr(peptides2,"MSqRob_exp_annotation")

# Select columns with "Intensities." and assign them as indices to quant_cols
quant_cols <- which(grepl("Intensity.", colnames(peptides2)))
quant_cols

# Change the col names of peptides2 to match with the exp_annotation
colnames(peptides2) <- sub(pattern = "Intensity.", replacement = "", x = colnames(peptides2))
colnames(peptides2)

saveRDS(peptides2, file = "peptides_norm.rds")

#Call df2protdata function
system.time(proteins <- df2protdata(peptides2, acc_col = "Proteins",
                                    quant_cols = quant_cols, quant_name = "quant_value",
                                    run_name = "run", annotations = c("Gene.names", "Proteins", "Protein.names")))
# user  system elapsed 
# 28.44    0.00   29.06   
#save(proteins, file = "proteins.RData")

# Extract the quant_value from protdata object and turn the df's inside into data.tables ####
proteinsQuant <- attr(proteins, "data")

for (i in 1:length(proteinsQuant)) {
  # Add Proteins accordingly
  proteinsQuant[[i]] <- proteinsQuant[[i]][ , .(meanQuant = mean(quant_value),
                                                Proteins = attr(proteins, "accession")[[i]],
                                                Gene.names = attr(proteins, "annotation")[i]),
                                            by = "run"]
  # Mean of quant_value by each of 28 runs
  proteinsQuant[[i]] <- dcast.data.table(data = proteinsQuant[[i]], formula = Proteins+Gene.names ~ run, value.var = "meanQuant")
}

proteinsQuant <- rbindlist(proteinsQuant, use.names = TRUE, fill = TRUE)

# SAVE!!! ####
saveRDS(proteinsQuant, file = "proteinsQuant.rds")

# # If want to add more columns to the protdata object => "addVarFromVar()"
# techrep <- as.factor(paste0("t_",1:18))
# names(techrep) <- levels(exp_annotation$run)
# techrep
# proteins <- addVarFromVar(proteins, basecol="run", name="techrep", vector=techrep)
# proteins

### 4. Fit Robust Ridge Models with Huber Weights for each protein (TAKES TIME!) ####
system.time(protLM <- fit.model(proteins, response = "quant_value",
                                fixed = c("cell_treatment"), random = c("run", "Sequence"),
                                add.intercept = TRUE, weights = "Huber"))
# user  system elapsed 
# 1692.91    4.66 1703.11
# SAVE!!! ####
saveRDS(protLM, file = "RobRidgeModel.rds")

# ## Inspect the model ####
# getAccessions(protLM)[1:5]
# getModels(protLM)[1:5]
# getAnnotations(protLM)[1:5,]
# 
# # Models 24 and 50
# protLM[c(24,50)]
# 
# # Model for accession "WP_003033643"
# protLMF["WP_003033643"]
# 
# # Extract and store ONE model for close inspection
# modelWP_003040227 <- getModels(protLMFranc["WP_003040227"], simplify=TRUE)
# betaBVcovDf <- getBetaVcovDf(modelWP_003040227)
# 
# betaBVcovDf$beta
# betaBVcovDf$vcov
# betaBVcovDf$df
# betaBVcovDf$df_exp
# betaBVcovDf$sigma
# 
# # Extract and store >1 models for close inspection
# betaVcovDfList <- getBetaVcovDfList(protLMFranc[1:2])
# str(betaVcovDfList)

### 5. INFERENCES ####
attr(getModels(protLM[3]), "MSqRob_levels")

## Contrasts ####

# Treatment Only model
# L <- makeContrast(contrasts = c(""), levels)

# GSK591 vs SGC2096
L1 <- makeContrast(contrasts = c("(cell_treatmentG561_GSK591 + cell_treatmentG564_GSK591 + cell_treatmentG583_GSK591)/3 -
                                (cell_treatmentG561_SGC2096 + cell_treatmentG564_SGC2096 + cell_treatmentG583_SGC2096)/3", # Overall effect
                                 
                                 # Effect on a single cell line
                                 "cell_treatmentG561_GSK591 - cell_treatmentG561_SGC2096",
                                 "cell_treatmentG564_GSK591 - cell_treatmentG564_SGC2096",
                                 "cell_treatmentG583_GSK591 - cell_treatmentG583_SGC2096",
                                 
                                 # Difference between 2 cell lines (is there an interaction effect?)
                                 "(cell_treatmentG561_GSK591 - cell_treatmentG561_SGC2096) -
                                 (cell_treatmentG564_GSK591 - cell_treatmentG564_SGC2096)",
                                 
                                 "(cell_treatmentG561_GSK591 - cell_treatmentG561_SGC2096) -
                                 (cell_treatmentG583_GSK591 - cell_treatmentG583_SGC2096)",
                                 
                                 "(cell_treatmentG564_GSK591 - cell_treatmentG564_SGC2096) -
                                 (cell_treatmentG583_GSK591 - cell_treatmentG583_SGC2096)"),
                   
                   levels = c("cell_treatmentG561_GSK591", "cell_treatmentG561_LLY283", "cell_treatmentG561_SGC2096",
                              "cell_treatmentG564_GSK591", "cell_treatmentG564_LLY283", "cell_treatmentG564_SGC2096",
                              "cell_treatmentG583_GSK591", "cell_treatmentG583_LLY283", "cell_treatmentG583_SGC2096"))


# LLY283 vs SGC2096
L2 <- makeContrast(contrasts = c("(cell_treatmentG561_LLY283 + cell_treatmentG564_LLY283 + cell_treatmentG583_LLY283)/3 -
                                (cell_treatmentG561_SGC2096 + cell_treatmentG564_SGC2096 + cell_treatmentG583_SGC2096)/3", # Overall effect
                                 
                                 # Effect on a single cell line
                                 "cell_treatmentG561_LLY283 - cell_treatmentG561_SGC2096",
                                 "cell_treatmentG564_LLY283 - cell_treatmentG564_SGC2096",
                                 "cell_treatmentG583_LLY283 - cell_treatmentG583_SGC2096",
                                 
                                 # Difference between 2 cell lines (is there an interaction effect?)
                                 "(cell_treatmentG561_LLY283 - cell_treatmentG561_SGC2096) -
                                 (cell_treatmentG564_LLY283 - cell_treatmentG564_SGC2096)",
                                 
                                 "(cell_treatmentG561_LLY283 - cell_treatmentG561_SGC2096) -
                                 (cell_treatmentG583_LLY283 - cell_treatmentG583_SGC2096)",
                                 
                                 "(cell_treatmentG564_LLY283 - cell_treatmentG564_SGC2096) -
                                 (cell_treatmentG583_LLY283 - cell_treatmentG583_SGC2096)"),
                   
                   levels = c("cell_treatmentG561_GSK591", "cell_treatmentG561_LLY283", "cell_treatmentG561_SGC2096",
                              "cell_treatmentG564_GSK591", "cell_treatmentG564_LLY283", "cell_treatmentG564_SGC2096",
                              "cell_treatmentG583_GSK591", "cell_treatmentG583_LLY283", "cell_treatmentG583_SGC2096"))

# LLY283 vs GSK591
L3 <- makeContrast(contrasts = c("(cell_treatmentG561_LLY283 + cell_treatmentG564_LLY283 + cell_treatmentG583_LLY283)/3 -
                                 (cell_treatmentG561_GSK591 + cell_treatmentG564_GSK591 + cell_treatmentG583_GSK591)/3", # Overall effect
                                 
                                 # Effect on a single cell line
                                 "cell_treatmentG561_LLY283 - cell_treatmentG561_GSK591",
                                 "cell_treatmentG564_LLY283 - cell_treatmentG564_GSK591",
                                 "cell_treatmentG583_LLY283 - cell_treatmentG583_GSK591",
                                 
                                 # Difference between 2 cell lines (is there an interaction effect?)
                                 "(cell_treatmentG561_LLY283 - cell_treatmentG561_GSK591) -
                                 (cell_treatmentG564_LLY283 - cell_treatmentG564_GSK591)",
                                 
                                 "(cell_treatmentG561_LLY283 - cell_treatmentG561_GSK591) -
                                 (cell_treatmentG583_LLY283 - cell_treatmentG583_GSK591)",
                                 
                                 "(cell_treatmentG564_LLY283 - cell_treatmentG564_GSK591) -
                                 (cell_treatmentG583_LLY283 - cell_treatmentG583_GSK591)"),
                   
                   levels = c("cell_treatmentG561_GSK591", "cell_treatmentG561_LLY283", "cell_treatmentG561_SGC2096",
                              "cell_treatmentG564_GSK591", "cell_treatmentG564_LLY283", "cell_treatmentG564_SGC2096",
                              "cell_treatmentG583_GSK591", "cell_treatmentG583_LLY283", "cell_treatmentG583_SGC2096"))

## Results ####
# Choose contrast to output result ####
contrast <- L3
result_file <- "result_LLY283_GSK591.rds"
result_file2 <- "result_fixed_LLY283_GSK591.rds"

# Output ####
system.time(result <- test.protLMcontrast(protLM, contrast))
# user  system elapsed 
# 106.19    2.53  110.94

result <- prot.p.adjust(result, method = "fdr")
result <- prot.signif(result, level = 0.05)

result <- lapply(result, function(df) {setDT(df, key = "Proteins")})

# Join the result with proteinsQuant data
result <- lapply(result, function(dt) {dt[proteinsQuant, on = "Proteins"]})

# Are the 2 cols of Gene.names matched? If yes, remove one of them
sapply(result, function(dt) {table(dt[["Gene.names"]] == dt[["i.Gene.names"]])})

result <- lapply(result, function(dt) {dt[, c("i.Gene.names") := NULL]})

# SAVE!!! ####
saveRDS(result, file = result_file)

# Filtering ####
for (i in 1:length(result)) {
  
  # Keep rows with Gene.names
  result[[i]] <- result[[i]][Gene.names != ""]
  
  # Remove rows with NA in q-value
  result[[i]] <- result[[i]][!is.na(qval)]
  
  # Convert factors to characters (all factors cols)
  update_cols <- sapply(result[[i]], is.factor)
  result[[i]][, names(result[[i]])[update_cols] := lapply(.SD, as.character), .SDcols = update_cols]
  
  # Only keep the first gene if multiple genes detected in one row (separated by ";")
  update_genes <- result[[i]][["Gene.names"]][grep(pattern = ";", x = result[[i]][["Gene.names"]], perl = TRUE)]
  update_genes <- sapply(
    strsplit(update_genes, split = ";", fixed = TRUE),
    function(genes) {genes <- genes[1]}
  )
  result[[i]][["Gene.names"]][grep(pattern = ";", x = result[[i]][["Gene.names"]], perl = TRUE)] <- update_genes
  
  # If > 1 genes shown, select one with the smallest q-value (if need only 1 min then use "which.min(qval)")
  result[[i]] <- result[[i]][
    result[[i]][ , .(row_index = .I[qval == min(qval)]), by = "Gene.names"][["row_index"]]
  ]

}

# Checking q-values of those who still have same gene names
# If all of them > 0.05 => just keep one of each
qval_inspect <- list()
for (i in 1:length(result)) {
  qval_inspect[[i]] <- result[[i]][["qval"]][which(result[[i]][, count := .N, by = "Gene.names"][["count"]] > 1)]
}
sapply(qval_inspect, table)

# All qval > 0.05 => remove duplicates
result <- lapply(result, function(dt) {dt <- dt[!duplicated(Gene.names)]})
#sapply(result, function(dt) table(dt[["signif"]]))

# SAVE!!! ####
saveRDS(result, file = result_file2)
