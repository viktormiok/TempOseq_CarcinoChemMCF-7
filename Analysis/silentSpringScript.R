########################################################################################################################################################################
#
# Filename: silentSpringScript.R
# Description: Data quality assesment and visualisation
# Last modified: 20-10-2017
#
########################################################################################################################################################################


#remove all the objects
rm(list = ls())

# maybe not
# options(stringsAsFactors = FALSE)

### Set working directory
setwd("/Users/viktorian/Documents/SilentSpring/SilentSpring Data")

### Loading packages
library(plyr)
library(DESeq2)
library(ggplot2)

###import condition table
colData <- read.csv("Condition_table_detailed.csv", header=TRUE)
row.names(colData) <- colData$well_ID

###reorder
colData <- colData[order(rownames(colData)),] 

names(colData) <- c("well_ID", "lane", "chemical", "concentration", 
                    "estrogen", "biological_replicate", "plate_assignment",
                    "technical_replicate", "plate")

# turn additional variables into factors, only keeping concentration as numeric
# for plotting purposes for now
factor_cols <- c("lane", "estrogen", "technical_replicate")
colData[factor_cols] <- lapply(colData[factor_cols], as.factor)

colData$chemical <- relevel(colData$chemical, ref = "control")

# if 0 estrogen is starved then
levels(colData$estrogen) <- c("starved", "stimulated")

#Create a "mega-variable" that combines chemical, concentration (as a factor) and estrogen
colData$mega_var <- factor(paste0(colData$chemical, "_", colData$concentration, 
                                  "uM_Estrogen_", colData$estrogen))

###import biospyder files
###replicate 1
PL1A_raw <- read.csv("PL1A.csv", header=TRUE)
PL1B_raw <- read.csv("PL1B.csv", header=TRUE)
PL1C_raw <- read.csv("PL1C.csv", header=TRUE)
PL1D_raw <- read.csv("PL1D.csv", header=TRUE)
###replicate 2
PL2A_raw <- read.csv("PL2A.csv", header=TRUE)
PL2B_raw <- read.csv("PL2B.csv", header=TRUE)
PL2C_raw <- read.csv("PL2C.csv", header=TRUE)
PL2D_raw <- read.csv("PL2D.csv", header=TRUE)
###replicate 3
PL3A_raw <- read.csv("PL3A.csv", header=TRUE)
PL3B_raw <- read.csv("PL3B.csv", header=TRUE)
PL3C_raw <- read.csv("PL3C.csv", header=TRUE)
PL3D_raw <- read.csv("PL3D.csv", header=TRUE)

###combine plates into one table
plates_all <- join_all(list(PL1A_raw, PL1B_raw, PL1C_raw, PL1D_raw,
                            PL2A_raw, PL2B_raw, PL2C_raw, PL2D_raw,
                            PL3A_raw, PL3B_raw, PL3C_raw, PL3D_raw), 
                       by = "Gene", type = "full")

###put in alphabetical order
plates_all <- plates_all[order(plates_all$Gene), ]

### Turn count dataframe into a matrix
countData_matrix <- as.matrix(plates_all[, 2:(length(plates_all))])
rownames(countData_matrix) <- plates_all$Gene

###Biospyder included a list of genes that showed no expression in our data
no_signal_genes <- c("CTSD_1642", "CTTN_1650", "CYP1A2_1699",
                     "CYP2A6_1709", "CYP2E1_1722", "DHX16_1862",
                     "E2F1_2008", "EDN1_2029", "ESRRA_2238")

# sanity check - find all zero signal genes
noexp <- countData_matrix[no_signal_genes, ]
rowMeans(noexp)
ind <- which(rowMeans(countData_matrix) == 0)
zeros <- countData_matrix[ind, ]
# there's 14 zero signal genes
dim(zeros)

dds <- DESeqDataSetFromMatrix(countData = countData_matrix,
                              colData = colData,
                              design= ~ estrogen + chemical)

### collapsing technical replicates- FYI- "rep" refers to biological replicate
dds$id <- factor(paste0(dds$mega_var, "_Plate_", dds$plate))
table(dds$id)
colData(dds)[dds$id == "Tam_10uM_Estrogen_stimulated_Plate_PL2_A",]
dds <- collapseReplicates(dds, dds$id, renameCols=TRUE)
colData(dds)[dds$id == "Tam_10uM_Estrogen_stimulated_Plate_PL2_A",]
# don't do this!
# dds$well_ID <- NULL
# dds$Technical_rep <- NULL

# This would filter out 14 zero count genes
dds <- dds[rowSums(counts(dds)) >  1, ]

# Estimate size factor of the data
dds <- estimateSizeFactors(dds)

#  Estimate despersion of the data
dds <- estimateDispersions(dds)

# Plot the despersion
plotDispEsts(dds)

###############################################################################
#   Data quality assesment and visualisation
###############################################################################

# Transform the data using variance stabilization trabsformation
vst.data <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local")

### Visualize the sample distances (euclidian distance) in the heatmaps using pheatmap package
sampleDistMatrix <- as.matrix(dist(t(assay(vst.data))))
rownames(sampleDistMatrix) <- dds$estrogen
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = dist(t(assay(vst.data))),
         clustering_distance_cols = dist(t(assay(vst.data))),
         col = colors)

# Principal component analysis of the data
pc <- prcomp(t(assay(vst.data)))
pc <- as.data.frame(pc$x)
pc$lane <- colData(dds)$lane
pc$chemical <- colData(dds)$chemical
pc$concentration <- colData(dds)$concentration
pc$estrogen <- colData(dds)$estrogen
pc$biological_replicate <- colData(dds)$biological_replicate 
pc$plate_assignment <- colData(dds)$plate_assignment

# Plot the principal component analysis
p1 <- ggplot(pc, aes(PC1, PC2, shape = estrogen, 
                     size = concentration, color = chemical)) +
  geom_point() +
  theme(legend.position = "right", legend.text = element_text(size = 8))
p1

###############################################################################
#  The PCA and sample distance plots indicate that we can split the data,     #
#  based on estrogen factor. Further statistical analysis is performed        #
#  separately for estorgen stimulated and starved data                        #
###############################################################################

# Split data set on E2 starved and stimulated

ddsEstarv <- dds[,dds$estrogen=="starved"]
ddsEstimul <- dds[,dds$estrogen=="stimulated"]


# Transform the data using variance stabilization trabsformation
vst.data0 <- varianceStabilizingTransformation(ddsEstarv, blind = TRUE, fitType = "local")

### Visualize the sample distances (euclidian distance) in the heatmaps using pheatmap package
sampleDistMatrix <- as.matrix(dist(t(assay(vst.data0 ))))
rownames(sampleDistMatrix) <- paste(ddsEstarv$chemical, ddsEstarv$concentration, sep="")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = dist(t(assay(vst.data0))),
         clustering_distance_cols = dist(t(assay(vst.data0))),
         col = colors)

select = order(apply(counts(ddsEstarv), 1, sd), decreasing=T)[1:250]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(250)
#
colnames(vst.data0) <- paste(ddsEstarv$chemical, ddsEstarv$concentration, ddsEstarv$plate, sep="_")
heatmap.2(assay(vst.data0)[select,], col = hmcol, trace="none", margin=c(10, 6), dendrogram = "column")

# Principal component analysis of the data
pc <- prcomp(t(assay(vst.data0)))
pc <- as.data.frame(pc$x)
pc$lane <- colData(ddsEstarv)$lane
pc$chemical <- colData(ddsEstarv)$chemical
pc$concentration <- colData(ddsEstarv)$concentration
pc$estrogen <- colData(ddsEstarv)$estrogen
pc$biological_replicate <- colData(ddsEstarv)$biological_replicate 
pc$plate_assignment <- colData(ddsEstarv)$plate_assignment

# Plot the principal component analysis
p1 <- ggplot(pc, aes(PC1, PC2, color = chemical, 
                     size = concentration)) +
  geom_point() +
  theme(legend.position = "right", legend.text = element_text(size = 8))
p1


