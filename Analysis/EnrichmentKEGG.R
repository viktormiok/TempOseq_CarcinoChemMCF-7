########################################################################################################################################################################
#
# Filename: silentSpringScript.R
# Description: Data enrichement analysis for KEGG and GO pathways
# Last modified: 26-10-2017
#
########################################################################################################################################################################


#remove all the objects
rm(list = ls())

### Set working directory
setwd("/Users/viktorian/Documents/SilentSpring/SilentSpring Data")

### Loading packages
library(plyr)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(UpSetR)
library(EnrichmentBrowser)
library(biomaRt)

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

colData$chem_dose <- factor(paste0(colData$chemical, "_", colData$concentration, "uM"))
colData$chem_dose <- relevel(colData$chem_dose, ref = "control_0uM")

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

no_signal_genes <- c("CTSD_1642", "CTTN_1650", "CYP1A2_1699",
                     "CYP2A6_1709", "CYP2E1_1722", "DHX16_1862",
                     "E2F1_2008", "EDN1_2029", "ESRRA_2238")

ind <- which(rowMeans(countData_matrix) == 0)
zeros <- countData_matrix[ind, ]
dim(zeros)

rownames(zeros)


dds <- DESeqDataSetFromMatrix(countData = countData_matrix,
                              colData = colData,
                              design= ~ plate + estrogen + chem_dose + chem_dose:estrogen)


dds$id <- factor(paste0(dds$chemical, "_", dds$concentration, "uM_Estrogen_", dds$estrogen, "_Plate_", dds$plate))
dds <- collapseReplicates(dds, dds$id, renameCols=TRUE)


ordered <- colData(dds)[order(colMeans(counts(dds)), decreasing = TRUE), ]
ordered$meanCnt <- sort(colMeans(counts(dds)), decreasing = TRUE)
ordered[1:20, c(1,2,3,4,5,11)]

dds <- dds[rowSums(counts(dds)) >  0, ]

dds <- estimateSizeFactors(dds)

dds <- estimateDispersions(dds, fitType = "local")
plotDispEsts(dds)


# this function will need to be changed when the model changes
getResults <- function(ch_dose, est) {
  contrast <- integer(77)
  if (ch_dose == "control") {
    contrast[13] <- 1
  } else {
    ind <- grep(ch_dose, resultsNames(dds))
    contrast[ind[1]] <- 1
    if (est) contrast[ind[2]] <- 1
  }
  # assign(paste0(ch_dose, "_Estrogen_", est), results(dds, contrast = contrast))
  return(results(dds, contrast = contrast))
}

# added a single line to this function from the UpSetR package to keep 
# the gene names in the output matrix
fromList <- function(input){
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x){x <- as.vector(match(elements, x))}))
  data[is.na(data)] <- as.integer(0); data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) !=0), ]
  names(data) <- names(input)
  rownames(data) <- elements
  return(data)
}

findDEG <- function(chem = NULL, dose = NULL, estrogen = NULL, alpha = 0.05) {
  if (is.null(chem) & is.null(dose)) {
    stop("You need to pass either chemical or dose!")
  }
  if (chem!="control" & !is.null(estrogen)) { # not if youre comparing controls - change!
    if (!(estrogen %in% c(0,1))) stop("Estrogen parameter needs to be either 1 or 0!")
  }
  if (chem == "control" & !is.null(dose)) {
    stop("You can only compare controls in the presence and absence of estrogen!")
  }
  
  if (!is.null(chem) & is.null(dose) & chem != "control") {
    ch_dose <- paste0(chem, c("_0.001uM", "_0.1uM", "_1uM", "_10uM"))
    results <- lapply(ch_dose, getResults, est = estrogen)
    sigG <- lapply(results, function(x) rownames(x)[x$padj < alpha & !is.na(x$padj)])
    names(sigG) <- paste0(ch_dose, "_Vs_Control_at_E", estrogen)
  } else if (chem == "control") {
    results <- getResults("control") 
    aux <- rownames(results)[results$padj < alpha & !is.na(results$padj)]
    sigG <- list("Estrogen_stimulated_Vs_Estrogen_starved" = aux)
  }
  # system("say Your differentialy expressed genes list is ready!")
  return(sigG)
}



dds <- nbinomWaldTest(dds, maxit = 500)

resultsNames(dds)

# Estrogen stimulated Vs starved - controls only test

Estrogen_1vs0 <- findDEG(chem = "control")
length(Estrogen_1vs0[[1]])

human=useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
ann <- getBM(attribute=c("external_gene_name", "entrezgene"), mart=human)
colnames(ann) <- c("GeneName","gene_id")


genenames <- rownames(countData_matrix)
genenames <- strsplit(genenames, "_")
symbols <- integer(501)
names(symbols) <- unique(sapply(genenames, function(x) x[1]))


estDE <- sapply(strsplit(Estrogen_1vs0[[1]], "_"), "[", 1)

for (gene in estDE) {
  if (gene %in% names(symbols)) {
    ind <- which(names(symbols) == gene)
    symbols[ind] <- 1
  }
}

SSgen <- cbind(names(symbols),symbols)
colnames(SSgen) <- c("GeneName","diffExp")

dat <- merge(ann,SSgen,by="GeneName")
dat[,3] <- as.numeric(as.character(dat[,3]))
sigGen <- dat[which(dat[,3]==1),]
noSigGen <- dat[which(dat[,3]==0),]

###############################################################################
#  Function 
###############################################################################
pathwayTable <- function(db, allG, sigG, noSigG){
  fin <- data.frame(matrix(ncol = 6, nrow = length(db)))
  for(i in 1:length(db)){
    fin[i,1] <- sub("_.*", "", names(db[i]))
    fin[i,2] <- gsub("_", " ", sub("^[^_]*", "", names(db[i])))
    fin[i,3] <- length(db[[i]])
    fin[i,4] <- length(intersect(as.numeric(db[[i]]),allG))
    fin[i,5] <- length(intersect(as.numeric(db[[i]]),sigG))
    CrossTab <- matrix(c(length(intersect(sigG,as.numeric(db[[i]]))), 
                       length(intersect(noSigG,as.numeric(db[[i]]))),
                       length(setdiff(sigG, as.numeric(db[[i]]))),
                       length(setdiff(noSigG, as.numeric(db[[i]]))), nrow = 2,)
    fin[i,6] <- fisher.test(CrossTab, alternative = "greater")$p.value
    }
  fin$X7 <- p.adjust(fin[,6], method ="fdr")
  colnames(fin) <- c("keggID", "PathName", "PathSize", "OurGen",
                   "NrSigGen", "Pval","adjPval")

  return(fin[order(fin[,7], fin[,6], decreasing = F),])
}

# KEGG pathways
kegg.gs <- get.kegg.genesets("hsa")
res <- pathwayTable(kegg.gs, dat[,2], sigGen[,2], noSigGen[,2])
datafin <- write.table(res, file="KEGG.txt", sep="\t", row.names=FALSE)

# Gene Ontology pathways
go.gs <- get.go.genesets(org="hsa", onto="BP", mode="GO.db")
res <- pathwayTable(go.gs, dat[,2], sigGen[,2], noSigGen[,2])
datafin <- write.table(res, file="KEGG.txt", sep="\t", row.names=FALSE)

###############################################################################
#   KEGG
###############################################################################
kegg.gs <- get.kegg.genesets("hsa")

fin <- data.frame(matrix(ncol = 6, nrow = length(kegg.gs)))
for(i in 1:length(kegg.gs)){
  fin[i,1] <- sub("_.*", "", names(kegg.gs[i]))
  fin[i,2] <- gsub("_", " ", sub("^[^_]*", "", names(kegg.gs[i])))
  fin[i,3] <- length(kegg.gs[[i]])
  fin[i,4] <- length(intersect(as.numeric(kegg.gs[[i]]),dat[,2]))
  fin[i,5] <- length(intersect(as.numeric(kegg.gs[[i]]),sigGen[,2]))
  fin[i,6] <- fisher.test (matrix(c(length(intersect(sigGen[,2],as.numeric(kegg.gs[[i]]))), 
                                    length(intersect(noSigGen[,2],as.numeric(kegg.gs[[i]]))),
                                    length(setdiff(sigGen[,2],as.numeric(kegg.gs[[i]]))),
                                    length(setdiff(noSigGen[,2],as.numeric(kegg.gs[[i]])))), nrow = 2,
                                    dimnames = list(c("DE", "Not.DE"), c("In.pathway", "Not.in.pathway"))),
                           alternative = "greater")$p.value
}
fin$X7 <- p.adjust(fin[,6], method ="fdr")
colnames(fin) <- c("keggID", "Pathway name", "Pathway size", "Our genes in Pathway",
                   "Nr. Significat genes", "PvalChi2","adjPvalChi2")

res <-fin[order(fin[,6], decreasing = F),]

datafin <- write.table(res, file="EnrichmentAnalysisKEGG.txt", sep="\t", row.names=FALSE, col.names=TRUE)

###############################################################################
#   GO
###############################################################################

go.gs <- get.go.genesets(org="hsa", onto="BP", mode="GO.db")


fin <- data.frame(matrix(ncol = 6, nrow = length(go.gs)))
for(i in 1:length(go.gs)){
  fin[i,1] <- sub("_.*", "", names(go.gs[i]))
  fin[i,2] <- gsub("_", " ", sub("^[^_]*", "", names(go.gs[i])))
  fin[i,3] <- length(go.gs[[i]])
  fin[i,4] <- length(intersect(as.numeric(go.gs[[i]]),dat[,2]))
  fin[i,5] <- length(intersect(as.numeric(go.gs[[i]]),sigGen[,2]))
  fin[i,6] <- fisher.test(matrix(c(length(intersect(sigGen[,2],as.numeric(go.gs[[i]]))), 
                                    length(intersect(noSigGen[,2],as.numeric(go.gs[[i]]))),
                                    length(setdiff(sigGen[,2],as.numeric(go.gs[[i]]))),
                                    length(setdiff(noSigGen[,2],as.numeric(go.gs[[i]])))), nrow = 2,
                                  dimnames = list(c("DE", "Not.DE"), c("In.pathway", "Not.in.pathway"))),
                           alternative = "greater")$p.value
}
fin$X7 <- p.adjust(fin[,6], method ="BH")
colnames(fin) <- c("keggID", "Pathway name", "Pathway size", "Our genes in Pathway",
                   "Nr. Significat genes", "PvalChi2","adjPvalChi2")

res <-fin[order(fin[,6], decreasing = F),]

datafin <- write.table(res, file="EnrichmentAnalysisKEGG.txt", sep="\t", row.names=FALSE, col.names=TRUE)

###############################################################################
#   Pathway analysis
###############################################################################



datax <- counts(dds)
rownames(datax) <- sapply(strsplit(rownames(datax), "_"), function(x) x[1])

datax <- aggregate(datax, list(rownames(datax)), median)
rownames(datax) <- datax[,1]
datax <- round(datax[,-1],0)

kegg.gs <- get.kegg.genesets("hsa")
PathGen <- dat[dat[,2]%in%as.numeric(kegg.gs$hsa05224),1]
brestDat <- datax[rownames(datax)%in%PathGen,]

library(huge)
library(rags2ridges)

X.npn = huge.npn(t(brestDat))
Cx <- covML(X.npn)
ridgeEst <- ridgeP(Cx, lambda = 1, type = "Alt")


## Obtain regularized precision under optimal penalty
#OPT <- optPenalty.LOOCV(X.npn, lambdaMin = .5, lambdaMax = 30, step = 100)

## Determine support regularized (standardized) precision under optimal penalty
ridgeEst1 <- sparsify(ridgeEst, threshold = "localFDR", FDRcut = 0.9)
p <- abs(ridgeEst1[[2]])
diag(p) <- 0
edgeHeat(p)
Ugraph(p, type = "fancy", cut = 0.07)

nodCol <- ifelse(PathGenSig == 1, "red", "green")
Ugraph(p, type = "fancy", lay = "layout_with_fr",
       cut = 0.07, Vcolor =nodCol, prune = TRUE)

# remove all the objects
rm(list = ls())

# Set working directory
setwd("/Users/viktorian/Documents/SilentSpring/SilentSpring Data")

library(plyr)
library(DESeq2)
library(EnrichmentBrowser)
library(biomaRt)

human=useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
ann <- getBM(attribute=c("external_gene_name", "entrezgene"), mart=human)
colnames(ann) <- c("GeneName","gene_id")


genenames <- rownames(countData_matrix)
genenames <- strsplit(genenames, "_")
symbols <- integer(501)
names(symbols) <- unique(sapply(genenames, function(x) x[1]))


estDE <- sapply(strsplit(Estrogen_1vs0[[1]], "_"), "[", 1)

for (gene in estDE) {
  if (gene %in% names(symbols)) {
    ind <- which(names(symbols) == gene)
    symbols[ind] <- 1
  }
}

SSgen <- cbind(names(symbols),symbols)
colnames(SSgen) <- c("GeneName","diffExp")

dat <- merge(ann,SSgen,by="GeneName")
dat[,3] <- as.numeric(as.character(dat[,3]))
sigGen <- dat[which(dat[,3]==1),]
noSigGen <- dat[which(dat[,3]==0),]

###############################################################################
#  Function 
###############################################################################
pathwayTable <- function(db, allG, sigG, noSigG){
  fin <- data.frame(matrix(ncol = 6, nrow = length(db)))
  for(i in 1:length(db)){
    fin[i,1] <- sub("_.*", "", names(db[i]))
    fin[i,2] <- gsub("_", " ", sub("^[^_]*", "", names(db[i])))
    fin[i,3] <- length(db[[i]])
    fin[i,4] <- length(intersect(as.numeric(db[[i]]),allG))
    fin[i,5] <- length(intersect(as.numeric(db[[i]]),sigG))
    CrossTab <- matrix(c(length(intersect(sigG,as.numeric(db[[i]]))), 
                         length(intersect(noSigGen[,2],as.numeric(db[[i]]))),
                         length(setdiff(sigG,as.numeric(db[[i]]))),
                         length(setdiff(sigG,as.numeric(db[[i]])))), nrow = 2,)
    fin[i,6] <- fisher.test(CrossTab, alternative = "greater")$p.value
  }
  fin$X7 <- p.adjust(fin[,6], method ="fdr")
  colnames(fin) <- c("keggID", "PathName", "PathSize", "OurGen",
                     "NrSigGen", "Pval","adjPval")
  
  return(fin[order(fin[,7], fin[,6], decreasing = F),])
}

# KEGG pathways
kegg.gs <- get.kegg.genesets("hsa")
res <- pathwayTable(kegg.gs, dat[,2], sigGen[,2], noSigGen[,2])
datafin <- write.table(res, file="KEGG.txt", sep="\t", row.names=FALSE)

# Gene Ontology pathways
go.gs <- get.go.genesets(org="hsa", onto="BP", mode="GO.db")
res <- pathwayTable(go.gs, dat[,2], sigGen[,2], noSigGen[,2])
datafin <- write.table(res, file="KEGG.txt", sep="\t", row.names=FALSE)

