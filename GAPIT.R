#----
# Dependencies
source("https://zzlab.net/GAPIT/GAPIT.library.R")
source("https://zzlab.net/GAPIT/gapit_functions.txt")


install.packages("pkgload")
install.packages("usethis")
install.packages("LDheatmap")
install.packages("EMMREML")
install.packages("scatterplot3d")
install.packages("ape")
install.packages("lmerTest")
install.packages("emmeans")
install.packages("tidyverse")
install.packages("pbkrtest")
install.packages("tidyr")
install.packages("broom")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("multtest")
BiocManager::install("snpStats")
BiocManager::install("qvalue")
??qvalue
setwd("C:\\Users\\joshu\\OneDrive\\Josh\\Summer 2021\\John Innes Centre\\GAPIT")

library(lmerTest)
library(emmeans)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(qvalue)
#----
# GWAS analysis function
BethFunction <- function(CSVPath, SNPPath, MatrixPath, ContainsIDs){
  CSV <- read.csv(CSVPath) # Phenotype (Y)
  SNP <- read.delim(SNPPath, header=T) # Genotype in hapmap format (G)
  Covariate <- read.delim(MatrixPath, header=T) # Kinship Table (KI)
  ContainsIDs <- read.csv(ContainsIDs)[,1:2]
  
  # Normality testing
  PValues <- c(as.numeric(shapiro.test(CSV$trait)[2]), as.numeric(shapiro.test(log10(CSV$trait))[2]), as.numeric(shapiro.test(sqrt(CSV$trait))[2]))
  CSV_Copy <- CSV
  
  if(PValues[1] > 0.05 & PValues[2] > 0.05 & PValues[3] > 0.05){
    Phenotype_Vals <- CSV$trait # Use raw data
  } else {
    Phenotype_Vals <- list(CSV$trait, log10(CSV$trait), sqrt(CSV$trait))[[which.min(replace((0.05-PValues), PValues<=0, NA))]] # Use p-value closest to but below 0.5
  }
  CSV_Copy$trait <- Phenotype_Vals
  
  # LMM to estimate means
  Estimate_mean_results <- lmer(formula=trait~(1|location)*genotype_name, data=CSV_Copy)
  Emmeans <- emmeans(object=Estimate_mean_results, ~genotype_name)
  Trait <- tidy(Emmeans)
  
  # Change colnames and genotype_name to ID value:
  # Without NA Values
  Trait <- merge(Trait, ContainsIDs, by = "genotype_name")
  #print(Trait$genotype_name[Trait$genotype_name %in% Trait2$genotype_name == F])
  
  # With NA Values
  #Trait$X.Trait. <- ContainsIDs$X.Trait[match(Trait$genotype_name, ContainsIDs$genotype_name, nomatch=NA)]
  
  LMMMeans <- data.frame(X.Trait.=Trait$X.Trait., Autoground=Trait$estimate)
  colnames(LMMMeans)[1] <- "<Trait>"
  write.table(LMMMeans, file="Josh's Magical LMMMeans.txt", sep="\t", quote=F, row.names=F)
  
  # GAPIT run with any PCA number, function for working out best pca total
  GAPITPCA <- GAPIT(Y=LMMMeans, G=SNP, CV=Covariate, PCA.total=3, Model.selection=T)
  
  # Value should be corresponding PC number to highest BIC value in GAPIT.MLM.Autoground.BIC...
  BICValues <- read.csv("GAPIT.MLM.Autoground.BIC.Model.Selection.Results.csv")
  PCofBestBIC <- BICValues[which.max(BICValues[,2]),1]
  
  # Use this PC number in the final model
  ModelsVector <- c("GLM", "MLM", "BLINK", "FarmCPU")
  GAPITFINAL <- GAPIT(Y=LMMMeans, G=SNP, CV=Covariate, PCA.total=PCofBestBIC, Multiple_analysis = T, model = ModelsVector, Inter.Plot=T) #Emma?
}
# Example Data
BethFunction(CSVPath="RAW-AutoGround.csv", SNPPath="SNP-350k.hmp.txt", MatrixPath="Q-matrix-q1-q4.txt", ContainsIDs="20210303_curated_means_for_GWAS-post_outliers_FINAL.csv")
# fOne
BethFunction(CSVPath="RAW-fOne.csv", SNPPath="SNP-350k.hmp.txt", MatrixPath="Q-matrix-q1-q4.txt", ContainsIDs="20210303_curated_means_for_GWAS-post_outliers_FINAL.csv")
# rTwo
BethFunction(CSVPath="RAW-rTwo.csv", SNPPath="SNP-350k.hmp.txt", MatrixPath="Q-matrix-q1-q4.txt", ContainsIDs="20210303_curated_means_for_GWAS-post_outliers_FINAL.csv")
# rThree
BethFunction(CSVPath="RAW-fThree.csv", SNPPath="SNP-350k.hmp.txt", MatrixPath="Q-matrix-q1-q4.txt", ContainsIDs="20210303_curated_means_for_GWAS-post_outliers_FINAL.csv")

#----
# Work out which model is best fit using sum of least squares
ModelSelection <- function(ModelsVector=c("GLM", "MLM", "Blink", "FarmCPU"), directory, QLambda=0.05){
  setwd(directory)
  SumSquares <- c()
  #par(mfrow=c(2,2))
  for(Model in ModelsVector){
    TempFile <- read.csv(paste("GAPIT.", Model, ".Autoground.GWAS.Results.csv", sep=""))
    LogActualPvals <- -log10(TempFile$P.value[order(TempFile$P.value)])
    LogExpectedPvals <- -log10((1:length(TempFile$P.value))/(length(TempFile$P.value)+1))
    #PvalsDF <- data.frame(LogActualPvals, LogExpectedPvals)
    #ModelVals <- lm(LogActualPvals~LogExpectedPvals)
    #SumSquaresResid <- sum(as.numeric(ModelVals$residuals)^2)
  
    # Sequence which matches the actual pvals on the x axis and the y axis
    SumSquaresResid <- sum((LogActualPvals-LogExpectedPvals)^2)
  
    #print(mean(LogActualPvals-LogExpectedPvals))
    #plot(LogExpectedPvals, LogExpectedPvals, type="l", main=Model, ylab="Log10 Actual P-Values", xlab="Log10 Expected P-Values")
    #lines(y=LogActualPvals, x=LogExpectedPvals, col="red")
    #hist(LogActualPvals, main=Model)
  
    SumSquares[match(Model, ModelsVector)] <- SumSquaresResid
  }
  #print(SumSquares)
  #LowestSumSquares <- SumSquares[which.min(SumSquares)]
  BestModel <- ModelsVector[which.min(SumSquares)]
  #print(paste("The best model is", BestModel))
  
  GWASResults <- read.csv(paste("GAPIT.", BestModel, ".Autoground.GWAS.Results.csv", sep=""))
  Qvalues <- qvalue(GWASResults$P.value) # These are pvalues adjusted for FDR
  SignifP <- subset(Qvalues$pvalues, Qvalues$qvalues < QLambda)
  Num_Signif <- length(SignifP)
  SignifGWAS <- GWASResults[1:Num_Signif,]
  return(SignifGWAS)
}  

Example_data <- ModelSelection(directory="C:\\Users\\joshu\\OneDrive\\Josh\\Summer 2021\\John Innes Centre\\GAPIT\\Example Data") # Blink
RAWF1 <- ModelSelection(directory="C:\\Users\\joshu\\OneDrive\\Josh\\Summer 2021\\John Innes Centre\\GAPIT\\RAWF1") # GLM
RAWR2 <- ModelSelection(directory="C:\\Users\\joshu\\OneDrive\\Josh\\Summer 2021\\John Innes Centre\\GAPIT\\RAWR2") # Blink
RAWF3 <- ModelSelection(directory="C:\\Users\\joshu\\OneDrive\\Josh\\Summer 2021\\John Innes Centre\\GAPIT\\RAWF3") # Blink

RAWF3[,c(1,2,3,4,9)]

# csi (sciverse) needs to be time period
# https://cyverseuk.org/applications/causal-structure-inference-csi/
# cytoscape with gml files
# look in arabidopsis for parent genes
