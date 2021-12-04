#----

# https://bioconductor.riken.jp/packages/3.0/bioc/vignettes/qvalue/inst/doc/qvalue.pdf !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# qvalue <- qvalue(pvalues)
# plot as histogram
# estimation of fdr function also - can be placed on manhattan plot also

# PLINK
# get q value and fdr from p value using package (also bonferroni threshold)
# Then put into Hugh's script which plots the Manhattan plot
# ggplot plotter for decay

install.packages("plink")
library(plink)

#----
# LD Decay Plotter in R - LD.decay
install.packages("sommer")
library(sommer)

# LD.decay(markers, map, silent=F, unlinked=F, gamma=0.95)

# markers = numeric matrix of markers (cols) by individuals (rows) in -1, 0, 1, format
# map = data frame with columns 1) name of marker, 2) linkage group or chrom, 3) Position in cM or bp

#----
# LD Link R https://www.frontiersin.org/articles/10.3389/fgene.2020.00157/full
# We can't use the LDlinkR package, because that is only for humans - one of the arguments is a population from 1000 Genomes Project
# You also have to apply for a user token to use online
install.packages("LDlinkR")
library(LDlinkR)
?LDproxy
?LDhap
?LDpair
?LDpop

#----
# TASSEL example https://avikarn.com/2019-07-22-GWAS/
# Need to make sure that the taxa column is the same in both files.
setwd("C:\\Users\\joshu\\OneDrive\\Josh\\Summer 2021\\John Innes Centre\\TASSEL")

# TASSEL phenotype and genotype input data need to have the same taxa values
phenotype <- read.csv("RAW-Autoground.csv")
names <- read.csv("names.csv")[,1:2]
Trait <- merge(phenotype, names, by = "genotype_name")
Trait[,1] <- Trait$X.Trait.
Trait <- Trait[,c(1, 3, 6)]
write.csv(Trait, file="TASSELPhenotype.csv", quote=F, row.names=F)

# TASSEL cannot have duplicated data for GLM, need to take a mean for each taxa - this didn't work, don't know why
Trait <- aggregate(Trait[, 3], list(Trait$genotype_name, Trait$location), mean)
colnames(Trait) <- c("Taxa", "Location", "Trait")
write.csv(Trait, file="TASSELPhenotype.csv", quote=F, row.names=F)

# MLM works though.
mlm_stats <- read.table("mlm_stats.txt", header = T, sep = "\t")
head(mlm_stats)

# GWAS Significance Threshold
GWAS_Signif_Thresh <- -log10(0.05/length(mlm_stats$Marker))

# Adjusting P-values for multiple comparisons with bonferroni and FDR
library(dplyr)

adjusted_mlm <- mlm_stats %>%
  transmute(Marker, Chr, Pos, p,
            p_Bonferroni = p.adjust(mlm_stats$p, "bonferroni"),
            p_FDR = p.adjust(mlm_stats$p, "fdr")
            )
head(adjusted_mlm)

write.csv(adjusted_mlm, file="adjusted_p_MLM.csv", quote = T, eol = "\n", na = "NA")

# QQ plots
install.packages("qqman")
library(qqman)

par(mfrow=c(1,3))
qqman::qq(adjusted_mlm$p, main = "non-adjusted P-value")
qqman::qq(adjusted_mlm$p_Bonferroni, main = "Bonferroni")
qqman::qq(adjusted_mlm$p_FDR, main = "FDR")

mlm_stats
mlm_test <- subset(mlm_stats, !is.na(mlm_stats$p) == T)

par(mfrow=c(1,1))
qqman::manhattan(mlm_test, chr="Chr", bp="Pos", snp="Marker", p="p", annotateTop = T, suggestiveline = F, genomewideline = GWAS_Signif_Thresh, ylim=c(0, 7))
# None greater than the significance threshold which is why there is no line

# Estimating and plotting LD decay over distance
# May need to split by chromosome in order to do 
library(ggplot2)

#----
setwd("C:\\Users\\joshu\\OneDrive\\Josh\\Summer 2021\\John Innes Centre\\Data")

# This is just subsetting by chromosome 19, if you change the chrom == 19 to another number it will change
File <- read.delim("SNP-350k.hmp.txt", header = T)
head(File)
ChromNineteen <- subset(File, File$chrom == 19)
write.table(ChromNineteen, file = "ChromNineteen.txt", row.names = F, sep = "\t", quote = F)

# This function produces an LD plot - needs to be given TASSEL LD output file, where the 13th row is dist_bp, and the 14th row is R2
# To produce the input file, load the .hmp into TASSEL and then run LD on it to give an output
LDPlot <- function(txtfile){
  LD <- read.table(txtfile, header = T)
  Dist <- as.vector(LD[[13]])
  R2 <- as.vector(LD[[14]])
  
  R2_2 <- subset(R2, ((Dist != "N/A") == T) & ((R2 != "NaN") == T) & ((R2 != 1) == T))
  Dist_2 <- subset(Dist, ((Dist != "N/A") == T) & ((R2 != "NaN") == T) & ((R2 != 1) == T))
  
  R2_2 <- as.numeric(R2_2)
  Dist_2 <- as.numeric(Dist_2)
  
  ggplot() + 
    geom_point(aes(Dist_2/1000000, R2_2), na.rm = T, color = "darkgreen", alpha = 0.05) + 
    labs(x="Physical distance (Mb)", y=expression(LD~(R^{2})))
}

LDPlot("LD.txt")
LDPlot("ChromNineteenLD.txt")

