# brassica_napus_gwas
Code from my JIC summer internship with Professor Richard Morris

EulerandRungeKutta.R and MonteCarloPiEstimation.R are scripts I wrote in the first week of the internship, which introduced me to some mathematical/statistical methods.

GAPIT.R is the main script of my project, involving loading in data, testing trait data for normality, using LMMs to estimate trait value means, selecting the best model fit, and running GWAS analysis with GAPIT on several phenotype data sets.

LD Applied.R includes LD analysis and plotting of the genotype data, as formatted in TASSEL. The first couple of sections are attempts made at using certain packages which did not work. The working code begins on line 95, where a .hmp.txt file generated by TASSEL from genotype data is loaded into R. After being subsetted by chromosome, this is passed to the LDPlot function to be plotted.
