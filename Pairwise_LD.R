##Adgenet Pairwise LD plot
library(adegenet)
library(dartR)
library(ggplot2)
library(openxlsx)
##
## Load the data

###LOCUS 7 
gassoc_gls <-  read.xlsx("./Supplement_Seed_color_V1.xlsx",sheet = "Chr7_Locus_Genotypes")
rownames(gassoc_gls) <- gassoc_gls$Locus
gassoc_gls <- gassoc_gls[,-1]
gassoc_gls <- t(gassoc_gls-1)
##
x1 <- new("genlight", gassoc_gls,ploidy=2)


##Calculate LD
ld.calc <- gl.report.ld(
  x1,
  name = NULL,
  save = F,
  nchunks = 2,
  ncores = 95,
  verbose = 3)

save(ld.calc, file = "./LD_Loc7.RData")