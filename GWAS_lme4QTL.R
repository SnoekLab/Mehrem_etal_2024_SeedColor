#  /$$       /$$      /$$ /$$$$$$$$ /$$   /$$             /$$     /$$        /$$$$$$  /$$      /$$  /$$$$$$   /$$$$$$ 
# | $$      | $$$    /$$$| $$_____/| $$  | $$            | $$    | $$       /$$__  $$| $$  /$ | $$ /$$__  $$ /$$__  $$
# | $$      | $$$$  /$$$$| $$      | $$  | $$  /$$$$$$  /$$$$$$  | $$      | $$  \__/| $$ /$$$| $$| $$  \ $$| $$  \__/
# | $$      | $$ $$/$$ $$| $$$$$   | $$$$$$$$ /$$__  $$|_  $$_/  | $$      | $$ /$$$$| $$/$$ $$ $$| $$$$$$$$|  $$$$$$ 
# | $$      | $$  $$$| $$| $$__/   |_____  $$| $$  \ $$  | $$    | $$      | $$|_  $$| $$$$_  $$$$| $$__  $$ \____  $$
# | $$      | $$\  $ | $$| $$            | $$| $$  | $$  | $$ /$$| $$      | $$  \ $$| $$$/ \  $$$| $$  | $$ /$$  \ $$
# | $$$$$$$$| $$ \/  | $$| $$$$$$$$      | $$|  $$$$$$$  |  $$$$/| $$      |  $$$$$$/| $$/   \  $$| $$  | $$|  $$$$$$/
# |________/|__/     |__/|________/      |__/ \____  $$   \___/  |__/       \______/ |__/     \__/|__/  |__/ \______/ 
#                                                  | $$                                                               
#                                                  | $$                                                               
#                                                  |__/     

#####################Dependencies####################################
## load needed packages
library(Matrix)
#library(colr)
library(MASS)
#library(ggplot2)
library(openxlsx)
# for Step 1: linear mixed model (no SNPs)
library(lme4qtl) 
# for Step 2: association tests
library(matlm) 
library(wlm) 
# for Step 3: explore GWAS results
#library(qq)
#library(cowplot)
#library(gplots)
#library(plyr)
#library(dplyr)
#library(parallel)
options(scipen = 999)
#library(scales)
#####################Functions#######################################

##Inputs
#genotypes: Matrix with columns being lines and rows being the SNPs
#trait: A vector with observances, binary or quantitative, with the same colnames and order as genotypes
#accession.type: A vector with respective Group of lines (e.g. virosa,saligna), same order as genotypes and trait cols
#phenotype.name: String, can be number or name
#snp.info: snp info in from of 3 columns (CHR,POS,QUAL)
#kinship: Kinship matrix (cov-var matrix of snp object)
#make.figures: Boolean, if TRUE the figures for kinship heatmap and covariance matrix are generated and saved
#out.dir: Directory in which to store the figures

##Outputs
#Within the function two figures: Kinship matrix and Covariance matrix
#List with GWAS results: gassoc_gls


GWAS <- function(genotypes, trait, phenotype.name,snp.info, kinship, out.dir) {
  phenotype <- toString(phenotype.name)
  letkin <- kinship
  usemat <- genotypes
  print(trait)
  selc <- !is.na(trait) #Selects the lines with an observation (removes lines that have an NA)
  trait.names <- trait[selc]
  use.trait <- trait[selc]
  print("Traits are selected.")
  print(paste("The phenotype ID is ", phenotype,".", sep=""))
 

  ## Filter in usemat and kinship object for mapping
  print(selc)
  print("THIS")
  print(nrow(usemat))
  usemat <- usemat[selc,]

  letkin <- letkin[selc,selc]
  ## Filter again MAF 5%
  threshold <- round(nrow(usemat)*0.05, digits=0)
  selc2 <- apply(usemat==1,2,sum) > threshold & apply(usemat==3,2,sum) > threshold
  print("SNPs falling within MAF >= 5%")
  print(table(selc2))
  usemat <- usemat[,selc2]
  #Save SNP positions from final SNP set
  
  chr <- snp.info[selc2,1]
  pos <- snp.info[selc2,2]
  
  
  
  print("Genotype matrix filtered and transformed.")
  
  ####Compute the PCS
  pca <- prcomp(usemat, scale = TRUE)
  pcs <- as.data.frame(pca$x[, 1:5])
  ### start mapping by making decomposition
  ID <- rownames(letkin) ; length(ID)
  cbind(ID,use.trait)
  mod <- lme4qtl::relmatLmer(use.trait ~ pcs$PC1 + pcs$PC2 + pcs$PC3 + pcs$PC4 + pcs$PC5 + (1|ID), relmat = list(ID = letkin))
  
  ##Calculate heritability
  herit.mod <- lme4qtl::VarProp(mod)
  V <- lme4qtl::varcov(mod, idvar = "ID")
  V_thr <- V
  V_thr[abs(V) < 1e-10] <- 0
  decomp <- wlm::decompose_varcov(V, method = "evd", output = "all")
  W <- decomp$transform
  print("Decomposition of covariance matrix was performed.")
 
  ## make data object for mapping without any extra factors
  nnn <- rep(1,ncol(letkin))

  ### GWAS with kinship
  gassoc_gls <- matlm::matlm(use.trait ~ nnn, nnn, pred =  usemat, ids = rownames(W), transform = W, batch_size = 2000, verbose = 2,cores = 1,stats_full = T)
  max.snp <- which.max(-log10(gassoc_gls$tab$pval))
  cofacsnp <- usemat[,max.snp]
  #GWAS with fixed
  gassoc_gls2 <- matlm::matlm(use.trait ~ nnn+cofacsnp, nnn, pred =  usemat, ids = rownames(W), transform = W, batch_size = 2000, verbose = 2,cores = 1,stats_full = T)
  
  gassoc_gls$trait.name <- phenotype
  gassoc_gls$trait.noobs <- length(use.trait)
  mrkno <- which.max(-log10(gassoc_gls$tab$pval))


  ###Save as integers
  gassoc_gls$tab$pval <- as.integer(round(-log10(gassoc_gls$tab$pval)*100000,0))
  gassoc_gls$tab$zscore <- as.integer(round(gassoc_gls$tab$zscore*100000,0))
  gassoc_gls$tab$se <- as.integer(round(gassoc_gls$tab$se*100000,0))
  gassoc_gls$tab$b <- as.integer(round(gassoc_gls$tab$b *100000,0))
  gassoc_gls <- as.data.frame(gassoc_gls$tab)
  save(gassoc_gls,file=paste(out.dir,"/GWAS_result_",phenotype,".out",sep=""))
  save(gassoc_gls2,file=paste(out.dir,"/GWAS2_result_",phenotype,".out",sep=""))
  

  save(chr,file=paste(out.dir,"/GWAS_chromosome.out",sep=""))
  save(pos,file=paste(out.dir,"/GWAS_position.out",sep=""))
  cofac <- usemat[,mrkno]
  save(cofac,file=paste(out.dir,"/GWAS_cofac.out",sep=""))
  print("Results saved.")

}

########################START SCRIPT#########################
#setwd("/home/sarah/sarah2/LettuceKnow/BGI_Metabo_GWAS/")
ext.dir <- "/home/sarah/sarah2/LettuceKnow/GWAS_SeedColor_2023/"
main.dir <- paste(ext.dir,"/GWAS_Results/",sep="")
dir.create(main.dir)


###INPUT


#### Load phenotype data
pheno <- load("/home/sarah/sarah2/LettuceKnow/GWAS_SeedColor_2023/data/seed.color.sat.out")
pheno <- eval(parse(text=pheno))
colnames(pheno) <- tolower(colnames(pheno))
base.dir <- paste(main.dir, "BGI_",sep="") ##Indicate which SNPs were used
print(colnames(pheno))
#Load genotype object for GWAS mapping
usemat <- load("/home/sarah/sarah2/LettuceKnow/BGI_GWAS/data/obj_bgi.sat.snps.out")
usemat <- eval(parse(text=usemat))
snp.info <- usemat[,c(1,2)]
usemat <- usemat[,c(-1,-2,-3)]
#usemat <- usemat[,!(colnames(usemat) %in% c("S053","S130"))]
#print(colnames(usemat))
usemat[usemat == 9] <- 2
usemat <- as.matrix(t(usemat))
print(rownames(usemat))
print(dim(usemat))
print("Above is before modifying")
#Load kinship matrix
load("/home/sarah/sarah2/LettuceKnow/BGI_GWAS/data/BGI_Sat_kinship.out")
###Input phenotype
i <- commandArgs(trailingOnly = T)
i <- as.numeric(i)
trait <- rownames(pheno)[i] 
pheno <- pheno[,!(colnames(pheno) %in% c("s053","s130","s131","s132"))]
new.dir <- paste(base.dir,trait,sep="")
dir.create(new.dir)

log.file.w <- file(paste(new.dir,"/","BGI_",trait,"_warning.log",sep=""),open="wt")
sink(file=log.file.w,type="message")
print("Pheno names here:")
print(names(pheno[i,]))
letkin_in <- letkin[tolower(names(pheno[i,])),tolower(names(pheno[i,]))] #In case we do not have information for all lines with this phenotype
usemat_in <- usemat[tolower(names(pheno[i,])),] #In case we do not have information for all lines with this phenotype
print(dim(usemat_in)) 
GWAS(genotypes = usemat_in, trait = as.vector(pheno[i,]), phenotype.name = trait, kinship=letkin_in, snp.info = snp.info, out.dir=new.dir)

sink(type="message")

close(log.file.w)

print(paste("GWAS finished. Phenotype is ",trait,sep=""))
#print(paste(nrow(pheno) - i," traits of ",nrow(pheno)," to go.",sep=""))




