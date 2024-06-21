###GAPIT Model for serriola GWAS

library(GAPIT)
library(data.table)



##Load serriola SNP data

load("pheno.ser.out") ##phenotype object
load("ser.snps.myGD.out")##snp object

##Make phenotypes look like MyY 
myY <- t(phenotypes.ser)
myY <- as.data.frame(cbind(rownames(myY),myY))
colnames(myY)[1] <- c("Taxa")
rownames(myY) <- NULL
##columns 2:5 need to be numeric
myY[,2:5] <- apply(myY[,2:5],2,as.numeric)
print(colnames(myY))

#GWAS
print("we run GWAS")
myGAPIT=GAPIT(
  Y=myY[,c(1,2)], #fist column is ID
  GD=myGD,
  GM=myGM,
  PCA.total=3,
  SNP.MAF=0.05,
  model=c("MLMM","FarmCPU", "Blink"),
  Multiple_analysis=TRUE)

myGAPIT=GAPIT(
  Y=myY[,c(1,2)], #fist column is ID
  GD=myGD,
  GM=myGM,
  PCA.total=3,
  SNP.MAF=0.05,
  model=c("SUPER"),
  Multiple_analysis=TRUE)