


### Load devtools
library(devtools)
### Install lme4qtl from github
install_github("variani/lme4qtl")


## load needed packages
library(Matrix)
library(MASS)
library(ggplot2)
library(openxlsx)
# for Step 1: linear mixed model (no SNPs)
library(lme4qtl) 
# for Step 2: association tests
library(matlm) # use custom bug fixed matlm_pred function(s)
library(wlm) 
# for Step 3: explore GWAS results
library(cowplot)
library(gplots)

################################################################################
##
## Data 
##
################################################################################

### load r object with Lettuce (L. sativa) SNPs
load(file="obj_unisnp_num.out") # just 7 MB so can be supplemental
## inspect unisnp
unisnp[1:5,]
unisnp[1:5,10:ncol(unisnp)]
dim(unisnp)
snp.info <- unisnp[,1:10]
snp.info
#### phenotype data [[ <-- link to supplement table]]

pheno <- read.xlsx(xlsxFile = "Supplement_Seed_color_draft.xlsx",startRow = 2,sheet = "Zhang_2017_data",rowNames = T)
pheno$Seed.color_CGN_add

## get the values for which also genotype data exists
trait <-pheno[colnames(unisnp)[10:ncol(unisnp)],"Seed.color_CGN_add"]
trait <- as.numeric(trait == "white")
selc <- !is.na(trait)

use.trait <- trait[selc]
names(use.trait) <- colnames(unisnp)[10:ncol(unisnp)][selc]


####################### GWAS ################################################################


# prep and check snp matrix

usemat <- t(apply(t(unisnp[,10:ncol(unisnp)]),1,as.numeric))[selc,]
usemat[,1:5]
usemat[!usemat==-1] <- 0
usemat[usemat == -1] <- 2
dim(usemat) # 138 geno ; 155 741 SNPs
length(trait)
### calculate covariance matrix to be used as kinship matrix

letkin <- cov(t(usemat))
rownames(letkin) 
colnames(letkin)

threshold <- round(nrow(usemat)*0.05, digits=0)
##Filter for MAF
#selc2 <- apply(usemat==0,2,sum) > threshold & apply(usemat==1,2,sum) > threshold
#usemat <- usemat[,selc2]
#snp.info <- snp.info[selc2,]

### start mapping by making decomposition 
ID <- rownames(letkin) ; length(ID)
cbind(ID,use.trait)

####Compute the PCS
pca <- prcomp(usemat+1, scale = F)
pcs <- as.data.frame(pca$x[, 1:5])

mod <- lme4qtl::relmatLmer(use.trait ~ pcs$PC1 + pcs$PC2 + pcs$PC3 + pcs$PC4 + pcs$PC5 + (1|ID), relmat = list(ID = letkin))
V <- lme4qtl::varcov(mod, idvar = "ID")
Matrix::image(V[1:20, 1:20], main = "Estimated V (with artifacts)") # some artifacts close to zero due to limited numeric precision
V_thr <- V
V_thr[abs(V) < 1e-10] <- 0
Matrix::image(V_thr[1:20, 1:20], main = "Estimated V (with artifacts removed)")
decomp <- wlm::decompose_varcov(V, method = "evd", output = "all")
W <- decomp$transform

## total no. of SNPs for GWAS
sum(table(unisnp$X.CHROM))
## make data object for mapping without any extra factors
nnn <- rep(1,ncol(letkin))

### GWAS w/o kinship
# gassoc_lm <- matlm::matlm(use.trait ~ nnn,nnn, pred = usemat, ids = ids,batch_size = 2000, verbose = 2)
# plot(-log10(gassoc_lm$tab$pval))

### GWAS with kinship
gassoc_gls <- matlm::matlm(use.trait ~ nnn, nnn, pred =  usemat, ids = rownames(W), transform = W,batch_size = 2000, verbose = 2,cores = 1)
lod <- -log10(gassoc_gls$tab$pval)

# save(gassoc_gls,file="Zhang_GWAS_seed_col_outp_2023.out")

##### Manhattan plots with ggplot

p <- lod
pos <- snp.info$POS/1e6
chr <- snp.info$X.CHROM
# adjp <- -log10(p.adjust(gassoc_gls$tab$pval,"BH"))
which.max(p)
to.pl <- data.frame(chr,pos,p)
dim(to.pl)

# to.pl[mrkno,]
## get max and min pos per chr
mm.chr <- (aggregate(snp.info$POS,list(snp.info$X.CHROM),range))
mm.chr <- data.frame(mm.chr$Group.1,mm.chr$x)
to.pl[snp.info$POS %in% c(mm.chr[,2],mm.chr[,3]),"p"] <- 2
to.pl <- to.pl[to.pl$p>2 | snp.info$POS %in% c(mm.chr[,2],mm.chr[,3]),]

bf <- -log10(0.05/nrow(snp.info))
manhat<- ggplot(to.pl,aes(pos,p))+
  geom_point(shape=19,alpha=0.4)+
  geom_hline(yintercept=bf, linetype='dotted', col = 'red',size=1)+
  scale_x_continuous(breaks = c(50,100,150,200,250,300,350),expand = c(0,0))+
  facet_grid(.~chr,space="free_x",scale="free_x")+
  xlab("Position (Mbp)") + ylab(bquote(-log[10] (p)))+
  theme_cowplot() +
  theme(panel.background = element_rect(size=0.2,color="black"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth =0.2),
        axis.text.x = element_text(angle = 90,size=12,hjust = 1,vjust = 0.5),
        axis.text.y = element_text(lineheight =8,size=12),
        panel.spacing.x = unit(1,"mm"),
        panel.spacing.y = unit(1,"mm"),
        strip.text.y = element_text(face="bold", size=9),
        panel.grid.major.y = element_line(color = "grey2",
                                          size = 0.1,
                                          linetype = 3))
manhat
ggsave("Figure_S4.pdf",manhat,width = 12,height = 5,dpi = 300,bg = "white")




########################## END ###################################################################################