####Script for Plotting Figure 2
####LIBS
library(ggplot2)
library(facetscales)
library(cowplot)
library(ggfastman)

#######################Figure 2
#Load Data Sativa
gassoc_gls <-  read.xlsx("./Supplement_Seed_color_V1.xlsx",sheet = "S3SeeMGrey_GWAS_Result_Sativa")
####Plot Sativa only
### without facet scales, with pseudo points.....

to.pl <- as.data.frame(cbind(gassoc_gls$pval,gassoc_gls$chr,gassoc_gls$pos/1e6))
colnames(to.pl) <- c("pval","chr","pos")
to.pl$pval <- as.numeric(to.pl$pval)
to.pl$pos <- as.numeric(to.pl$pos)
to.pl <- to.pl[complete.cases(to.pl),]
to.pl <- to.pl[to.pl$pval >2,]
bf <- -log10(0.05/1566489)
pseudo.points <- data.frame(cbind(rep(1:9,2),rep(2,18),c(rep(0,9),
                      214.8,217.1,257.8,377.4,339.6,193.1,195.5,309.6,203.9)))
colnames(pseudo.points) <- c("chr","pval","pos")
species <- rep("L. sativa",nrow(to.pl))
threshold <- rep(7.5,nrow(to.pl))
to.pl.all <- cbind(to.pl,species,threshold)
fig2a <- ggplot(to.pl,aes(pos,pval))+
  geom_point(shape=19,alpha=0.4)+
  geom_point(data = pseudo.points,aes(pos,pval),alpha=0)+
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

fig2a

#######Figure 2B
####Load data Sativa
phenotypes <-  read.xlsx("./Supplement_Seed_color_V1.xlsx",sheet = "S1PhenotypeCollection")
cofac.sat <-  read.xlsx("./Supplement_Seed_color_V1.xlsx",sheet = "Chr7_Sativa_TopSNP")
cofac <- as.numeric(cofac.sat$Genotype)
names(cofac) <- cofac.sat$BGI_ID
pheno <- phenotypes$SeeMeanGrey[match(toupper(names(cofac)),phenotypes$BGIID)]
to.pl <- as.data.frame(cbind(pheno,cofac))
to.pl[,2] <- as.character(to.pl[,2])
color_palette <- c( "#000000", "#8B7355","#F5F5DC")
to.pl$cofac[to.pl$cofac == 1] <- "REF:AA"
to.pl$cofac[to.pl$cofac == 2] <- "REF/ALT:AG"
to.pl$cofac[to.pl$cofac == 3] <- "ALT:GG"
to.pl$cofac <- factor(to.pl$cofac,levels=c("REF:AA","REF/ALT:AG","ALT:GG"))



fig2b <- ggplot(to.pl, aes(x=cofac, y=pheno)) + 
  geom_boxplot(notch=F,outlier.shape = NA,fill="grey",linewidth=0.5,)+
  geom_point(aes(y = pheno,fill=pheno,color=pheno), position = position_jitter(width = 0.1),size=4,shape=21)+
  scale_fill_gradientn(colors=color_palette)+
  scale_color_gradientn(colors="black",guide ="none" )+
  ylab("Seed coat color")+
  xlab("Allele status")+
  ggtitle("Chr 7| 50 317 920 bp")+
  theme_cowplot()+
  # annotate("label",
  # x = 1:length(table(to.pl$cofac)),
  # y = aggregate(pheno ~ cofac, to.pl, median)[ , 2],
  # label = table(to.pl$cofac),
  # col = "black",
  #fontface="bold")+
  theme(plot.title = element_text(face = "bold"),legend.position = "none")+
  labs(fill="Color")

fig2b

##Combine plots
fig2 <- plot_grid(fig2a,fig2b,labels = c("A","B"),nrow=1,rel_widths=c(3,1))
fig2
ggsave("./Figure_2.pdf",fig2,width=12,height=3,dpi = 300,bg = "white")
dev.off()


#################################Figure 3


# ####Load Serriola results

gassoc_gls <-  read.xlsx("./Supplement_Seed_color_V1.xlsx",sheet = "S7SeeMGrey_GWAS_Result_Serriola")


to.pl <- as.data.frame(cbind(gassoc_gls$pval,gassoc_gls$chr,gassoc_gls$pos/1e6))
colnames(to.pl) <- c("pval","chr","pos")
to.pl$pval <- as.numeric(to.pl$pval)
to.pl$pos <- as.numeric(to.pl$pos)
to.pl <- to.pl[complete.cases(to.pl),]
to.pl <- to.pl[to.pl$pval >2,]
bf <- -log10(0.05/3214826)
pseudo.points <- data.frame(cbind(rep(1:9,2),rep(2,18),c(rep(0,9),
                                                         214.8,217.1,257.8,377.4,339.6,193.1,195.5,309.6,203.9)))
colnames(pseudo.points) <- c("chr","pval","pos")
species <- rep("L. sativa",nrow(to.pl))
threshold <- rep(7.5,nrow(to.pl))
to.pl.all <- cbind(to.pl,species,threshold)
fig3 <- ggplot(to.pl,aes(pos,pval))+
  geom_point(shape=19,alpha=0.4)+
  geom_point(data = pseudo.points,aes(pos,pval),alpha=0)+
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

fig3

ggsave("Figure_3.pdf",fig3,width=12,height=3,dpi=300,bg="white")

#####




