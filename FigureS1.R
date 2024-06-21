##Script for manhattan plot with cofactor in L. sativa (Fig S1)
library(openxlsx)
library(ggplot2)
library(cowplot)
#Load Data Sativa
gassoc_gls <-  read.xlsx("../Supplement_Seed_color_V1.xlsx",sheet = "S3_SeeMGrey_GWAS_Result_Sativa")
gassoc_gls2 <-  read.xlsx("../Supplement_Seed_color_V1.xlsx",sheet = "SeeMGreyCofac_Result_Sat")
##vector of colors for gassoc and gassoc 2

####Plot Sativa only
### without facet scales, with pseudo points.....

to.pl <- as.data.frame(cbind(gassoc_gls$pval,gassoc_gls$chr,gassoc_gls$pos/1e6))
to.pl2 <- as.data.frame(cbind(gassoc_gls2$pval,gassoc_gls2$chr,gassoc_gls2$pos/1e6))
colnames(to.pl) <- c("pval","chr","pos")
colnames(to.pl2) <- c("pval","chr","pos")
to.pl$pval <- as.numeric(to.pl$pval)
to.pl2$pval <- as.numeric(to.pl2$pval)
to.pl$pos <- as.numeric(to.pl$pos)
to.pl2$pos <- as.numeric(to.pl2$pos)
to.pl <- to.pl[complete.cases(to.pl),]
to.pl2 <- to.pl2[complete.cases(to.pl2),]
to.pl <- to.pl[to.pl$pval >2,]
to.pl2 <- to.pl2[to.pl2$pval >2,]
bf <- -log10(0.05/1566489)
pseudo.points <- data.frame(cbind(rep(1:9,2),rep(2,18),c(rep(0,9),
                                                         214.8,217.1,257.8,377.4,339.6,193.1,195.5,309.6,203.9)))
colnames(pseudo.points) <- c("chr","pval","pos")
species <- rep("L. sativa",nrow(to.pl))
threshold <- rep(7.5,nrow(to.pl))
to.pl.all <- cbind(to.pl,species,threshold)
figs1 <- ggplot(to.pl,aes(pos,pval))+
  geom_point(shape=19,alpha=0.4)+
  geom_point(data = to.pl2,aes(pos,pval),shape=19,alpha=0.4,color="red")+
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

figs1

ggsave("../Figure_S1.pdf",figs1,width=12,height=3,dpi=300,bg="white")
