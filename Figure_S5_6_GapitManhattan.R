##Plot for GAPIT GWAS
library(ggplot2)
library(openxlsx)
library(cowplot)
library(patchwork)
#devtools::install_github("roman-tremmel/ggfastman", upgrade = "never")
library(ggfastman)

###Function make manhattan plot
make_manhat <- function (x,bf,pseudo.points){
  manhat<- ggplot(x,aes(pos,pval))+
    geom_point(shape=19,alpha=0.4)+
    geom_point(data = pseudo.points,aes(pos,pval),alpha=0)+
    geom_hline(yintercept=bf, linetype='dotted', col = 'red',size=1)+
    scale_x_continuous(breaks = c(50,100,150,200,250,300,350),expand = c(0,0))+
    facet_grid(.~chr,space="free_x",scale="free_x")+
    xlab("Position (Mbp)") + ylab(bquote(-log[10] (P)))+
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
  return(manhat)
}

##Load from GWAS results for QQs
super.r <- read.csv("Serriola_Results/GAPIT_SeeMeanGrey/GAPIT.Association.GWAS_Results.SUPER.SeeMeanGrey.csv")
farmcpu.r <- read.csv("Serriola_Results/GAPIT_SeeMeanGrey/GAPIT.Association.GWAS_Results.FarmCPU.SeeMeanGrey.csv")
mlmm.r <- read.csv("Serriola_Results/GAPIT_SeeMeanGrey/GAPIT.Association.GWAS_Results.MLMM.SeeMeanGrey.csv")
blink.r <- read.csv("Serriola_Results/GAPIT_SeeMeanGrey/GAPIT.Association.GWAS_Results.BLINK.SeeMeanGrey.csv")


####Load results from supplement for manhattans
gassoc_gls.mlmm <-  read.xlsx("./Supplement_Seed_color.xlsx",sheet = "S11SeeMGrey_Ser_MLMM")
gassoc_gls.farmcpu <-  read.xlsx("./Supplement_Seed_color.xlsx",sheet = "S8SeeMGrey_Ser_FarmCPU")
gassoc_gls.super <-  read.xlsx("./Supplement_Seed_color.xlsx",sheet = "S10SeeMGrey_Ser_SUPER")
gassoc_gls.blink <-  read.xlsx("./Supplement_Seed_color.xlsx",sheet = "S9SeeMGrey_Ser_Blink")


##Make a plot for each gassoc result
bf <- -log10(0.05/3633591) #Number of SNPs tested
pseudo.points <- data.frame(cbind(rep(1:9,2),rep(2,18),c(rep(0,9),
                                                         214.8,217.1,257.8,377.4,339.6,193.1,195.5,309.6,203.9)))
colnames(pseudo.points) <- c("chr","pval","pos")

to.pl1 <- as.data.frame(cbind(gassoc_gls.farmcpu$P.value,gassoc_gls.farmcpu$Chr,gassoc_gls.farmcpu$Pos/1e6))
colnames(to.pl1) <- c("pval","chr","pos")
farmcpu_manhat <- make_manhat(to.pl1,bf,pseudo.points)



to.pl2 <- as.data.frame(cbind(gassoc_gls.mlmm$P.value,gassoc_gls.mlmm$Chr,gassoc_gls.mlmm$Pos/1e6))
colnames(to.pl2) <- c("pval","chr","pos")
mlmm_manhat <- make_manhat(to.pl2,bf,pseudo.points)
mlmm_manhat <- mlmm_manhat + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
##Combine the plots and give each a title using patchwork in one column
# Add labels to the plots
farmcpu_manhat <- farmcpu_manhat +theme(
  axis.title.x = element_blank(),axis.text.x = element_blank())


mlmm_manhat <- mlmm_manhat + theme(axis.text.x = element_blank())
##Make QQ plot for each gwas result read from csv
farmqq <- fast_qq(pvalue =farmcpu.r$P.value, speed = "fast")+labs(title = "FarmCPU")
mlmmqq <- fast_qq(pvalue = mlmm.r$P.value, speed = "fast")+labs(title = "MLMM")
# Combine plots vertically and display


##Super and blink go to supplement
to.pl3 <- as.data.frame(cbind(gassoc_gls.blink$P.value,gassoc_gls.blink$Chr,gassoc_gls.blink$Pos/1e6))
colnames(to.pl3) <- c("pval","chr","pos")
blink_manhat <- make_manhat(to.pl3,bf,pseudo.points)
blink_manhat <- blink_manhat+theme(axis.title.x = element_blank(),axis.text.x = element_blank())

blinkqq <- fast_qq(pvalue = blink.r$P.value, speed = "fast",)+labs(title = "BLINK")
blinkqq

to.pl4 <- as.data.frame(cbind(gassoc_gls.super$P.value,gassoc_gls.super$Chr,gassoc_gls.super$Pos/1e6))
colnames(to.pl4) <- c("pval","chr","pos")
super_manhat <- make_manhat(to.pl4,bf,pseudo.points)
super_manhat <- super_manhat
superqq <- fast_qq(pvalue = super.r$P.value, speed = "fast")+labs(title = "SUPER")
superqq

combined_plots <- plot_grid(farmcpu_manhat,mlmm_manhat,blink_manhat,super_manhat, ncol = 1,
                            labels = c("A","B","C","D"),label_size = 15,
                            rel_heights = c(1,1.1,1,1.4))
combined_plots
ggsave("Figure_S5.png",combined_plots,width = 16,height = 8,dpi = 300,bg = "white")




combined_qqs <- plot_grid(farmqq,mlmmqq,blinkqq,superqq, ncol = 2,label_size=15,labels=c("A","B","C","D"))
combined_qqs
ggsave("Figure_S6.png",combined_qqs,width = 16,height = 8,dpi = 300,bg = "white")

