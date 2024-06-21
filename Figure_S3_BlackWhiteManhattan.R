##GWAS Manhattan white and black seed

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
  return(manhat)
}

###Load data and filter to add to supplement
white <- read.xlsx("Supplement_Seed_color_V1.xlsx",sheet="S5SeeMWhite_GWAS_Result_Sativa")
black <- read.xlsx("Supplement_Seed_color_V1.xlsx",sheet="S6SeeMBlack_GWAS_Result_Sativa")


##Make a plot for each gassoc result
bf <- -log10(0.05/1566489)
pseudo.points <- data.frame(cbind(rep(1:9,2),rep(2,18),c(rep(0,9),
                                                         214.8,217.1,257.8,377.4,339.6,193.1,195.5,309.6,203.9)))
colnames(pseudo.points) <- c("chr","pval","pos")

to.pl1 <- as.data.frame(cbind(white$P.value,white$Chr,white$Pos/1e6))
colnames(to.pl1) <- c("pval","chr","pos")
white.manhat <- make_manhat(to.pl1,bf,pseudo.points)
white.manhat <- white.manhat + theme(axis.title.x= element_blank())+ggtitle("White seed coat (Grey mean > 100|n=80)")
white.manhat
to.pl2 <- as.data.frame(cbind(black$P.value,black$Chr,black$Pos/1e6))
colnames(to.pl2) <- c("pval","chr","pos")
black.manhat <- make_manhat(to.pl2,bf,pseudo.points)
black.manhat <- black.manhat +ggtitle("Black seed coat (Grey mean < 100|n=49)")


combined_plots <- plot_grid(white.manhat,black.manhat, ncol = 1,labels = c("A","B"),
                            rel_heights = c(1,1.1))
combined_plots
ggsave("Figure_S3.pdf",combined_plots,width = 12,height = 5,dpi = 300,bg = "white")

