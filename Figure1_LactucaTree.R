#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("ggtree")

##Script for Figure 1A-B
library(openxlsx)
library(ggplot2)
library(viridis)
library(ggtree)
library(cowplot)
library(ape)
library(gplots)
library(stringr)
library(png)
##Figure 1A
##Load png of seed colors
fig1a <- readPNG("./Figure_1A.png")
fig1a <- ggplot() + ggplot2::annotation_custom(grid::rasterGrob(fig1a,
                                            width=ggplot2::unit(1,"npc"),
                                            height=ggplot2::unit(1,"npc")),
                           -Inf, Inf, -Inf, Inf)

fig1a <- fig1a+theme(plot.margin = unit(c(0,0,0,0), "cm"))
fig1a

##Figure 1B

lactree <- read.xlsx("./Supplement_Seed_color_V1.xlsx", sheet = "Lactuca_tree",colNames = F)
lactree <- read.tree(text=lactree[1,1])


leaves_to_remove <- c("L.siberica", "L.tartarica", "L.viminea","L.perennis","L.tenerrima")
lactree <- ape::drop.tip(lactree,leaves_to_remove)
lactree$tip.label
lactree$tip.label
lactree$edge
lactree$Nnode
lactree$node.label
treeorder<- lactree$tip.label
str(lactree)
attr(lactree,"phylo")
attr(lactree,"cladewise")


use.col <- c("blue1","blue1","blue1","blue1","purple","purple",
             "purple","purple","red3",
             "red3","red3","red3",
             "black","black",
             "black","black")
use.nuc <- c(rep("NucA",12),rep("NucB",4))
mytree <- ggtree(lactree)+
  #geom_hilight(mapping = aes(subset = node %in% c(18),fill=c("NucA")), ###RECTANGLES show up but dont cover whole tree
               #alpha=.3,type = "auto",align="both")+
  #geom_hilight(mapping = aes(subset = node %in% c(28),fill=c("NucB")),
              # alpha=.3,type = "auto",align="both")+
  annotate(geom = "text",x=4,y=5,label="NucA")+
  annotate(geom = "text",x=4,y=4,label="NucB")+
  geom_tiplab(color="black",hjust = 0,align = F,offset = 1.25,fontface="italic")+
  geom_tippoint(color=use.col,size=5)+
  geom_strip('L.serriola', 'L.sagittata', barsize=1, color='blue1', label="Pool I",
             offset = 0.25, offset.text=0.4,angle = 90,hjust = 0.5) +
  geom_strip('L.georgica', 'L.saligna', barsize=1, color='purple', label="Pool II",
             offset = 0.25, offset.text=0.4,angle = 90,hjust = 0.5) + 
  geom_strip('L.quercina', 'L.indica', barsize=1, color='red1', label="Pool III",
             offset = 0.25, offset.text=0.4,angle = 90,hjust = 0.5) + 
  scale_x_continuous(expand = c(0,5))+
  #scale_y_continuous(expand = c(0,0))+
  
  theme(legend.position = "none")

mytree
treeorder <- c("L.serriola","L.altaica",
               "L.sativa","L.sagittata","L.georgica",
               "L.virosa",
               "L.aculeata","L.saligna","L.quercina",
               "L.raddeana","L.homblei","L.indica",
               "L.canadensis","L.biennis","L.floridana","L.palmensis") 
mytree <- mytree+
  geom_hilight(mapping = aes(subset = node %in% c(18,28),fill=c("NucA","NucB")),alpha=.3,type = "auto",align="left")+
   scale_fill_manual(values=c("steelblue", "darkgreen"))+
  ##no plot margings
  theme(plot.margin = unit(c(0,0,0,-10), "cm"))
mytree

#### complete figure can be found all the way down!


#############################################################
## Add Phenotypes
lkphe <-  read.xlsx("./Supplement_Seed_color_V1.xlsx",sheet="S1PhenotypeCollection")


lkphe
colnames(lkphe)
lkphe <- lkphe[!lkphe$Species_SB=="L. sativa x L. serriola",]

### Mean values and sd per species group
table(lkphe$Species_SB)


## set 1: Sativa, Serriola, Saligna, Virosa, Wild
## set 2: Sativa types
## set 3: all species 

unique(lkphe$Species_SB)
let.group <- lkphe$Species_SB


## Seed traits......-----------------------------------------------------------
## Seed color ; Seed Size ; Seed length 

## seed color (grey int)
aggregate(lkphe$SeeMeanGrey,list(lkphe$Species_SB),mean,na.rm=T)


colnames(lkphe)
to.pl <- data.frame(lkphe[,c("Species_SB","SeeMeanGrey")])
colnames(to.pl) <- c("Species","Seed_color")
to.pl$Species <- str_remove(to.pl$Species," ")
to.pl <- to.pl[to.pl$Species %in% lactree$tip.label,]
to.pl$Species <- factor(to.pl$Species,levels = rev(treeorder))
min.obs <- min(to.pl[,2],na.rm=T)
max.obs <- max(to.pl[,2],na.rm=T)
use.seq <- seq(min.obs,max.obs,length.out=301)
use.y <- rep(use.seq,each=length(treeorder))
use.species <- rep(treeorder,301)
bgrnd <- data.frame(use.species,use.y)
bgrnd$use.species <- factor(bgrnd$use.species,levels = rev(treeorder))

#add.spe <- c("L. sibirica","L. tatarica","L. viminea","L. perennis","L. tenerrima")
#aprox.seedcol <- c(NA,80,60,60,60)
#ads <- data.frame(add.spe,aprox.seedcol)
use.col <- rev(c(rep("blue",4),rep("purple",4),rep("red",4),rep("black",4)))
mytree
to.pl <- to.pl[!is.na(to.pl$Seed_color),]



seed_color_fig <- ggplot(to.pl,aes(Species,Seed_color))+
  geom_tile(data =  bgrnd,aes(x=use.species,y=use.y,fill=use.y))+
  geom_vline(xintercept = 1:21-0.5,col="white")+
  geom_violin(scale="width",col="grey",size=0.2)+
  stat_summary(fun = "mean", colour = use.col, size = 3, geom = "point")+
  #geom_point(data=ads,aes(add.spe,aprox.seedcol),col="grey70")+
  scale_fill_gradientn(colors = c("black","tan","ivory"))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0),limits=rev(treeorder))+
  coord_flip()+
  theme(legend.position = "none")+
  ylab("Seed coat color")

seed_color_fig

dev.off()

####Change Tree labels
summary_data<- as.data.frame(table(to.pl$Species))
summary_data <- summary_data[match(summary_data$Var1,treeorder),]
colnames(summary_data)<- c("Species","Freq")
summary_data$Label.new <- paste(summary_data$Species," (n=",summary_data$Freq,")",sep="")
mytree$data$label <- summary_data$Label.new[match(mytree$data$label,summary_data$Species)]
mytree
########################################################################3
##### Complete figure 



mf.theme <- theme(axis.title.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.ticks.x = element_line(size=0.2),
                  panel.border = element_rect(fill = NA,color = "black",size=0.2))
mf.theme
plot_grid(mytree,
          seed_color_fig+mf.theme,
          ncol=2,align = "h",rel_widths = c(1,1))

mf.theme <- theme(axis.title.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.ticks.x = element_line(size=0.2),
                  panel.border = element_rect(fill = NA,color = "black",size=0.2),
                  plot.margin = unit(c(0,0,0,0), "cm"))

fig1b <-plot_grid(mytree,
          seed_color_fig+mf.theme,
          ncol=2,align = "hv",axis ="tblr", rel_widths = c(3,1.5))

fig1b 

######Figure 1C


library(ggplot2)
library(openxlsx)

##LoadData

pheno <- read.xlsx("Supplement_Seed_color_V1.xlsx", sheet = "S1PhenotypeCollection")
pheno <- pheno[!is.na(pheno$Subgroup_SB),]

##Plot dots and violin based on pheno$seedMgray

pheno$Subgroup_SB <- factor(pheno$Subgroup_SB)
pheno <- pheno[!is.na(pheno$SeeMeanGrey),]
pheno <- pheno[,c(9,11)]
fig1c <- ggplot(pheno, aes(x=Subgroup_SB, y=SeeMeanGrey)) +
  geom_violin(scale = "count",fill = "grey", color = "grey", alpha = 0.7) +
  geom_jitter(width=0.3,aes(fill=SeeMeanGrey),pch=21,size=2)+
  scale_fill_gradientn(colors = c("black","tan","ivory"),name="Seed coat color")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=10),
        legend.position = "none") +
  labs(x="Subgroup", y="Seed coat color") +
  ##Remove grey lines from background
  theme(panel.background = element_rect(fill="white")) +
  ##Keep outline around plot 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),plot.margin = unit(c(0,0,0,0), "cm")) +
  guides(alpha=F)
fig1c 


###Use patchwork to combine all numbers and label them a to c
library(patchwork)


patch2 <- fig1a +fig1b +fig1c +plot_layout(design = "
              12
              32
              ",widths = c(1,2), height = c(1,1,4))+
  plot_annotation(tag_levels = 'A',
      theme = theme(plot.tag= element_text(size = 16),plot.tag.position =c ))
  
patch2

# Create the layout using cowplot's plot_grid function

ggsave(patch2,file="Figure_1.pdf",width=16,height=12,dpi=300)
## END #####################################################################################################
############################################################################################################