###Script for LD plot in Sativa Locus 7
library(ggplot2)
library(ldsep)
library(openxlsx)
library(dplyr)
library(stringr)
library(ape)
library(dartR)
library(ggplot2)
library(hscovar)
library(openxlsx)
library(ggborderline)
library(patchwork)


#######################
#Load gene annotation file

gene.anno <- read.gff("./GCF_002870075.3_Lsat_Salinas_v8_genomic.gff") #https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_002870075.3/
#filter gene anno for the first 9 unique seqIDs from .
gene.anno <- gene.anno[gene.anno$seqid %in% unique(gene.anno$seqid)[1:9],]
chr.conv <- cbind(unique(as.character(gene.anno$seqid)),1:9)
gene.anno$chromosome.number <- chr.conv[match(gene.anno$seqid,chr.conv[,1]),2]
gene.anno <- gene.anno[!gene.anno$type == "region",]

##Extract the gene name from info only and cbind to gene.anno

gassoc_gls <-  read.xlsx("./Supplement_Seed_color_V1.xlsx",sheet = "S3SeeMGrey_GWAS_Result_Sativa")

locus.chr7 <- read.xlsx("Supplement_Seed_color_V1.xlsx",sheet="Chr7_Locus_Genotypes",rowNames = T)

locus.chr7 <- locus.chr7[rownames(locus.chr7) %in% paste(paste(gassoc_gls$chr,gassoc_gls$pos,sep="_")),]
gassoc_gls.7 <- gassoc_gls[match(rownames(locus.chr7),paste(gassoc_gls$chr,gassoc_gls$pos,sep="_")),]
chromosome <- 7
pos.max <- gassoc_gls$pos[which.max(gassoc_gls$pval)]
boundaries <- 1e6
start.seq <- pos.max/1e6 - 1
end.seq <- pos.max/1e6 + 1

snp.to.pl <- as.data.frame(cbind(gassoc_gls.7$pval,gassoc_gls.7$chr,gassoc_gls.7$pos))
colnames(snp.to.pl) <- c("pval","CHR","POS")

library(viridis)
library(cowplot)

mhpl <- ggplot(snp.to.pl,aes(POS/1e6,pval))+
  geom_point(alpha=0.7)+
  geom_hline(yintercept=7.52, linetype='dotted', col = 'red',size=0.2)+
  
  xlab("Position (Mbp)") + ylab(bquote(-log[10](p))) +
  theme_cowplot()+
  xlim(start.seq,end.seq)+
  theme(panel.border = element_rect(colour = "black",linetype = "solid",linewidth = 0.2),
        axis.line = element_blank(),
        axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(linewidth = 0.2),
        text=element_text(size=15),
        panel.grid.major.x = element_line(linewidth = 0.2,color = "grey70"),
        legend.key.height = unit(0.8, "cm"),
        legend.position = c(0,1),
        legend.justification = c(-0.05,1),
        plot.margin = unit(c(0,0,0,0), "cm"))

mhpl

###Gene annotation plot
##Filter gene.anno by locus
locus.info <- gene.anno[gene.anno$chromosome.number == as.character(chromosome) & 
                          between(gene.anno$start,start.seq*1e6,end.seq*1e6),]

locus.info.topl <- locus.info[locus.info$type =="gene",]

anno.plot <- ggplot(locus.info.topl)+
  geom_segment(aes(x=start/1e6, xend=end/1e6, y=chromosome.number, yend=chromosome.number,color=type), 
               size=15,alpha=0.8)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5),limits = c(start.seq,end.seq)) +
  theme(text=element_text(size=15),
        legend.position = "none",
        legend.key.size=unit(0.8,"cm"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y =element_blank(),
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=10),
        plot.margin = unit(c(0,0,0,0), "cm"))

anno.plot


###Add LD plot

load("LD_Loc7.RData")
head(test)
###
#Filter SNP file
locus.chr7 <- read.xlsx("./Supplement_Seed_color_V1.xlsx",sheet="Chr7_Locus_Genotypes",rowNames = T)

sat.snps <- locus.chr7
##Separate rownames into CHR and POS of sat.snp with splitting at _
snp.info <- str_split(rownames(sat.snps),"_")
chr <- sapply(snp.info, function(x) x[1])
pos <- sapply(snp.info, function(x) x[2])
sat.snps$CHR <- as.numeric(chr)
sat.snps$POS <- as.numeric(pos)

# Get the unique loc1 and loc2 values
locs <- sort(unique(c(sat.snps$POS[test$loc1], sat.snps$POS[test$loc2])))

# Create an empty matrix to store the R^2 values
ld_matrix <- matrix(NA, nrow = length(locs), ncol = length(locs))

ld_matrix[cbind(test$loc1, test$loc2)] <- test$R2

dimnames(ld_matrix) <- list(locs, locs)


###Prune matrix
test_grouped <- tagSNP(ld_matrix, threshold = 0.90)
test_grouped_unlist <- as.numeric(unlist(rlist::list.select(test_grouped, tagsnp)))
##subset LD matrix using test grouped unlist
ld_matrix <- ld_matrix[as.character(sort(test_grouped_unlist)), as.character(sort(test_grouped_unlist))]


##Make symmetric matrix

# Create symmetric_matrix by adding ld_matrix to its transpose
symmetric_matrix <- ld_matrix + t(ld_matrix)

diag(symmetric_matrix) <- 0.9
snp.pos <- sat.snps$POS[match(sort(test_grouped_unlist),sat.snps$POS)]

library(reshape2)
ld_mat <- ld_matrix
ld <- melt(ld_mat,na.rm=TRUE)
#filter LD matrix such that x and y are between start.seq and end.seq
ld <- ld[ld$Var1 > start.seq*1e6 & ld$Var1 < end.seq*1e6 &
           ld$Var2 > start.seq*1e6 & ld$Var2 < end.seq*1e6 ,]
names(ld) <- c("x","y","r2")

ld$x1 <- ld$x+((ld$y-ld$x)/2)
ld$y1 <- ld$y-ld$x
ld$r2c <- cut(ld$r2,breaks=seq(0,1,0.2),labels=c("0-0 - 0.2","0.2 - 0.4",
                                                 "0.4 - 0.6","0.6 - 0.8",
                                                 "0.8 - 1.0"),include.lowest = T)
ld$r2c <- factor(ld$r2c,levels=rev(c("0-0 - 0.2","0.2 - 0.4",
                                     "0.4 - 0.6","0.6 - 0.8",
                                     "0.8 - 1.0")))


ld_plt <- ggplot(ld,aes(x=x1/1e6,y=y1/1e6,col=r2c))+
  geom_point(shape=20,size=0.1,alpha=1)+
  scale_color_manual(values=c("#ca0020","#f4a543","#d1e5f0","#67a9cf","#2166ac"))+
  scale_x_continuous(position = "top",limits = c(start.seq,end.seq))+
  scale_y_reverse(limits=rev(c(min(ld$y1/1e6),max(ld$y1/1e6))))+
  guides(colour=guide_legend(title="R2",override.aes=list(shape=20,size=8)))+
  theme_bw(base_size=14)+
  theme(panel.border=element_blank(),
        axis.ticks=element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x= element_blank(),
        axis.title = element_blank(),
        ##Remove margins
        plot.margin = margin(0,0,0,0),
        panel.grid.major.x = element_line(linewidth = 0.2,color = "grey70"),
        panel.grid.major.y = element_blank()
        
        
        
  )

ld_plt

##Combine both plots anno.plot and mhpl using patchwork

# Combine plots vertically
combined_plots <- mhpl/anno.plot/ld_plt

# Adjust the margin to make space for the legend
combined_plots <- combined_plots + plot_layout(ncol = 1,heights = c(2, 0.3,3))
combined_plots
# Print the combined plots
ggsave(combined_plots,file="FigureS2.pdf",width=10,height=8,dpi=300)
