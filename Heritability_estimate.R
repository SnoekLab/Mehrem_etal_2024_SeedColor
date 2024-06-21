##This scripts calculates the broad-sense heritability of seed coat color

library(openxlsx)
library(data.table)

rawdata <- read.xlsx("./Supplement_Seed_color_V1.xlsx",sheet = "S2Seed_Color_Measurements")
##Take the average of every two rows (length and width of one seed)
rawdata <- data.frame(
  TKIID = rawdata$TKIID[seq(1, nrow(rawdata), by = 2)],
  Mean = rowMeans(matrix(rawdata$Mean, ncol = 2, byrow = TRUE)))
  
  
acc <- read.xlsx("./Supplement_Seed_color_V1.xlsx",sheet = "S1PhenotypeCollection")

# sativa
sat.acc <- acc$TKIID[acc$Species_SB == "L. sativa"] # TKI-001 to TKI-141
sat.acc <- sat.acc[!is.na(sat.acc)]
seed.col <- rawdata$Mean[rawdata$TKIID %in% sat.acc] ; 
raw.acc <- rawdata$TKIID[rawdata$TKIID %in% sat.acc] ; unique(raw.acc)
sat.acc[!sat.acc %in% raw.acc]
meansq.collect <- anova(lm(seed.col~raw.acc))
meansq.collect

within <- meansq.collect$`Mean Sq`[2]
dfRun <- meansq.collect$Df[1]
dfWithin <- meansq.collect$Df[2]
geno <- meansq.collect$`Mean Sq`[1]
between <- (geno-within)/((dfWithin/(dfRun+1))+1) # (S1^2-S2^2)/J
total <- between+within
between # Between genotype Variance
within # Within genotype Variance

h2.sat <- round((between/total)*100,2)
h2.sat # 93.11


# serriola

ser.acc <- acc$TKIID[acc$Species_SB == "L. serriola"] # TKI-142 to TKI-340
ser.acc <- ser.acc[!is.na(ser.acc)]
seed.col <- rawdata$Mean[rawdata$TKIID %in% ser.acc]
raw.acc <- rawdata$TKIID[rawdata$TKIID %in% ser.acc] ; unique(raw.acc)
meansq.collect <- anova(lm(seed.col~raw.acc))
meansq.collect$Df
meansq.collect

within <- meansq.collect$`Mean Sq`[2]
dfRun <- meansq.collect$Df[1]
dfWithin <- meansq.collect$Df[2]
geno <- meansq.collect$`Mean Sq`[1]
between <- (geno-within)/((dfWithin/(dfRun+1))+1) # (S1^2-S2^2)/J
total <- between+within # (S1^2+(S2^2)
between # Between genotype Variance
within # Within genotype Variance

h2.ser <- round((between/total)*100,2)
h2.ser # 40.39


