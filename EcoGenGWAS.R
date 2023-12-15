# Doing the the correlation matrix for BLUPS in VT and NC
# Jasmine 12/06/23

# reading in files

######## bash ##################

# modified GEA script to run the GWAS final

#mkdir ~/myresults/RSGWAS

REF="/netfiles/ecogen/PopulationGenomics/ref_genome/Pabies1.0-genome_reduced.fa"

SUFFIX="nc_s"

#_s means that pop structure is accounted for  

# Let's start with bio10; can do others as time permits...
BIOVAR="budbreak_blups_nc.txt"

# path to the red spruce bam files
#INPUT="/netfiles/ecogen/PopulationGenomics/fastq/red_spruce/cleanreads/bam340"

OUTPUT=~/myresults/ANGSD_GWAS2

# make the bamlist files
#ls ${INPUT} >${OUTPUT}/${SUFFIX}_bam.list


# Run ANGSD to estimate the genotype probabilities and perform the GEA:

ANGSD -b ${OUTPUT}/G2_budbreak_bam.list \
-ref ${REF} -anc ${REF} \
-out ${OUTPUT}/${SUFFIX}  \
-nThreads 1 \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 20 \
-minQ 20 \
-GL 1 \
-doCounts 1 \
-minInd 47 \
-setMinDepthInd 1 \
-setMaxDepthInd 40 \
-skipTriallelic 1 \
-doMajorMinor 1 \
-doMaf 1 \
-SNP_pval 1e-6 \
-minMaf 0.05 \
-doPost 1 \
-doAsso 5 \
-yQuant /data/users/j/h/jhart12/mydata/${BIOVAR} \
-cov /data/project_data/budbreak/allGL_genPC1_2.txt

# gwas with pop structure accounted PCangsd and input at -cov
######## bash ##################

vt <- read.delim(file.choose("budbreak_blups_vt.txt"), header = FALSE)

nc <- read.delim(file.choose("budbreak_blups_nc.txt"), header = FALSE)

# merging the two BLUPS
cb <- cbind(vt, nc)

# renaming the blups columns
names(cb)[1] <- "VT"
names(cb)[2] <- "NC"

# computing correlation matrix 
res <- cor(cb)

res

# importing latitude and longitude

library(ggplot2) # plotting
library(ggpubr) # plotting
library(dplyr)


setwd("~/Documents/GitHub/EcologicalGenomics23/PopGenomics/data") # set the path to where you saved the pcANGSD results on your laptop

## First, let's work on the genetic PCA:

meta <- read.csv("RedSpruce_340_metadata.csv")

COV <- as.matrix(read.table("allGL.cov")) # read in the genetic covariance matrix

PCA <- eigen(COV)

lo <- select(meta, Latitude, Longitude)


library(Hmisc)

res2 <- rcorr(as.matrix(cb))

res2

res3 <- as.numeric(res2)

heatmap(res2)

# returns just the vt vs nc

list.files()

#importing full dataset
bb <- read.csv(file.choose("budbreak_blups_per_garden.csv"), header = FALSE)

names(bb) <- bb[1,]
bb <- bb[-1,]

bb$Budbreak_2020 <- as.numeric(bb$Budbreak_2020)

#removing MD
clbb <- subset(bb, Garden=="Vermont" | Garden=="North_Carolina")

# making a categorical variable
#clbb$Garden <- as.factor(clbb$Garden)

library(tidyr)
library(dplyr)
library(ggpubr)
library(ggplot2)

# getting two seperate blup infos
cbv <- subset(clbb, Garden=="Vermont")
cbn <- subset(clbb, Garden=="North_Carolina")

# renaming column names
names(cbn)[1] <- "NC_Fam"
names(cbn)[2] <- "NC_bb"

names(cbv)[1] <- "VT_Fam"
names(cbv)[2] <- "VT_bb"

c2 <- cbind(cbv, cbn)

cres <- cor.test(c2$VT_bb, c2$NC_bb,
                 method = "pearson")

cres 
# t = 8.4984, df = 332, p-value = 6.595e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
  #0.3303642 0.5070091
#sample estimates:
 # cor = 0.4226931 
#weakly positive association between BLUPs

#checking normality
shapiro.test(c2$VT_bb) # W = 0.99358, p-value = 0.1674
shapiro.test(c2$NC_bb) #W = 0.97841, p-value = 6.526e-05
# not normal should use spearman?


# vt bb
ggqqplot(c2$VT_bb, ylab = "VT BB") # looks normal
ggqqplot(c2$NC_bb, ylab = "NC BB")

cor(c2)

# making the numeric factors

c3 <- c2

c3$VT_Fam <- as.factor(c3$VT_Fam)
c3$NC_Fam <- as.factor(c3$NC_Fam)

library(stringr)

c3$VT_Fam <- str_split(c3$VT_Fam, "\\_", simplify=T)[,1]

# delete other and pool/sort by average family blup score

c4<- subset(c3, select = c(VT_Fam,VT_bb,NC_bb))

c5 <- t(c4)

colnames(c5) <- c5[1,]

c5 <- c5[-1,]

cc <- c5

rc4 <- rcorr(t(c7))

View(cc)

c5 <- as.numeric(c5$VT_bb)

c11 <- t(c10)

colnames(c11) <- c11[1,]

c11 <- c11[-1,]

View(c11)

c12 <- t(c11)

View(c12)

c12 <- unlist(c12)

colnames(c12)[[1]] <- "VT_bb"
colnames(c12)[[2]] <- "NC_bb"

c11[1:65] <- sapply(c11[1:65],as.numeric)

c12 %>% mutate_if(is.character,as.numeric)

c12[[1]] <- as.numeric(c12[[1]])
c12[[2]] <- as.numeric(c12[[2]])

library(dplyr)

c112 <- matrix(as.numeric(c12),    # Convert to numeric matrix
                  ncol = ncol(c12))

c113 <- as.data.frame(c112)

colnames(c113)[1] <- "VT_bb"
colnames(c113)[2] <- "NC_bb"

rownames(c113) <- c10$VT_Fam

c114 <- t(c113)

c115 <- as.matrix(abs(dist(c114['VT_bb', ] - c114['NC_bb', ])))

c117 <- melt(c115)

c113$dif <- c113$dif

c113$dif <- c113$NC_bb - c113$VT_bb

c116<- as.data.frame(as.table(cor(c115)))

cc <- c5

# sorting and summarizing via fam id

c10 <- c5 %>%  
  group_by(VT_Fam) %>%  
  summarise(VT_bb = mean(VT_bb), NC_bb = mean(NC_bb))

c5 <- t(c5)


c5[[VT_bb]] <- as.numeric(c5[[VT_bb]])
c5$NC_bb <- as.numeric(c5$NC_bb)

c6 <- subset(c5, select = "VT_bb","NC_bb")

# convert tibble to data frame

c5 <- data.frame(c5)

ctrans <- t(c5)

c6 <- data.frame(cbind(names(c5), t(c5)))

#changing first row to header c6
c9 <- c11
names(c9) <- c9[1,]

c9 <- c9[,-1]
c10 <- as.matrix(c9)
c10 <- c10[-1,]
cnum1 <- matrix( 
  as.numeric(c10), ncol = 65, dimnames = c10)

c13 <- as.list(c9)
#deleting first and first column row bc not needed
c6 <- data.frame(cbind(names(c5), t(c5)))
names(c6) <- c6[1,]
c6 <- c6[-1,]
c6 <- c6[,-1]

library(tidyverse)
library(reshape2)
#making them all numeric 

c7 <- sapply(c6, as.numeric,USE.NAMES = TRUE)

rownames(c7)<- c("VTb", "NCb")

c8 <- t(c7) #transposing long way
str(c8)

library(ggplot2)
c9 <- as.matrix(dist(c7['VTb', ]))
c99 <- as.matrix(dist(c7['NCb', ]))

c20 <- as.matrix(dist(abs(c7['VTb', ] - c7['NCb', ])))
c40 <- as.matrix(abs(dist(c7['VTb', ] - c7['NCb', ])))
c31 <- as.matrix(abs(diff(c7['NCb', ] - c7['VTb', ])))

c555 <- as.matrix(dist(cc[,'VT_bb']))

# getting the plasticity scores

pp <- as.matrix(abs(dist(c7['VTb', ] - c7['NCb', ])))

p3 <- subset(pp, select = "AB_05")


# convert bb to numeric 

cc <- t(cc)

library(dplyr)

cc2 <- as.data.frame(cc)
cc2 <- matrix(as.numeric(cc), ncol = 2)

rownames(cc2) <- rownames(cc)
colnames(cc2) <- colnames(cc)


# across all families 

dim(cc3)

cc3 <- t(cc2)

cc4 <- as.data.frame(cc3)

cc5 <- t(cc4)

cc5 <- as.data.frame(cc5)


fin1 <- as.matrix(abs(dist(cc5[, 'NC_bb'] - cc5[, 'VT_bb'])))



cfin <- as.matrix(dist(cc2['NCb', ] - cc2['VTb', ]))


c23 <- as.matrix(abs(diff(c7['NCb', ] - c7['VTb', ])))

colnames(c23) <- "Plasticity Score by Locality"

write.table(c23,
            "",
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F)

write.table(c23
c19 <- melt(c9)

colnames(c19) <- c("x", "y", "value")

#showing similarities between families VT
ggplot(c19, aes(x=x, y=y, fill=value))+
  geom_tile() +
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  guides(fill = guide_colourbar(label = FALSE,
                                ticks = FALSE))


c33 <- melt(c31)
c22 <- melt(c20)

colnames(c22) <- c("VT", "NC", "value")

colnames(c33) <- c("x", "y", "value")

c34 <- as.data.frame(c33)

########################
#### NC vs VT BLUP by localities #####
########################

#showing similarities between localities at NC and VT
ggplot(c22, aes(x=VT, y=NC, fill=value))+
  geom_tile() +
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  guides(fill = guide_colourbar(label = TRUE,
                                ticks = TRUE))+
  labs(fill = "Plasticity Scores", title = "Early Budbreak Plasticity by Locality" )

fin2 <- melt(fin1)

colnames(fin2) <- c("VT", "NC", "value")

# heat map of families 
ggplot(fin2, aes(x=VT, y=NC, fill=value))+
  geom_tile() +
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  guides(fill = guide_colourbar(label = TRUE,
                                ticks = TRUE))+
  labs(fill = "Plasticity Scores", title = "Early Budbreak Plasticity by Individual")



ggplot(c34, aes(x=x, y=y, fill=value))+
  geom_tile() +
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  guides(fill = guide_colourbar(label = FALSE,
                                ticks = FALSE))


heatmap(c9, Colv = NA, Rowv = NA)
heatmap(c99)

library(Hmisc)


# how to do placisity measurement scatter plot

# import the ancestry scores (these are the .Q files)

q <- read.table("allGL.admix.3.Q", sep=" ", header=F)

colnames(q) <- c("V1", "V2", "V3")

qq = subset(q, select = -c(MDAd))

cb$VT <- as.numeric(cb$VT)
cb$NC <- as.numeric(cb$NC)

fin <- cbind(q, cb)

d1 <- cb

d2 <- d1

fin$dif <- fin$dif

fin$dif <- (abs(fin$NC - fin$VT))

d2$dif <- abs(d2$VT - d2$NC)

din <- cbind()

ggplot(fin, aes(V1, dif)) +
         geom_point()

ggplot(fin, aes(V2, dif)) +
  geom_point()

ggplot(fin, aes(V3, dif)) +
  geom_point()

fin2 <- data.frame(Admixture = c(fin$V1, fin$V2, fin$V3), Location = c(rep(c("VT"), 333), rep(c("NC"), 333)), BLUP = fin$dif))


fin2 <- data.frame(Admixture = c(fin$V1, fin$V2, fin$V3), Source = c(rep(c("Core"), 333), rep(c("Margin"), 333), rep(c("Edge"), 333)), BLUP = fin$dif)

                   
scat <- ggplot(fin2, aes(x=Admixture,y=BLUP, color=Source)) +
         geom_point()

scat +labs(x = "Admixture Score", y = "Plasticity Score", title = "Admixture vs. Budbreak Plasticity of Red Spruce", color= "Source Population")


# importing latitude and longitude

library(ggplot2) # plotting
library(ggpubr) # plotting
library(dplyr)


setwd("~/Documents/GitHub/EcologicalGenomics23/PopGenomics/data") # set the path to where you saved the pcANGSD results on your laptop

## First, let's work on the genetic PCA:

meta <- read.csv("RedSpruce_340_metadata.csv")

COV <- as.matrix(read.table("allGL.cov")) # read in the genetic covariance matrix

PCA <- eigen(COV)

lo <- dplyr::select(meta, Latitude, Longitude, Family)

total <- merge(lo, bb, by = "Family")

t2 <- subset(total, Garden=="Vermont" | Garden=="North_Carolina")

t3 <- t2 %>%
  group_by(Family) %>%
  summarise(dif = abs(first(Budbreak_2020) - last(Budbreak_2020)), Latitude = first(Latitude), Longitude = first(Longitude))

library(wesanderson)
library(viridis)

ggplot(t3, aes(x=Latitude, y=Longitude, color=dif)) +
  geom_point(alpha=0.5, size=3) +
  scale_color_viridis(discrete = FALSE)

library(colors3d)
library(raster)

lol <- colors2d(lo[, c("Latitude", "Longitude")])

load <- dplyr::select(meta, Latitude, Longitude, Family, Elevation, Pop)

law <- load %>% slice(1:334)

p4 <- cbind(law, p3)

colnames(p4)[6] <- "Plasticity Score"

p5 <- subset(p4, select = c("Latitude", "Longitude", "Elevation", "Plasticity Score"))

p4$Elevation <- as.numeric(p4$Elevation)

p6 <- melt(p4)

ggplot(p4, aes(x=Latitude, y=`Plasticity Score`, color=Pop)) +
  geom_point()+
  labs(title = "Plasticity Scores by Latitude", color="Locality")
  

ggplot(p4, aes(x=Elevation, y=`Plasticity Score`, color=Pop)) +
  geom_point()+
  labs(title = "Plasticity Scores by Elevation", color="Locality")

library(RColorBrewer)
# lat by long 
ggplot(p4, aes(x=Latitude, y=Longitude, color=`Plasticity Score`)) +
  geom_point(alpha=0.5)+
  scale_color_gradient(low="blue",high="red")+
  labs(title = "Plasticity Score by Location", color="Plasticity Score")

ggplot(t3, aes(x=Latitude, y=dif)) +
  geom_point(alpha=)

plot(p4$Latitude, p4$`Plasticity Score`)

#Create a custom color scale
library(RColorBrewer)
myColors <- brewer.pal(5,"Set1")
names(myColors) <- levels(p4$Family)
colScale <- scale_colour_manual(name = "Family",values = myColors)

tr <- ggplot(p4, aes(Latitude, `Plasticity Score`, color=Family)) +
  geom_point()

