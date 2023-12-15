#Ashley Lantigua/JH
#12/06/23
#eigen value 1

library(ggplot2) # plotting
library(ggpubr) # plotting
library(dplyr)


setwd("~/Documents/GitHub/EcologicalGenomics23/PopGenomics/data") # set the path to where you saved the pcANGSD results on your laptop

## First, let's work on the genetic PCA:

meta <- read.csv("RedSpruce_340_metadata.csv")

COV <- as.matrix(read.table("allGL.cov")) # read in the genetic covariance matrix

PCA <- eigen(COV)

lo <- select(meta, Latitude, Longitude)

data=as.data.frame(PCA$vectors)
data=data[,c(1:2)] # the second number here is the number of PC axes you want to keep

write.table(data,
            "allGL_genPC1_2.txt",
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F)

# extract the principal components from the COV matrix

## How much variance is explained by the first few PCs?

var <- round(PCA$values/sum(PCA$values),3)

var[1:3] #0.025 0.023 0.003

# A "screeplot" of the eigenvalues of the PCA:

barplot(var, 
        xlab="Eigenvalues of the PCA", 
        ylab="Proportion of variance explained",
        main="Scree Plot of Budbreak in Red Spruce")

## Bring in the bam.list file and extract the sample info:

names <- read.table("G2_budbreak_bam.list")
names <- unlist(strsplit(basename(as.character(names[,1])), split = ".sorted.rmdup.bam"))

split = strsplit(names, "_")
pops <- data.frame(names[1:333], do.call(rbind, split[1:333]))
names(pops) = c("Ind", "Pop", "Row")

## A quick and humble PCA plot:

plot(PCA$vectors[,1:2],
     col=as.factor(pops[,2]),
     xlab="PC1",ylab="PC2", 
     main="Genetic PCA")

## A more beautiful PCA plot using ggplot :)

data=as.data.frame(PCA$vectors)
data=data[,c(1:3)]
data= cbind(data, pops)

data2 <- unlist(strsplit(data$Ind,".", fixed = TRUE))
data2 <- data2[!data2 == "bam"]
data2 <- data2[!data2 == "final"]


str(data2)


data3 <- cbind(data, data2)

str(data3)

names(data3)[names(data3) == "data2"] <- "Family"

str(data3)

str(meta)
finaldata <- dplyr::left_join(data3, meta, by = "Family") 

finald <- dplyr::select(finaldata, Pop.x, Family, Region, Row, V1, V2)




cols=c("#377eB8","#EE9B00","#0A9396")

ggscatter(finald, x = "V1", y = "V2",
          color = "Region",
          mean.point = TRUE,
          star.plot = TRUE) +
  theme_bw(base_size = 13, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.7)), 
        axis.text = element_text(size=13), 
        legend.position = "right") +
  labs(x = paste0("PC1: (",var[1]*100,"%)"), y = paste0("PC2: (",var[2]*100,"%)")) +
  scale_color_manual(values=c(cols), name="Source population", labels = c("Core", "Margin", "Edge")) +
  guides(colour = guide_legend(nrow = 2))


h## Next, we can look at the admixture clustering:

# import the ancestry scores (these are the .Q files)


q <- read.table("allGL.admix.3.Q", sep=" ", header=F)

str(q)


K = dim(q)[2] #Find the level of K modeled

# k = 3
## order according to population code

ord<-order(pops[,2])

ord

# make the plot:
barplot(t(q)[,ord],
        col=cols[1:K],
        space=0,border=NA,
        xlab="Populations",ylab="Admixture proportions",
        main=paste0("Red spruce K=",K))
text(tapply(1:nrow(pops),pops[ord,2],mean),-0.05,unique(pops[ord,2]),xpd=T) 
abline(v=cumsum(sapply(unique(pops[ord,2]),function(x){sum(pops[ord,2]==x)})),col=1,lwd=1.2)

list.files()
