# Ecogen GWAS Manhattan Plot
#12/10.23


###################################
#  Manhattan 4 Red Spruce #
###################################

library(RcppCNPy) # for reading python numpy (.npy) files

setwd("~/Documents/GitHub/EcologicalGenomics23/PopGenomics/data")

list.files()

### read in selection statistics (these are chi^2 distributed)

#commented out to upload s<-npyLoad("allGL_poly.selection.npy")

# convert test statistic to p-value
pval <- as.data.frame(1-pchisq(s,1))
names(pval) = "p_PC1"

## read positions
# commented out to upload p <- read.table("allGL_poly_mafs.sites",sep="\t",header=T, stringsAsFactors=T)
dim(p)
# [1] 1923974       8


p_filtered = p[which(p$kept_sites==1),]
dim(p_filtered)
# [1] 1393931       8

# How many sites got filtered out when testing for selection? Why?

## make manhattan plot
mh <- manhattanr(plot(-log10(pval$p_PC1),
     col=p_filtered$chromo,
     xlab="CHR",
     ylab="-log10(p-value)",
     main="Manhattan Plot of Budbreak in Red Spruce (K3)"))

# We can zoom in if there's something interesting near a position...

plot(-log10(pval$p_PC1[2e05:2.01e05]),
     col=p_filtered$chromo, 
     xlab="Position", 
     ylab="-log10(p-value)", 
     main="Selection outliers: pcANGSD e=1 (K2)")

#other manhattan plot attempt

library(qqman)

manhattan(pval$p_PC1, p_filtered$chromo)

library(lattice)

manhattan.plot<-function(chr, pos, pvalue, 
                         sig.level=NA, annotate=NULL, ann.default=list(),
                         should.thin=T, thin.pos.places=2, thin.logp.places=2, 
                         xlab="Chromosome", ylab=expression(-log[10](p-value)),
                         col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {
        
        if (length(chr)==0) stop("chromosome vector is empty")
        if (length(pos)==0) stop("position vector is empty")
        if (length(pvalue)==0) stop("pvalue vector is empty")
        
        #make sure we have an ordered factor
        if(!is.ordered(chr)) {
                chr <- ordered(chr)
        } else {
                chr <- chr[,drop=T]
        }
        
        #make sure positions are in kbp
        if (any(pos>1e6)) pos<-pos/1e6;
        
        #calculate absolute genomic position
        #from relative chromosomal positions
        posmin <- tapply(pos,chr, min);
        posmax <- tapply(pos,chr, max);
        posshift <- head(c(0,cumsum(posmax)),-1);
        names(posshift) <- levels(chr)
        genpos <- pos + posshift[chr];
        getGenPos<-function(cchr, cpos) {
                p<-posshift[as.character(cchr)]+cpos
                return(p)
        }
        
        #parse annotations
        grp <- NULL
        ann.settings <- list()
        label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
                            col=NULL, fontface=NULL, fontsize=NULL, show=F)
        parse.label<-function(rawval, groupname) {
                r<-list(text=groupname)
                if(is.logical(rawval)) {
                        if(!rawval) {r$show <- F}
                } else if (is.character(rawval) || is.expression(rawval)) {
                        if(nchar(rawval)>=1) {
                                r$text <- rawval
                        }
                } else if (is.list(rawval)) {
                        r <- modifyList(r, rawval)
                }
                return(r)
        }
        
        if(!is.null(annotate)) {
                if (is.list(annotate)) {
                        grp <- annotate[[1]]
                } else {
                        grp <- annotate
                } 
                if (!is.factor(grp)) {
                        grp <- factor(grp)
                }
        } else {
                grp <- factor(rep(1, times=length(pvalue)))
        }
        
        ann.settings<-vector("list", length(levels(grp)))
        ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
        
        if (length(ann.settings)>1) { 
                lcols<-trellis.par.get("superpose.symbol")$col 
                lfills<-trellis.par.get("superpose.symbol")$fill
                for(i in 2:length(levels(grp))) {
                        ann.settings[[i]]<-list(pch=pch, 
                                                col=lcols[(i-2) %% length(lcols) +1 ], 
                                                fill=lfills[(i-2) %% length(lfills) +1 ], 
                                                cex=cex, label=label.default);
                        ann.settings[[i]]$label$show <- T
                }
                names(ann.settings)<-levels(grp)
        }
        for(i in 1:length(ann.settings)) {
                if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
                ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                                                      parse.label(ann.settings[[i]]$label, levels(grp)[i]))
        }
        if(is.list(annotate) && length(annotate)>1) {
                user.cols <- 2:length(annotate)
                ann.cols <- c()
                if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
                        ann.cols<-match(names(annotate)[-1], names(ann.settings))
                } else {
                        ann.cols<-user.cols-1
                }
                for(i in seq_along(user.cols)) {
                        if(!is.null(annotate[[user.cols[i]]]$label)) {
                                annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
                                                                            levels(grp)[ann.cols[i]])
                        }
                        ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
                                                                annotate[[user.cols[i]]])
                }
        }
        rm(annotate)
        
        #reduce number of points plotted
        if(should.thin) {
                thinned <- unique(data.frame(
                        logp=round(-log10(pvalue),thin.logp.places), 
                        pos=round(genpos,thin.pos.places), 
                        chr=chr,
                        grp=grp)
                )
                logp <- thinned$logp
                genpos <- thinned$pos
                chr <- thinned$chr
                grp <- thinned$grp
                rm(thinned)
        } else {
                logp <- -log10(pvalue)
        }
        rm(pos, pvalue)
        gc()
        
        #custom axis to print chromosome names
        axis.chr <- function(side,...) {
                if(side=="bottom") {
                        panel.axis(side=side, outside=T,
                                   at=((posmax+posmin)/2+posshift),
                                   labels=levels(chr), 
                                   ticks=F, rot=0,
                                   check.overlap=F
                        )
                } else if (side=="top" || side=="right") {
                        panel.axis(side=side, draw.labels=F, ticks=F);
                }
                else {
                        axis.default(side=side,...);
                }
        }
        
        #make sure the y-lim covers the range (plus a bit more to look nice)
        prepanel.chr<-function(x,y,...) { 
                A<-list();
                maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+.5;
                A$ylim=c(0,maxy);
                A;
        }
        
        xyplot(logp~genpos, chr=chr, groups=grp,
               axis=axis.chr, ann.settings=ann.settings, 
               prepanel=prepanel.chr, scales=list(axs="i"),
               panel=function(x, y, ..., getgenpos) {
                       if(!is.na(sig.level)) {
                               #add significance line (if requested)
                               panel.abline(h=-log10(sig.level), lty=2);
                       }
                       panel.superpose(x, y, ..., getgenpos=getgenpos);
                       if(!is.null(panel.extra)) {
                               panel.extra(x,y, getgenpos, ...)
                       }
               },
               panel.groups = function(x,y,..., subscripts, group.number) {
                       A<-list(...)
                       #allow for different annotation settings
                       gs <- ann.settings[[group.number]]
                       A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
                       A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
                       A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
                       A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
                       A$x <- x
                       A$y <- y
                       do.call("panel.xyplot", A)
                       #draw labels (if requested)
                       if(gs$label$show) {
                               gt<-gs$label
                               names(gt)[which(names(gt)=="text")]<-"labels"
                               gt$show<-NULL
                               if(is.character(gt$x) | is.character(gt$y)) {
                                       peak = which.max(y)
                                       center = mean(range(x))
                                       if (is.character(gt$x)) {
                                               if(gt$x=="peak") {gt$x<-x[peak]}
                                               if(gt$x=="center") {gt$x<-center}
                                       }
                                       if (is.character(gt$y)) {
                                               if(gt$y=="peak") {gt$y<-y[peak]}
                                       }
                               }
                               if(is.list(gt$x)) {
                                       gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
                               }
                               do.call("panel.text", gt)
                       }
               },
               xlab=xlab, ylab=ylab, 
               panel.extra=panel.extra, getgenpos=getGenPos, ...
        );
}

pvv <- ( -log10(pval$p_PC1[2e05:2.01e05]))
pf2 <- merge(p_filtered, pvv)

manhattan.plot(p_filtered$chromo, p_filtered$position, -log10(pval$p_PC1[2e05:2.01e05]))

# get the contig with the lowest p-value for selection
sel_contig <- p_filtered[which(pval==min(pval$p_PC1)),c("chromo","position")]
sel_contig

#           chromo  position lowest p-value
#529962 MA_856442     2564

write.table(unique(sel_contig$chromo),
            "topcon.txt", 
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F)

topcontig <- (p_filtered[which(pval==min(pval$p_PC1)),c("chromo","position")])



# get all the outliers with p-values below some cutoff
cutoff=1e-3   # equals a 1 in 5,000 probability
outlier_contigs <- p_filtered[which(pval<cutoff),c("chromo","position")]
outlier_contigs

# get all the outliers with p-values below some cutoff
cutoff6=1e-6   # equals a 1 in 5,000 probability
outlier_contigs6 <- p_filtered[which(pval<cutoff6),c("chromo","position")]
outlier_contigs6

# get all the outliers with p-values below some cutoff
cutoff7=1e-7   # equals a 1 in 5,000 probability
outlier_contigs7 <- p_filtered[which(pval<cutoff7),c("chromo","position")]
outlier_contigs7

write.table(unique(outlier_contigs7$chromo),
            "top50contigs.txt", 
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F)

dim(outlier_contigs7)[1]

cutoff8=1e-8   # equals a 1 in 5,000 probability
outlier_contigs8 <- p_filtered[which(pval<cutoff8),c("chromo","position")]
outlier_contigs8

# top 6 most relevant contigs 
#   chromo            positon
#529961  MA_856442     2550
#529962  MA_856442     2564
#1385462  MA_17819    86835
#1448203   MA_3745    11224
#1523923  MA_62267    12494
#1920306  MA_19598     1209


write.table(unique(outlier_contigs8$chromo),
            "top6contigs.txt", 
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F)



# how many outlier loci < the cutoff?
dim(outlier_contigs)[1]

# how many unique contigs harbor outlier loci?
length(unique(outlier_contigs$chromo))

write.table(unique(outlier_contigs$chromo),
            "allGL_poly_PC1_outlier_contigs.txt", 
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F)

### Day 8 Code for outliers

outliers_PC1 <- p_filtered[which(pval$p_PC1<cutoff),c("chromo","position")]

# how many outlier loci < the cutoff?
dim(outliers_PC1)[1]
# [1] 9260

# write them out to a file
write.table(outliers_PC1,
            "allGL_poly_outliers_PC1.txt", 
            sep=":",
            quote=F,
            row.names=F,
            col.names=F)

COV <- as.matrix(read.table("allGL_poly.cov"))

PCA <- eigen(COV)

data=as.data.frame(PCA$vectors)
data=data[,c(1:2)] # the second number here is the number of PC axes you want to keep

write.table(data,
            "allGL_poly_genPC1_2.txt",
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F)


library(manhattanly)

mh <- manhattanr(x = p_filtered, chr = p_filtered$chromo, bp = p_filtered$position, pval$p_PC1)

manhattanly(p_filtered, snp = "position", )
