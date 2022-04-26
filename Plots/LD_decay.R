library(dplyr)
library(stringr)
library(ggplot2)
library(DataCombine)

##### Using tomahawk 
# https://www.biostars.org/p/347796/

# define plotting functions
#' @title plotPairwiseLD
#' @description Plots R2 heatmap across the chromosome (like Haploview)
#' @param dfr A data.frame with minimum CHROM_A, POS_A, CHROM_B, POS_B and R2.
#' An output from tomahawk works.
#' @param chr A chromosome name.
#' @param xlim A two number vector specifying min and max x-axis limits. Any one or both can be defaulted by specifying NA.
#' @param ylim A two number vector specifying min and max y-axis limits. Any one or both can be defaulted by specifying NA.
#' @param minr2 A value between 0 and 1. All SNPs with R2 value below this 
#' value is excluded from plot.
#' 
plotPairwiseLD <- function(dfr,chr,xlim=c(NA,NA),ylim=c(NA,NA),minr2) {
  if(missing(dfr)) stop("Input data.frame 'dfr' missing.")
  
  if(!missing(chr)) {
    ld <- filter(ld,CHROM_A==get("chr") & CHROM_B==get("chr"))
  }
  ld <- filter(ld,POS_A<POS_B)
  
  if(!missing(minr2)) {
    ld <- filter(ld,R2>get("minr2"))
  }
  
  ld <- ld %>% arrange(R2)
  
  ld$x <- ld$POS_A+((ld$POS_B-ld$POS_A)/2)
  ld$y <- ld$POS_B-ld$POS_A
  ld$r2c <- cut(ld$R2,breaks=seq(0,1,0.2),labels=c("0-0 - 0.2","0.2 - 0.4",
                                                   "0.4 - 0.6","0.6 - 0.8",
                                                   "0.8 - 1.0"))
  ld$r2c <- factor(ld$r2c,levels=rev(c("0-0 - 0.2","0.2 - 0.4",
                                       "0.4 - 0.6","0.6 - 0.8",
                                       "0.8 - 1.0")))
  
  ggplot(ld,aes(x=x,y=y,col=r2c))+
    geom_point(shape=20,size=0.1,alpha=0.8)+
    scale_color_manual(values=c("#ca0020","#f4a582","#d1e5f0","#67a9cf","#2166ac"))+
    scale_x_continuous(limits=xlim)+
    scale_y_continuous(limits=ylim)+
    guides(colour=guide_legend(title="R2",override.aes=list(shape=20,size=8)))+
    labs(x="Chromosome (Bases)",y="")+
    theme_bw(base_size=14)+
    theme(panel.border=element_blank(),
          axis.ticks=element_blank()) %>%
    return()
}

#' @title plotDecayLD
#' @description Plots R2 heatmap across the chromosome (like Haploview)
#' @param dfr A data.frame with minimum CHROM_A, POS_A, CHROM_B, POS_B and R2.
#' An output from tomahawk works.
#' @param chr A chromosome name.
#' @param xlim A two number vector specifying min and max x-axis limits. Any one or both can be defaulted by specifying NA.
#' @param ylim A two number vector specifying min and max y-axis limits.  Any one or both can be defaulted by specifying NA.
#' @param avgwin An integer specifying window size. Mean R2 is computed within windows.
#' @param minr2 A value between 0 and 1. All SNPs with R2 value below this 
#' value is excluded from plot.
#' 
plotDecayLD <- function(dfr,chr,xlim=c(NA,NA),ylim=c(NA,NA),avgwin=0,minr2) {
  if(missing(dfr)) stop("Input data.frame 'dfr' missing.")
  
  if(!missing(chr)) {
    ld <- filter(ld,CHROM_A==get("chr") & CHROM_B==get("chr"))
  }
  ld <- filter(ld,POS_A<POS_B)
  
  if(!missing(minr2)) {
    ld <- filter(ld,R2>get("minr2"))
  }
  
  ld <- ld %>% arrange(R2)
  
  ld$dist <- ld$POS_B-ld$POS_A
  
  if(avgwin>0) {
    ld$distc <- cut(ld$dist,breaks=seq(from=min(ld$dist),to=max(ld$dist),by=avgwin))
    ld <- ld %>% group_by(distc) %>% summarise(dist=mean(dist),R2=mean(R2))
  }
  
  ggplot(ld,aes(x=dist,y=R2))+
    geom_point(shape=20,size=0.15,alpha=0.7)+
    geom_smooth()+
    scale_x_continuous(limits=xlim)+
    scale_y_continuous(limits=ylim)+
    labs(x="Distance (Bases)",y=expression(LD~(r^{2})))+
    theme_bw(base_size=14)+
    theme(panel.border=element_blank(),
          axis.ticks=element_blank()) %>%
    return()
}

# Then we read the file.
ld <- read.delim("~/Dropbox/paperultima/LD decay/pel.chr1.ld",sep="\t",comment.char="#")
colnames(ld) <- c("FLAG","CHROM_A","POS_A","CHROM_B","POS_B","REF_REF","REF_ALT","ALT_REF","ALT_ALT","D","Dprime","R","R2","P","ChiSqMode","ChiSqTable")

plotPairwiseLD(ld) #Plot pairwise LD plot.

#Fig: LD decay plot with a min R2 limit of 0.25 and R2 values averaged within 10 Kb windows.
plotDecayLD(ld,minr2=0.25,avgwin=100000)

# Simulated population
ld2 <- read.delim("~/Dropbox/paperultima/LD decay/bottleneck125i50.ld",sep="\t",comment.char="#")
colnames(ld2) <- c("FLAG","CHROM_A","POS_A","CHROM_B","POS_B","REF_REF","REF_ALT","ALT_REF","ALT_ALT","D","Dprime","R","R2","P","ChiSqMode","ChiSqTable")

plotDecayLD(ld2,minr2=0.25,avgwin=100000)







###============================================================================
######### Alternatively, more manual way
# See https://www.biostars.org/p/300381/
## Look below for an alternative (from https://www.biostars.org/p/347796/)

dfr <- read.delim("~/Dropbox/paperultima/BR.summary",sep="",header=,check.names=F,stringsAsFactors=F)
ldd <- read.delim("/Users/Pedro/Dropbox/input_chr16/LDfigures/Lddecay.BR.stat",sep="\t",
                  header=T,check.names=F,stringsAsFactors=F)
colnames(ldd) <- c("dist","rsq", "mean_D", "sum", "sumD", "pairs")

ldd$distc <- cut(ldd$dist,breaks=seq(from=min(ldd$dist)-1,to=max(ldd$dist)+1,by=10000))

ldd1 <- ldd %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))

ldd1 <- ldd1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))

ggplot()+
  geom_point(data=ldd1,aes(x=start,y=mean),size=0.4,colour="grey20")+
  geom_line(data=ldd1,aes(x=start,y=mean),size=0.3,alpha=0.5,colour="grey40")+
  labs(x="Distance (Megabases)",y=expression(LD~(r^{2})))+
  #  scale_x_continuous(breaks=c(0,2*10^6,4*10^6,6*10^6,8*10^6),labels=c("0","2","4","6","8"))+
  theme_bw()



