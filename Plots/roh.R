# 21 April 2019
# Runs of homozygosity
########################################

setwd("~/Dropbox/paperultima/ROH/")

br.indiv <- read.table("BR_n171.thin6.1Mb.ROHs.hom.indiv", header=TRUE, stringsAsFactors=FALSE)

## plot(NSEG~KB, data=br.indiv)

# mean density: mean(br.hom$DENSITY) # 17.82682 kb and others: mean(pops.hom$DENSITY[!pops.hom$pop=="BR"])
br.hom <- read.table("BR_n171.thin6.1Mb.ROHs.hom", header=TRUE, stringsAsFactors=FALSE)

# br.hom <- read.table("BR_n171.thin6.ROHs.hom", header=TRUE, stringsAsFactors=FALSE) # Galina's

## boxplot(br.hom$KB)
## boxplot(br.hom[br.hom$KB>450 & br.hom$KB<1500, ]$KB)


## Populations
## read *.hom data for all the populations

## poporder <- c(FIN, CEU, GBR, IBS, TSI, NAM, PEL, MXL, CLM, PUR, BR, ASW, ACB, LWK, ESN, YRI, MSL, GWD)

popnames <- c("GBR", "FIN", "PUR", "CLM", "IBS", "PEL", "ACB", "GWD", "ESN", "MSL", "CEU", "YRI", "LWK", "ASW", "MXL", "TSI", "NAM")

pops.hom <- cbind( pop="BR", br.hom )

## Merge all populations into one table
 
for(i in 1: length(popnames)) {    
    myfile <- paste(popnames[i], ".thin6.1Mb.ROHs.hom", sep="" )
    tmp <- read.table( myfile, header=TRUE, stringsAsFactors=FALSE )
    tmp <- cbind( pop = popnames[i], tmp )
    pops.hom <- rbind( pops.hom, tmp )
}

boxplot(KB ~ pop, pops.hom) 
boxplot(KB[KB>450 & KB<1500] ~ pop[KB>450 & KB<1500], pops.hom)
boxplot(KB[KB<1500] ~ pop[KB<1500], pops.hom)
boxplot(KB[KB<5000] ~ pop[KB<5000], pops.hom)
boxplot(KB[KB<15000] ~ pop[KB<15000], pops.hom, ylab="Kb") # I think this can go to supplements

          
## Add an extra column with convertion of IID from character to factor
pops.hom$IIDasfactor <- as.factor(pops.hom$IID)

## Add an extra column with numeric ID to each individual

pops.hom$IIDasnumeric <- as.numeric(pops.hom$IIDasfactor)


## read *.indiv data for all the populations

pops.indiv <- cbind( pop="BR", br.indiv )

## merge all populations together
for(i in 1: length(popnames)) {    
    myfile <- paste( popnames[i], ".thin6.1Mb.ROHs.hom.indiv", sep="" )
    tmp <- read.table( myfile, header=TRUE )
    tmp <- cbind( pop = popnames[i], tmp )
    pops.indiv <- rbind( pops.indiv, tmp )
}

## Averages (kb)                                          Median
mean(pops.hom$KB[pops.hom$pop=="BR"]) #  [1] 1798.214   1293.272 ± 2014.296
mean(pops.hom$KB[pops.hom$pop=="PUR"]) # [1] 1931.319   1320.979 ± 2420.271
mean(pops.hom$KB[pops.hom$pop=="CLM"]) # [1] 2157.274   1327.207 ± 3039.986 
mean(pops.hom$KB[pops.hom$pop=="MXL"]) # [1] 1601.101   1277.03  ± 1573.584
mean(pops.hom$KB[pops.hom$pop=="PEL"]) # [1] 1578.037   1319.637 ± 1321.919
mean(pops.hom$KB[pops.hom$pop=="NAM"]) # [1] 1690.565   1334.639 ± 1622.895

## Test: lets plot some populations
par(tck = -.01)
plot( NSEG~KB, data=pops.indiv[pops.indiv$pop=="BR", ], col="green", pch=20, ylim=c(1,150) )
points(NSEG~KB, data=pops.indiv[pops.indiv$pop=="PEL", ], col="red", pch=20, ylim=c(1,150))
points(NSEG~KB, data=pops.indiv[pops.indiv$pop=="CLM", ], col="grey", pch=20, ylim=c(1,150))
# points(NSEG~KB, data=pops.indiv[pops.indiv$pop=="CEU", ], col="darkgreen", pch=20)
# points(NSEG~KB, data=pops.indiv[pops.indiv$pop=="GBR", ], col="pink", pch=20)

points(NSEG~KB, data=pops.indiv[pops.indiv$pop=="GWD", ], col="darkgreen", pch=1, cex=0.8)

##############################
## I will use color palette from "pophelper" package
standard_12=c("#2121D9", "#9999FF", "#DF0101","#04B404", "#FFFB23", "#FF9326", "#A945FF", "#0089B2", "#B26314", "#610B5E", "#FE2E9A", "#BFF217")

## Add column switching KB to Mb values

pops.indiv$Mb <- pops.indiv$KB/1000

## Load the required packages
library(Hmisc)
library(FactoMineR)

##################################################

## Plot European populations
tiff(filename = "ROH_NSEG_EUR.tiff", width = 200, height = 200, units = "mm",
     compression = "lzw", res = 1500, type = "cairo", antialias = "gray")

par(tck = -.01)
plot(NSEG~Mb, data=pops.indiv[pops.indiv$pop=="FIN", ], col="blue", pch=20, ylim=c(70, 150), xlim=c(50, 155))
points(NSEG~Mb, data=pops.indiv[pops.indiv$pop=="TSI", ], col="red", pch=20)
points(NSEG~Mb, data=pops.indiv[pops.indiv$pop=="IBS", ], col="pink", pch=20)
points(NSEG~Mb, data=pops.indiv[pops.indiv$pop=="GBR", ], col="green", pch=20)
title("European populations")
legend(50,150, legend=c("FIN","TSI","IBS","GBR"), col=c("blue","red","pink","green"), pch=20, pt.cex=2)

dev.off()

############
## Plot Latin Americans populations

tiff(filename = "ROH_NSEG_LatAm.tiff", width = 200, height = 200, units = "mm",
     compression = "lzw", res = 1500, type = "cairo", antialias = "gray")

par(tck = -.01)

## Brazilian: BR

#latin <-  pops.indiv[pops.indiv$pop==c("PUR","BR","CLM"), c(1,5,8)]
# with(latin,  dataEllipse(Mb, NSEG,
#                               groups = latin$pop,
#                               id=list(n=3, labels=rownames(latin)),
#                               ylim=c(1,150), xlim=c(1,300), 
#                               center.pch="+", level=.50, fill=TRUE, 
#                               fill.alpha=0.1, add = T))

plot( NSEG~Mb, data=pops.indiv[pops.indiv$pop=="BR", ], 
      col=adjustcolor(standard_12[3], alpha.f=1), pch=".", 
      cex=2, ylim=c(1,150), xlim=c(1,300), xlab="Size (Mb)" ,ylab="Number of segments" )
dataEllipse(pops.indiv$Mb[pops.indiv$pop=="BR" ], pops.indiv$NSEG[pops.indiv$pop=="BR"], 
            ylim=c(1,150), xlim=c(1,300),level=.5, col=standard_12[3],  
            center.pch=FALSE, 
            fill=TRUE, fill.alpha=0.1,plot.points=F, add = T)
text( 48, 8, "BR", col=standard_12[3])

# br <- pops.indiv[pops.indiv$pop=="BR", c(1, 5, 8)]
# br.2 <- coord.ellipse(br, level.conf = 0.5, npoint = 700, bary=FALSE)
# points(br.2$res[, 2:3], pch=".", cex=1, col=standard_12[3])
# polygon(br.2$res[, 2:3], border=NA, col=adjustcolor(standard_12[3], alpha.f=0.2))
title("Latin Americans")

## Colombians: CLM
points( NSEG~Mb, data=pops.indiv[pops.indiv$pop=="CLM", ], 
        col=adjustcolor(standard_12[2], alpha.f=1), pch=".", 
        cex=2, ylim=c(1,150), xlim=c(1,300), xlab="Size (Mb)",ylab="Number of segments"  )
dataEllipse(pops.indiv$Mb[pops.indiv$pop=="CLM" ], pops.indiv$NSEG[pops.indiv$pop=="CLM"], 
            ylim=c(1,150), xlim=c(1,300),level=.5, col=standard_12[2],  
            center.pch=FALSE, ellipse.label ="CLM", 
            fill=TRUE, fill.alpha=0.1,plot.points=F, add = T)
# clm <- pops.indiv[pops.indiv$pop=="CLM", c(1, 5, 8)]
# clm.2 <- coord.ellipse(clm, level.conf = 0.5, npoint = 700, bary=FALSE)
# points(clm.2$res[, 2:3], pch=".", cex=1, col=standard_12[2])
# polygon(clm.2$res[, 2:3], border=NA, col=adjustcolor(standard_12[2], alpha.f=0.2))
text( 85, 135, "CLM", col=standard_12[2], font=2)

## Native Americans: NAM
points( NSEG~Mb, data=pops.indiv[pops.indiv$pop=="NAM", ], 
        col=adjustcolor(standard_12[6], alpha.f=1), pch=".", 
        cex=2, ylim=c(1,150), xlim=c(1,300), xlab="Size (Mb)",ylab="Number of segments"  )
dataEllipse(pops.indiv$Mb[pops.indiv$pop=="NAM" ], pops.indiv$NSEG[pops.indiv$pop=="NAM"], 
            ylim=c(1,150), xlim=c(1,300),level=.5, col=standard_12[6],  
            center.pch=FALSE, ellipse.label ="NAM",
            fill=TRUE, fill.alpha=0.1,plot.points=F, add = T)

# nam <- pops.indiv[pops.indiv$pop=="NAM", c(1, 5, 8)]
# nam.2 <- coord.ellipse(nam, level.conf = 0.5, npoint = 2700, bary=FALSE)
# points(nam.2$res[, 2:3], pch=".", cex=1, col=standard_12[6])
# polygon(nam.2$res[, 2:3], border=NA, col=adjustcolor(standard_12[6], alpha.f=0.2))
text( 290, 303, "NAM", col=standard_12[6], font=2)

## Mexicans: MXL
points( NSEG~Mb, data=pops.indiv[pops.indiv$pop=="MXL", ], 
        col=adjustcolor(standard_12[4], alpha.f=1), pch=".",
        cex=2, ylim=c(1,150), xlim=c(1,300), xlab="Size (Mb)",ylab="Number of segments"  )
dataEllipse(pops.indiv$Mb[pops.indiv$pop=="MXL" ], pops.indiv$NSEG[pops.indiv$pop=="MXL"], 
            ylim=c(1,150), xlim=c(1,300),level=.5, col=standard_12[4], 
            center.pch=FALSE,
            fill=TRUE, fill.alpha=0.1,plot.points=F, add = T)
text( 30, 50, "MXL", col=standard_12[4])
# mxl <- pops.indiv[pops.indiv$pop=="MXL", c(1, 5, 8)]
# mxl.2 <- coord.ellipse(mxl, level.conf = 0.5, npoint = 700, bary=FALSE)
# points(mxl.2$res[, 2:3], pch=".", cex=1, col=standard_12[4])
# polygon(mxl.2$res[, 2:3], border=NA, col=adjustcolor(standard_12[4], alpha.f=0.2))

## Peruvians: PEL
points( NSEG~Mb, data=pops.indiv[pops.indiv$pop=="PEL", ], 
        col=adjustcolor(standard_12[1], alpha.f=1), pch=".", 
        cex=2, ylim=c(1,150), xlim=c(1,300), xlab="Size (Mb)",ylab="Number of segments"  )
dataEllipse(pops.indiv$Mb[pops.indiv$pop=="PEL" ], pops.indiv$NSEG[pops.indiv$pop=="PEL"], 
            ylim=c(1,150), xlim=c(1,300),level=.5, col=standard_12[1],  
            center.pch=FALSE, 
            fill=TRUE, fill.alpha=0.1,plot.points=F, add = T)
text( 150, 55, "PEL", col=standard_12[1])
# pel <- pops.indiv[pops.indiv$pop=="PEL", c(1, 5, 8)]
# pel.2 <- coord.ellipse(pel, level.conf = 0.5, npoint = 700, bary=FALSE)
# points(pel.2$res[, 2:3], pch=".", cex=1, col=standard_12[1])
# polygon(pel.2$res[, 2:3], border=NA, col=adjustcolor(standard_12[1], alpha.f=0.2))

## Puerto Ricans: PUR
points( NSEG~Mb, data=pops.indiv[pops.indiv$pop=="PUR", ], 
        col=adjustcolor(standard_12[11], alpha.f=1), pch=".",
        cex=2, ylim=c(1,150), xlim=c(1,300), xlab="Size (Mb)", ylab="Number of segments" )
dataEllipse(pops.indiv$Mb[pops.indiv$pop=="PUR" ], pops.indiv$NSEG[pops.indiv$pop=="PUR"], 
            ylim=c(1,150), xlim=c(1,300),level=.5, col=standard_12[11],  
            center.pch=FALSE, ellipse.label ="PUR",
            fill=TRUE, fill.alpha=0.1,plot.points=F, add = T)
# pur <- pops.indiv[pops.indiv$pop=="PUR", c(1, 5, 8)]
# pur.2 <- coord.ellipse(pur, level.conf = 0.5, npoint = 700, bary=FALSE)
# points(pur.2$res[, 2:3], pch=".", cex=1, col=standard_12[11])
# polygon(pur.2$res[, 2:3], border=NA, col=adjustcolor(standard_12[11], alpha.f=0.2))
text( 40, 105, "PUR", col=standard_12[11], font=2)

# add grid
grid()

dev.off()

##################################################

## Plot Brazilian vs.Europeean populations

tiff(filename = "ROH_NSEG_BRvsEUR.tiff", width = 200, height = 200, units = "mm",
     compression = "lzw", res = 1500, type = "cairo", antialias = "gray")

par(tck = -.01)

## Brazilian: BR
plot( NSEG~Mb, data=pops.indiv[pops.indiv$pop=="BR", ], col=adjustcolor(standard_12[3], alpha.f=1), pch=".", cex=2, ylim=c(1,350), xlim=c(1,400), xlab="Mb" )

br <- pops.indiv[pops.indiv$pop=="BR", c(1, 5, 8)]
br.2 <- coord.ellipse(br, level.conf = 0.5, npoint = 700, bary=FALSE)
points(br.2$res[, 2:3], pch=".", cex=1, col=standard_12[3])
polygon(br.2$res[, 2:3], border=NA, col=adjustcolor(standard_12[3], alpha.f=0.2))
# text( 48, 35, "BR", col=standard_12[3], font=2)

## Finnish: FIN
points( NSEG~Mb, data=pops.indiv[pops.indiv$pop=="FIN", ], col=adjustcolor("blue", alpha.f=1), pch=".", cex=2, ylim=c(1,350), xlim=c(1,400), xlab="Mb" )

fin <- pops.indiv[pops.indiv$pop=="FIN", c(1, 5, 8)]
fin.2 <- coord.ellipse(fin, level.conf = 0.5, npoint = 700, bary=FALSE)
points(fin.2$res[, 2:3], pch=".", cex=1, col="blue")
polygon(fin.2$res[, 2:3], border=NA, col=adjustcolor("blue", alpha.f=0.2))
# text( 48, 35, "BR", col=standard_12[3], font=2)

## British: GBR
points( NSEG~Mb, data=pops.indiv[pops.indiv$pop=="GBR", ], col=adjustcolor("blue", alpha.f=1), pch=".", cex=2, ylim=c(1,350), xlim=c(1,400), xlab="Mb" )

gbr <- pops.indiv[pops.indiv$pop=="GBR", c(1, 5, 8)]
gbr.2 <- coord.ellipse(gbr, level.conf = 0.5, npoint = 700, bary=FALSE)
points(gbr.2$res[, 2:3], pch=".", cex=1, col="blue")
polygon(gbr.2$res[, 2:3], border=NA, col=adjustcolor("blue", alpha.f=0.2))

## Iberians from Spain: IBS
points( NSEG~Mb, data=pops.indiv[pops.indiv$pop=="IBS", ], col=adjustcolor("blue", alpha.f=1), pch=".", cex=2, ylim=c(1,350), xlim=c(1,400), xlab="Mb" )

ibs <- pops.indiv[pops.indiv$pop=="IBS", c(1, 5, 8)]
ibs.2 <- coord.ellipse(ibs, level.conf = 0.5, npoint = 700, bary=FALSE)
points(ibs.2$res[, 2:3], pch=".", cex=1, col="blue")
polygon(ibs.2$res[, 2:3], border=NA, col=adjustcolor("blue", alpha.f=0.2))

## Toscani in Italia: TSI
points( NSEG~Mb, data=pops.indiv[pops.indiv$pop=="TSI", ], col=adjustcolor("blue", alpha.f=1), pch=".", cex=2, ylim=c(1,350), xlim=c(1,400), xlab="Mb" )

tsi <- pops.indiv[pops.indiv$pop=="TSI", c(1, 5, 8)]
tsi.2 <- coord.ellipse(tsi, level.conf = 0.5, npoint = 700, bary=FALSE)
points(tsi.2$res[, 2:3], pch=".", cex=1, col="blue")
polygon(tsi.2$res[, 2:3], border=NA, col=adjustcolor("blue", alpha.f=0.2))

grid()

title("Brazilian vs. European populations")
legend(0,350, legend=c("BR","EUR"), col=c( adjustcolor(standard_12[3], alpha.f=1), adjustcolor("blue", alpha.f=0.2) ), pch=20, pt.cex=2)

dev.off()

##################################################

## Plot Brazilian vs. African populations

tiff(filename = "ROH_NSEG_BRvsAFR.tiff", width = 200, height = 200, units = "mm",
     compression = "lzw", res = 1500, type = "cairo", antialias = "gray")

par(tck = -.01)

## Brazilian: BR
plot( NSEG~Mb, data=pops.indiv[pops.indiv$pop=="BR", ], col=adjustcolor(standard_12[3], alpha.f=1), pch=".", cex=2, ylim=c(1,350), xlim=c(1,400), xlab="Mb" )

br <- pops.indiv[pops.indiv$pop=="BR", c(1, 5, 8)]
br.2 <- coord.ellipse(br, level.conf = 0.5, npoint = 700, bary=FALSE)
points(br.2$res[, 2:3], pch=".", cex=1, col=standard_12[3])
polygon(br.2$res[, 2:3], border=NA, col=adjustcolor(standard_12[3], alpha.f=0.2))

## African Caribbeans in Barbados: ACB
points( NSEG~Mb, data=pops.indiv[pops.indiv$pop=="ACB", ], col=adjustcolor(standard_12[4], alpha.f=1), pch=".", cex=2, ylim=c(1,350), xlim=c(1,400), xlab="Mb" )

acb <- pops.indiv[pops.indiv$pop=="ACB", c(1, 5, 8)]
acb.2 <- coord.ellipse(acb, level.conf = 0.5, npoint = 700, bary=FALSE)
points(acb.2$res[, 2:3], pch=".", cex=1, col=standard_12[4])
polygon(acb.2$res[, 2:3], border=NA, col=adjustcolor(standard_12[4], alpha.f=0.2))

## Gambian in Western Divisions in the Gambia: GWD
points( NSEG~Mb, data=pops.indiv[pops.indiv$pop=="GWD", ], col=adjustcolor(standard_12[4], alpha.f=1), pch=".", cex=2, ylim=c(1,350), xlim=c(1,400), xlab="Mb" )

gwd <- pops.indiv[pops.indiv$pop=="GWD", c(1, 5, 8)]
gwd.2 <- coord.ellipse(gwd, level.conf = 0.5, npoint = 700, bary=FALSE)
points(gwd.2$res[, 2:3], pch=".", cex=1, col=standard_12[4])
polygon(gwd.2$res[, 2:3], border=NA, col=adjustcolor(standard_12[4], alpha.f=0.2))

## Esan in Nigeria: ESN
points( NSEG~Mb, data=pops.indiv[pops.indiv$pop=="ESN", ], col=adjustcolor(standard_12[4], alpha.f=1), pch=".", cex=2, ylim=c(1,350), xlim=c(1,400), xlab="Mb" )

esn <- pops.indiv[pops.indiv$pop=="ESN", c(1, 5, 8)]
esn.2 <- coord.ellipse(esn, level.conf = 0.5, npoint = 700, bary=FALSE)
points(esn.2$res[, 2:3], pch=".", cex=1, col=standard_12[4])
polygon(esn.2$res[, 2:3], border=NA, col=adjustcolor(standard_12[4], alpha.f=0.2))

## Mende in Sierra Leone: MSL
points( NSEG~Mb, data=pops.indiv[pops.indiv$pop=="ESN", ], col=adjustcolor(standard_12[4], alpha.f=1), pch=".", cex=2, ylim=c(1,350), xlim=c(1,400), xlab="Mb" )

msl <- pops.indiv[pops.indiv$pop=="MSL", c(1, 5, 8)]
msl.2 <- coord.ellipse(msl, level.conf = 0.5, npoint = 700, bary=FALSE)
points(msl.2$res[, 2:3], pch=".", cex=1, col=standard_12[4])
polygon(msl.2$res[, 2:3], border=NA, col=adjustcolor(standard_12[4], alpha.f=0.2))

## Luhya in Webuye, Kenya: LWK
points( NSEG~Mb, data=pops.indiv[pops.indiv$pop=="LWK", ], col=adjustcolor(standard_12[4], alpha.f=1), pch=".", cex=2, ylim=c(1,350), xlim=c(1,400), xlab="Mb" )

lwk <- pops.indiv[pops.indiv$pop=="ESN", c(1, 5, 8)]
lwk.2 <- coord.ellipse(lwk, level.conf = 0.5, npoint = 700, bary=FALSE)
points(lwk.2$res[, 2:3], pch=".", cex=1, col=standard_12[4])
polygon(lwk.2$res[, 2:3], border=NA, col=adjustcolor(standard_12[4], alpha.f=0.2))

grid()

title("Brazilian vs. African populations")
legend(0,350, legend=c("BR","AFR"), col=c( adjustcolor(standard_12[3], alpha.f=1), adjustcolor(standard_12[4], alpha.f=0.2) ), pch=20, pt.cex=2)

dev.off()


