# 15/11/2015
# CLUMPP plot
########################################

## define working directory
# setwd("C:/Users/Galina/Documents/KROTIK/R-projects/20140219_mydata")
setwd("/home/krotik/Dropbox/R-projects_current/mydata_current")

## read data files
# tbl <- read.table("Clumpp/clumpp_runs/MY_2.run1.outfile")
tbl <- read.table("/home/krotik/Backup/Krolik_2016-03-01/KROTIK/R-projects/20140219_mydata/Clumpp/clumpp_runs/MY_2.run1.outfile")

# pops <- read.table("Clumpp/geno002-filt_recoded_1.ids")
pops <- read.table("/home/krotik/Backup/Krolik_2016-03-01/KROTIK/R-projects/20140219_mydata/Clumpp/geno002-filt_recoded_1.ids")

## remove a bunch of extraneous columns
tbl.2 <- within(tbl, rm(V1,V2,V3,V4,V5))

## re-order columns (customize the order)
tbl.2<-tbl.2[, c("V6","V7","V9","V10","V8","V11","V12")]

## re-order column by mean
## tbl.2<-tbl.2[,order(colMeans(tbl.2), decreasing=TRUE)]

## change population names for Brazilian population
## BRcamp (and others) -> BR
pops$V2 <- gsub(pattern="BR.+" , replacement="BR" ,x=pops$V2)

## add to dataframe a column with populations' names which is factor 
## the levels are organized in the order, in which they first occur in the data 
tbl.2$Pop <- with(tbl.2, factor(pops$V2, levels=unique(pops$V2)))
#--------------------------------------------------------
## re-order rows separately for each population (V8) in decrasing order
## from Pedro:
tbl.2 <- tbl.2[order(tbl.2$Pop, -(tbl.2$V6)), ] # tbl.2$V6 is South Europe component in this dataset
rownames(tbl.2) <- NULL
#--------------------------------------------------------
## constructing graph
mp <- barplot(t(as.matrix(tbl.2[,1:7])), col=rainbow(7), 
        xlab="Individual #", ylab="Ancestry", border=NA, plot=FALSE)

#--------------------------------------------------------
## Discovering Components (initial) order 
## V6 - South European
## V7 - North European
## V8 - MKK main component (African, non-Yorubian)
## V9 - Amerindian
## V10 - Yorubian component
## V11 - Asian component (jpt, chinese)
## V12 - Indian component (India)

#--------------------------------------------------------
## Cusomize colors
col=colorRampPalette(c(rgb(0.1, 0.5, 0.8), rgb(0.411764705882353,0.823529411764706,0.905882352941176), 
                       rgb(0.8,0.2,0),rgb(0.87843137254902,0.894117647058824,0.8), 
                       rgb(0.980392156862745,0.411764705882353,0)))

C1=c(168,206,191)
C1=C1/255
C2=c(235,172,203)
C2=C2/255
C3=c(99,168,162)
C3=C3/255
C4=c(75,48,67)
C4=C4/255
C5=c(201,172,168)
C5=C5/255
C6=c(229,73,112)
C6=C6/255
C7=c(15,118,114)
C7=C7/255

col=colorRampPalette(c(rgb(C1[1],C1[2],C1[3]), rgb(C2[1],C2[2],C2[3]),
                       rgb(C3[1],C3[2],C3[3]), rgb(C4[1],C4[2],C4[3]),
                       rgb(C5[1],C5[2],C5[3]),rgb(C6[1],C6[2],C6[3]),
                       rgb(C7[1],C7[2],C7[3])))

mp <- barplot(t(as.matrix(tbl.2[,1:7])), 
        col=col(7), 
        xlab="Individual #", ylab="Ancestry", border=NA)

## add to data frame a column of meadspoints of each bar
tbl.2$mp <- with(tbl.2, mp)

## define coordinates of a meadpoint for each popuation group
pop.mp <- with (tbl.2, tapply(tbl.2$mp, tbl.2$Pop, FUN = mean))

## define coordinates to draw lines to the right of eacn population group
## pop.lp <- with (tbl.2, tapply(tbl.2$mp+(tbl.2$mp[2]-tbl.2$mp[1]), tbl.2$Pop, FUN = max))

## Write text on the bottom of the figure
mtext(side=1, at=pop.mp, text=names(summary(tbl.2$Pop)), line=0, cex=0.7)

## define coordinates to draw lines to the right of eacn population group
pop.lp <- with (tbl.2, tapply(tbl.2$mp+(1.5*(tbl.2$mp[2]-tbl.2$mp[1]) - tbl.2$mp[1]), tbl.2$Pop, FUN = max))

## ad vertical lines delimiting each population group to the right
abline(v=pop.lp, untf = FALSE, col="black", lty=1, lwd=0.5)

#########################################################
## Customizing plot (adding black lines at the end of each popultion) 

## Prepare new data frame (iqual, but a bit different)
tbl.3 <- within(tbl, rm(V1,V2,V3,V4,V5))
tbl.3<-tbl.3[, c("V6","V7","V9","V10","V8","V11","V12")]

## Add a dummy column to the created data frame
dummy <- rep(0, times=dim(tbl.3)[1])
tbl.3$dummy <- with(tbl.3, dummy)

## Add a column with population names as a factor
tbl.3$Pop <- with(tbl.3, factor(pops$V2, levels=unique(pops$V2)))

## Order by the first column values inside each popuation (without changing population order)
tbl.3 <- tbl.3[order(tbl.3$Pop, -(tbl.3$V6)), ]
rownames(tbl.3) <- NULL


## add extra lines to the table after each population
### Create a row with equal column names
Row <- tbl.3[1, ]
### Change all values, but population name
Row[, 1:8] <- c(rep(0,7), 1)
### Repeat this row n times (n iqual to the number of populations) and append one to each other
df <- do.call("rbind", replicate(16, Row, simplify = FALSE))
levels=unique(tbl.3$Pop)
df$Pop <- levels

### Now, we will merge the principal table with the created data frame
tbl.3 <- rbind(tbl.3, df)

### Create one more column wtih index
tbl.3$index <- 1:nrow(tbl.3)

### Order by population name
tbl.3 <- tbl.3[order(tbl.3$Pop), ]
tbl.3 <- tbl.3[, -ncol(tbl.3)]
rownames(tbl.3) <- NULL
#--------------------------------------------------------
## constructing graph
mp <- barplot(t(as.matrix(tbl.3[,1:8])), col=c(rainbow(7), "black"),
        xlab="Individual #", ylab="Ancestry", border=NA, plot=FALSE)

#--------------------------------------------------------
## add to data frame a column of meadspoints of each bar
tbl.3$mp <- with(tbl.3, mp)

## define coordinates of a meadpoint for each popuation group
pop.mp <- with (tbl.3, tapply(tbl.3$mp, tbl.3$Pop, FUN = mean))

## Write text on the bottom of the figure
mtext(side=1, at=pop.mp, text=names(summary(tbl.3$Pop)), line=0, cex=0.7)

#--------------------------------------------------------

tbl.4 <- within(tbl, rm(V1,V2,V3,V4,V5))
tbl.4<-tbl.4[, c("V6","V7","V9","V10","V8","V11","V12")]
tbl.4$Pop <- with(tbl.4, factor(pops$V2, levels=unique(pops$V2)))

aggregated <- aggregate(tbl.4[,1:7], list(tbl.4$Pop), FUN=mean)
aggregated <- as.matrix(aggregated)[, -1]
barplot(t(aggregated), col=rainbow(7))

a <- (summary(tbl.4$Pop))
names.arg <- names(a)
names(a) <- NULL

barplot(t(aggregated), width=a, names.arg=names.arg, col=rainbow(7), xlab="")

################################################################
medians <- with (tbl.2, tapply(tbl.2, tbl.2[,8], FUN = median)) ##
## boxplots
grp=factor(pops$V2, levels=unique(pops$V2))
boxplot(tbl.4$V6~grp, width=a, names.arg=names.arg, col=rainbow(7)[1]) ## South European component
boxplot(tbl.4$V7~grp, width=a, names.arg=names.arg, col=rainbow(7)[2]) ## North European component
boxplot(tbl.4$V9~grp, width=a, names.arg=names.arg, col=rainbow(7)[3]) ## Amerindian component
boxplot(tbl.4$V10~grp, width=a, names.arg=names.arg, col=rainbow(7)[4]) ## Yorubian component

boxplot((tbl.4$V6+tbl.4$V7)~grp, width=a, names.arg=names.arg, col=rainbow(20)[2]) ## South+North European components

################################################################

## All graphs together

png("Admixture_Clumpp_plot.png",  width = 2600, height = 2600, res=300, type = "cairo", antialias = "none", family="serif")
par(mfrow = c(5, 1))
par(mar = c(2.1, 3.1, 1.1, 0), oma = c(1, 1, 1, 1))
par(tcl=-0.25)
par(xpd=TRUE)

## 1-st graph
mp <- barplot(t(as.matrix(tbl.3[,1:8])), col=c(rainbow(7), "black"),
        xlab="Individual #", ylab="Ancestry", border=NA)
## Write text on the bottom of the figure
mtext(side=1, at=pop.mp, text=names(summary(tbl.3$Pop)), line=0, cex=0.7)
legend("topright", inset=c(-0.005, 0), legend="(A)", cex=1.2, bty = "n")

## 2-nd graph
mp <- barplot(t(aggregated), width=a, names.arg=names.arg, col=rainbow(7), xlab="", mgp=c(3,0,0), axes=FALSE)
axis(side=2, at=seq(0,1,0.2), labels=c("0.0",0.2,0.4,0.6,0.8,"1.0"))
legend("topright", inset=c(-0.005, 0), legend="(B)", cex=1.2, bty = "n")

## 3-rd graph
## South+North European components
boxplot((tbl.4$V6+tbl.4$V7)~grp, width=a, names.arg=names.arg, col=rainbow(20)[2], at=(mp*16/mp[16])+0.5, axes=FALSE) # medlwd - thickness of midline
axis(1, at=(mp*16/mp[16])+0.5, labels=names.arg, mgp=c(3,0.1,0), cex.axis=1.1, tcl=0) 
axis(1, at=0.76:17.759, labels=FALSE, mgp=c(3,0.1,0), tcl=0)
axis(side=2, at=seq(0,1,0.2), labels=c("0.0",0.2,0.4,0.6,0.8,"1.0"))
legend("topright", inset=c(-0.005, 0), legend="(C)", cex=1.2, bty = "n")
# box(col = "grey60")

## 4-th graph
## Amerindian component
boxplot(tbl.4$V9~grp, width=a, names.arg=names.arg, col=rainbow(7)[3], at=(mp*16/mp[16])+0.5, axes=FALSE)
axis(1, at=(mp*16/mp[16])+0.5, labels=names.arg, mgp=c(3,0.1,0), cex.axis=1.1, tcl=0)
axis(1, at=0.76:17.759, labels=FALSE, mgp=c(3,0.1,0), tcl=0)
axis(side=2, at=seq(0,1,0.2), labels=c("0.0",0.2,0.4,0.6,0.8,"1.0"))
legend("topright", inset=c(-0.005, 0), legend="(D)", cex=1.2, bty = "n")
# box(col = "grey60")

## 5-th graph
## Yorubian component
boxplot(tbl.4$V10~grp, width=a, names.arg=names.arg, col=rainbow(7)[4], at=(mp*16/mp[16])+0.5, axes=FALSE) 
axis(1, at=(mp*16/mp[16])+0.5, labels=names.arg, mgp=c(3,0.1,0), cex.axis=1.1, tcl=0)
axis(1, at=0.76:17.759, labels=FALSE, mgp=c(3,0.1,0), tcl=0)
axis(side=2, at=seq(0,1,0.2), labels=c("0.0",0.2,0.4,0.6,0.8,"1.0"))
legend("topright", inset=c(-0.005, 0), legend="(E)", cex=1.2, bty = "n")
# box(col = "grey60")

dev.off() 
