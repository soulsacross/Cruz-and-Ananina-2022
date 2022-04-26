#------------------------ my data -----------------------------
# This contains iHS values in all populations
iHS <- read.table("/Users/Pedro/Dropbox/paperultima/iHS validation/all.iHS.ever.txt", header=F)

plot(sort(x = iHS$V1), ylim=c(-7.0, 7.0)) # this is empirical
par(new=TRUE)
plot(sort(x = iHS$V2[iHS$V2 > 1 | iHS$V2 < -1]), ylim=c(-7.0, 7.0), col="orange2", axes=FALSE) # this is simulated
par(new=TRUE)

# plot a double exponential curve (the distribution)
plot(x=sort(toy),ylim=c(-7.0, 7.0), type="l")


#####========= Now for individual populations =========
# set order, leaving first row and column empty to add Ld decay in gimp
layout.matrix <- matrix(c(0, 2, 4, 1, 3, 5), nrow = 3, ncol = 2)
layout(mat = layout.matrix,
               heights = c(1, 1), # Heights of the two rows
               widths = c(1, 1))
# layout.show(5) # check layout

# BR - 11066 empirical (normalized and suggestive, i.e. |iHS|>2) signals
br <- read.table("~/Dropbox/input_chr16/haplotype.scores/BR.signals.txt", header=F)
#br.sim <- as.data.frame(cbind(br$V7, iHS[10187:21252,2])) # df
br.sim <-list("Empirical"=br$V7, "Simulation"=iHS[4505:6996,2]) # list

#plot(sort(x=br.sim$V1), ylim=c(-5.0, 6.0)) # this is empirical
plot(sort(unlist(br.sim["Empirical"])), ylim=c(-5.0, 6.0), ylab = "iHS", xlab = "value count") # if list
par(new=TRUE)

#plot(sort(x=br.sim$V2), ylim=c(-5.0, 6.0), col="orange2", axes=FALSE) # this is simulated
plot(sort(unlist(br.sim["Simulation"])), ylim=c(-5.0, 6.0), col="orange2", axes=FALSE, ylab = "", xlab = "") # if list
title("BR iHS values vs. simulated", cex.main=1.1)

### PUR - 10186 signals
pr <- read.table("~/Dropbox/input_chr16/haplotype.scores/PUR.signals.txt", header=F)
# pr.sim <- as.data.frame(cbind(pr$V7, iHS[1:10186,2])) #df
pr.sim <- list("Empirical"=pr$V7, "Simulation"=iHS[6997:10133,2]) # list

#plot(sort(x=pr.sim$V1), ylim=c(-5.0, 6.0)) # this is empirical
plot(sort(unlist(pr.sim["Empirical"])), ylim=c(-5.0, 6.0), ylab = "iHS", xlab = "value count") # if list
par(new=TRUE)

#plot(sort(x=pr.sim$V2), ylim=c(-5.0, 6.0), col="orange2", axes=FALSE) # this is empirical
plot(sort(unlist(pr.sim["Simulation"])), ylim=c(-5.0, 6.0), col="orange2", axes=FALSE, ylab = "", xlab = "") # if list
title("PUR iHS values vs. simulated", cex.main=1.1)

##### CLM - 9804 signals
cl <- read.table("~/Dropbox/input_chr16/haplotype.scores/CLM.signals.txt", header=F)
# cl.sim <- as.data.frame(cbind(cl$V7, iHS[21253:31056,2])) #df
cl.sim <- list("Empirical"=cl$V7, "Simulation"=iHS[10134:13151,2]) # list

#plot(sort(x=cl.sim$V1), ylim=c(-5.0, 6.0)) # this is empirical
plot(sort(unlist(cl.sim["Empirical"])), ylim=c(-5.0, 6.0), ylab = "iHS", xlab = "value count") # if list
par(new=TRUE)

#plot(sort(x=cl.sim$V2), ylim=c(-5.0, 6.0), col="orange2", axes=FALSE) # this is empirical
plot(sort(unlist(cl.sim["Simulation"])), ylim=c(-5.0, 6.0), col="orange2", axes=FALSE, ylab = "", xlab = "") # if list
title("CLM iHS values vs. simulated", cex.main=1.1)

###### MXL - 9554 signals
mx <- read.table("~/Dropbox/input_chr16/haplotype.scores/MXL.signals.txt", header=F)
#mx.sim <- as.data.frame(cbind(mx$V7, iHS[5534:15087,2])) #df
mx.sim <- list("Empirical"=mx$V7, "Simulation"=iHS[2236:4504,2]) # list

# plot(sort(x=mx.sim$V1), ylim=c(-5.0, 6.0)) # this is empirical
plot(sort(unlist(mx.sim["Empirical"])), ylim=c(-5.0, 6.0), ylab = "iHS", xlab = "value count") # if list
par(new=TRUE)

#plot(sort(x=mx.sim$V2), ylim=c(-5.0, 6.0), col="orange2", axes=FALSE) # this is empirical
plot(sort(unlist(mx.sim["Simulation"])), ylim=c(-5.0, 6.0), col="orange2", axes=FALSE, ylab = "", xlab = "") # if list
title("MXL iHS values vs. simulated", cex.main=1.1)

#####PEL - 8606 signals
pl <- read.table("~/Dropbox/input_chr16/haplotype.scores/PEL.signals.txt", header=F)
#pl.sim <- as.data.frame(cbind(pl$V7, iHS[15089:23694,2])) # if df OR 
pl.sim <- list("Empirical"=pl$V7, "Simulation"=iHS[1:2235,2]) # list

#plot(sort(x=pl.sim$V1), ylim=c(-5.0, 6.0)) # this is empirical # if df OR
plot(sort(unlist(pl.sim["Empirical"])), ylim=c(-5.0, 6.0), ylab = "iHS", xlab = "value count") # if list
par(new=TRUE)

#plot(sort(x=pl.sim$V2), ylim=c(-5.0, 6.0), col="orange2", axes=FALSE) # this is empirical OR
plot(sort(unlist(pl.sim["Simulation"])), ylim=c(-5.0, 6.0), col="orange2", axes=FALSE, ylab = "", xlab = "") # if list
title("PEL iHS values vs. simulated", cex.main=1.1)


#####========= Now using ggplot =========
library(ggplot2)
library(gridExtra)
# Important note: to set colors and legend to colors manually, define color in aes() - see "plot
# a sole panel" (the colors <- c("...") is necessary in this case), if not plotting legend, set color
# outside aes() as "col="
## Observation: the best practice is to use tydr to rearrange values to a single column,
# then map to group in ggplot (see http://www.sthda.com/english/wiki/tidyr-crucial-step-reshaping-data-with-r-for-easier-analyses)

# Assinging colors
colnames(br.sim) <- c("Empirical", "Simulated")
colnames(pr.sim) <- c("Empirical", "Simulated")
colnames(mx.sim) <- c("Empirical", "Simulated")
colnames(pl.sim) <- c("Empirical", "Simulated")
colnames(cl.sim) <- c("Empirical", "Simulated")
colors <- c("Empirical" = "red", "Simulated" = "black") # this is used to plot the legend manually

##### to plot a sole panel (with legend)
ggplot() +
  geom_point(data = br.sim, aes(x=1:nrow(br.sim), y=sort(Simulated), color="Simulated"), size=1) +
  geom_point(data = br.sim, aes(x=1:nrow(br.sim), y=sort(Empirical), color="Empirical"), size=1) +
  labs(title="Empirical vs simulation values", 
       subtitle="Brazilians", 
       color = "iHS value",
       y="iHS",
       x="value count") +
  scale_color_manual(values = colors) +
  theme_classic()


########### to plot as a multi-panel plot

p1 <- ggplot() +
  geom_point(data = br.sim, aes(x=1:nrow(br.sim), y=sort(Simulated)), col="black", size=1) +
  geom_point(data = br.sim, aes(x=1:nrow(br.sim), y=sort(Empirical)), col="red", size=1) +
  labs(title="Empirical vs simulation values", 
       subtitle="Brazilians", 
       #       color = "iHS value",
       y="iHS",
       x="value count") +
  #  scale_color_manual(values = colors) +
  theme_classic()
#----------------------------------

p2 <- ggplot() +
  geom_point(data = pr.sim, aes(x=1:nrow(pr.sim), y=sort(Simulated)), col="black",size=1) +
  geom_point(data = pr.sim, aes(x=1:nrow(pr.sim), y=sort(Empirical)), col="red", size=1) +
  labs(
    #title="Empirical vs simulation values", 
    subtitle="Purto Ricans", 
    y="iHS",
    x="value count",
    color = "iHS value") +
  #  scale_color_manual(values = colors) +
  theme_classic()

#----------------------------------

p3 <- ggplot() +
  geom_point(data = cl.sim, aes(x=1:nrow(cl.sim), y=sort(Simulated)),col="black", size=1) +
  geom_point(data = cl.sim, aes(x=1:nrow(cl.sim), y=sort(Empirical)), col="red", size=1) +
  labs(
    #title="Empirical vs simulation values", 
    subtitle="Colombians", 
    y="iHS",
    x="value count",
    color = "iHS value") +
  scale_color_manual(values = colors) +
  theme_classic()

#----------------------------------

p4 <-  ggplot() +
  geom_point(data = mx.sim, aes(x=1:nrow(mx.sim), y=sort(Simulated), color="Simulated"), size=1) +
  geom_point(data = mx.sim, aes(x=1:nrow(mx.sim), y=sort(Empirical), color="Empirical"), size=1) +
  labs(
    #title="Empirical vs simulation values", 
    subtitle="Mexicans", 
    y="iHS",
    x="value count",
    color = "iHS value") +
  scale_color_manual(values = colors) +
  theme_classic()

#----------------------------------

p5 <- ggplot() +
  geom_point(data = pl.sim, aes(x=1:nrow(pl.sim), y=sort(Simulated)),col="black", size=1) +
  geom_point(data = pl.sim, aes(x=1:nrow(pl.sim), y=sort(Empirical)), col="red", size=1) +
  labs(
    #title="Empirical vs simulation values", 
    subtitle="Peruvians", 
    y="iHS",
    x="value count",
    color = "iHS value") +
  #  scale_color_manual(values = colors) +
  theme_classic()

grid.arrange(p1, p2, p3, p4, p5, nrow=3, ncol=2)


