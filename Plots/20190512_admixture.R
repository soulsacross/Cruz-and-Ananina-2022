# May 12th 2019 (Pedro)
# ADMIXTURE plot 
########################################
library(RColorBrewer)
library(randomcoloR)
library(viridis)
library(Cairo)

n <- 21
palette <- distinctColorPalette(n)
pie(rep(1, n), col=palette)

#setwd("/Users/Pedro/Dropbox/Documentos/OneDrive/Projeto_Doutorado/Análise/2013.04.03_plink/SNPs/Resultados/ADMIXTURE/")
setwd("/Users/Pedro/Documents/Sickle Paper/Original Files/updated/rev")

## tbl <- read.table("admixture/allpops-pca0.recoded.7.Q")
##tbl <- read.table("__allpops-pca0.LDfiltered.recoded.3.Q")
##tbl <- read.table("__ALL2016FINAL.4.Q")
#tbl <- read.table("__ALL.final12coded.6.Q")
# Without Joinville
#tbl <- read.table("/Users/Pedro/Documents/Sickle Paper/Original Files/updated/__ALL.final12coded.5.Q")
# Without Joinville and Recife PLUS JPT and CHB
tbl <- read.table("/Users/Pedro/Documents/Sickle Paper/Original Files/updated/rev/__ALL.final12coded.6.Q") 
#tbl <- read.table("/Users/Pedro/Documents/Sickle Paper/Original Files/updated/__ALL.final12coded.6.Q")

### OPTIONAL: order consistently to Galina's plot
## If done so, then now V1 is Native American
#				   V2 is still North European
#				   V3 is Bantu (LWK) 
#				   V4 is S. Eur.
#			 	   V5 is Mandé (GWD)
#tbl<-tbl[, c("V5","V2","V4","V1","V3")]
#tbl<-tbl[, c("V5","V2","V4","V1","V3", "V6")]

##pops <- read.table("PED_allpops-pca0.PEDRO")
##pops <- read.table("ALL.2016.indlist")
#pops <- read.table("ALL.wo.Jionv.indilst")
#pops <- read.table("ALL.wo.Jionv.Rcf.indlist")
# Also for vertical labels
#pops <- read.table("/Users/Pedro/Documents/Sickle Paper/Original Files/updated/rev/ALL.indlist.txt")
pops <- read.table("ALL.indlist.Vertical-spaced.txt")
pops$V6 <- gsub("\\."," ", pops$V6)

### OPTIONAL: order consistently to Galina's plot
# tbl2<-tbl[, c("V5","V2","V4","V1","V3")]  # VERY CLOSE
# tbl2<-tbl[, c("V5","V2","V1","V4","V3")]
# tbl2<-tbl[, c("V5","V2","V4","V3","V1")]
# tbl2<-tbl[, c("V5","V3","V2","V1","V4")]   # VERY CLOSE
# tbl2<-tbl[, c("V5","V2","V3","V1","V4")]
# tbl2<-tbl[, c("V5","V2","V3","V1","V4")]
# tbl2<-tbl[, c("V5","V4","V3","V1","V2")]
# tbl2<-tbl[, c("V5","V4","V2","V1","V3")]
# tbl2<-tbl[, c("V5","V4","V2","V1","V4")]
# tbl2<-tbl[, c("V5","V3","V4","V1","V2")]  # VERY CLOSE << chosen by now
# tbl2<-tbl[, c("V1","V3","V2","V5","V4")] 

tbl2<-tbl[, c("V6","V5","V3","V4","V1","V2")]

## add to dataframe a column with populations' names which is factor 
## the levels are organized in the order, in which they first occur in the data 
tbl2$Pop <- with(tbl2, factor(pops[,6], levels=unique(pops[,6])))

# Nice trick to reorder populations:
tbl2$Pop <- ordered(tbl2$Pop, levels = c("FIN", "CEU", "GBR","IBS", "TSI", "NAM", "PEL", "MXL", "CLM", "PUR", "BRZ", "SBR", "SUS", "AAM","ASW","ACB","YRI", "ESN","MSL" ,"GWD","LWK", "CHB", "JPT"))
tbl2$Pop <- ordered(tbl2$Pop, levels = c("Finnish", "N/W European", "Britain","Iberian", "Tuscan", "Aymara Quechua", "Peruvian", "Mexican", "Colombian", "Portorican", "Brazilian", "Sickle Brazil", "Sickle US", "African American","SW Afr American","Barbadian","Yoruban", "Esan","Mende" ,"Gambian","Luhya", "Japanese", "Han Chinese" ))

## Through this link (http://stackoverflow.com/questions/17374363/sorting-stacked-bar-plot-by-factor-in-dataframe) I have discovered how to 
## sort by ancestry component so the plot looks cleaner (of course this donnot mess with right populations' order). To do so, run the following line
tbl2 <- tbl2[order(tbl2$Pop,tbl2$V5),] # tbl$V5 is Native Arican component in this dataset
mp <- barplot(t(as.matrix(tbl2[,1:6])), col=rainbow(7), 
        xlab="Individual #", ylab="Ancestry", border=NA, axisnames=F)

## add to data frame a column of measpoints of each bar
tbl2$mp <- with(tbl2, mp)

## define coordinates of a meadpoint for each popuation group
pop.mp <- with (tbl2, tapply(tbl2$mp, tbl2$Pop, FUN = mean))

## define coordinates to drae lines to the right of eacn population group
pop.lp <- with (tbl, tapply(tbl2$mp+(tbl2$mp[2]-tbl2$mp[1]), tbl2$Pop, FUN = max))

## Write text on the bottom of the figure
mtext(side=1, at=pop.mp, text=names(summary(tbl2$Pop)), las=2, line=0, cex=0.5)

## define coordinates to drae lines to the right of eacn population group
pop.lp <- with (tbl2, tapply(tbl2$mp+(1.5*(tbl2$mp[2]-tbl2$mp[1]) - tbl2$mp[1]), tbl2$Pop, FUN = max))

## ad vertical lines delimiting each population group to the right
abline(v=pop.lp, untf = FALSE, col="black", lty=1, lwd=0.5)

#============= Lets now do the graph with the means ========================
#Add the following if saving to pdf
pdf("ADMIXTURE2018.MEAN.pdf", width=24, height=4)

aggregated <- aggregate(tbl2[,1:6], list(tbl2$Pop), FUN=mean)
aggregated <- as.matrix(aggregated)[, -1]
barplot(t(aggregated), col=rainbow(8))
a <- (summary(tbl2$Pop))
names.arg <- names(a)
names(a) <- NULL
barplot(t(aggregated),names.arg=names.arg, cex.names=0.595, las =2, col=c("#FF0000FF", "#FFDB00FF", "#49FF00FF", "#6FA661", "#0092FFFF", "#4900FFFF", "#FF00DBFF"), xlab="") # las=2 makes vertical labels
                                                    #col=rainbow(7)
dev.off()

#=========== Alternatively one could print directy on .pdf ==================
pdf("ADMIXTURE2018.pdf", width=10, height=4)
tbl <- tbl[order(tbl$Pop,tbl$V5),]
barplot(t(as.matrix(tbl[,1:5])), col=rainbow(5), 
        xlab="Individual #", ylab="Ancestry", border=NA, axisnames=F)

mtext(side=1, at=pop.mp, text=names(summary(tbl$Pop)), line=0, cex=0.5)

## define coordinates to drae lines to the right of eacn population group
pop.lp <- with (tbl, tapply(tbl$mp+(1.5*(tbl$mp[2]-tbl$mp[1]) - tbl$mp[1]), tbl$Pop, FUN = max))

## ad vertical lines delimiting each population group to the right
abline(v=pop.lp, untf = FALSE, col="black", lty=1, lwd=0.5)

dev.off()
#=========== testing ancestrality components in the populations==============
# Testing wo Joinville and wo Recife and JPT+CHB
media <- with (tbl2, tapply(tbl2$V1, tbl2$Pop, FUN = mean))   # FIN
media
Finnish     N/W.European          Britain          Iberian           Tuscan 
0.825846111      0.380971273      0.377401901      0.152137869      0.100329561 
Aymara.Quechua         Peruvian          Mexican        Colombian       Portorican 
0.005323256      0.022261365      0.053073656      0.051994255      0.056201519 
Brazilian    Sickle.Brazil        Sickle.US African.American  SW.Afr.American 
0.074415627      0.062929000      0.070552964      0.060521683      0.078192131 
Barbadian          Yoruban             Esan            Mende          Gambian 
0.049132635      0.008411222      0.007093889      0.008709706      0.003065655 
Luhya         Japanese      Han.Chinese 
0.002460889      0.011020337      0.016437184 

media <- with (tbl2, tapply(tbl2$V2, tbl2$Pop, FUN = mean)) #JPT+CHB
media

Finnish     N/W.European          Britain          Iberian           Tuscan 
0.037908293      0.006232505      0.004921473      0.014549972      0.022122271 
Aymara.Quechua         Peruvian          Mexican        Colombian       Portorican 
0.024025930      0.019977612      0.045878281      0.024161074      0.020034740 
Brazilian    Sickle.Brazil        Sickle.US African.American  SW.Afr.American 
0.016562333      0.016053067      0.018457179      0.010932967      0.014842459 
Barbadian          Yoruban             Esan            Mende          Gambian 
0.008489521      0.003471731      0.003503616      0.004907271      0.007409407 
Luhya         Japanese      Han.Chinese 
0.016286475      0.945994433      0.937824806 

media <- with (tbl2, tapply(tbl2$V3, tbl2$Pop, FUN = mean))    #NAM
media

Finnish     N/W.European          Britain          Iberian           Tuscan 
0.007025131      0.009452737      0.009946912      0.006115832      0.007095598 
Aymara.Quechua         Peruvian          Mexican        Colombian       Portorican 
0.956390721      0.752956882      0.453545359      0.260834989      0.140613846 
Brazilian    Sickle.Brazil        Sickle.US African.American  SW.Afr.American 
0.081753471      0.099869033      0.015228107      0.010539333      0.036201705 
Barbadian          Yoruban             Esan            Mende          Gambian 
0.005515073      0.003817657      0.002919586      0.004483682      0.003212814 
Luhya         Japanese      Han.Chinese 
0.003978525      0.013960135      0.010426563  

media <- with (tbl2, tapply(tbl2$V4, tbl2$Pop, FUN = mean))   # LWK (Bantu)
media
Finnish     N/W.European          Britain          Iberian           Tuscan 
0.005530848      0.006043111      0.004629681      0.013602981      0.013416346 
Aymara.Quechua         Peruvian          Mexican        Colombian       Portorican 
0.002247233      0.017590071      0.021453453      0.040578777      0.063466096 
Brazilian    Sickle.Brazil        Sickle.US African.American  SW.Afr.American 
0.250422765      0.280042878      0.394655786      0.404327100      0.387091262 
Barbadian          Yoruban             Esan            Mende          Gambian 
0.445815260      0.549448556      0.606105535      0.244408671      0.055870460 
Luhya         Japanese      Han.Chinese 
0.871888081      0.009228587      0.010604087

media <- with (tbl2, tapply(tbl2$V5, tbl2$Pop, FUN = mean))   #GWD (west)
media
Finnish     N/W.European          Britain          Iberian           Tuscan 
0.004346869      0.005535970      0.004801758      0.016510140      0.012856654 
Aymara.Quechua         Peruvian          Mexican        Colombian       Portorican 
0.002262651      0.016020682      0.032212578      0.052670862      0.087062087 
Brazilian    Sickle.Brazil        Sickle.US African.American  SW.Afr.American 
0.132302784      0.152515667      0.383222679      0.413184600      0.361116049 
Barbadian          Yoruban             Esan            Mende          Gambian 
0.420774219      0.431117611      0.378345697      0.734854929      0.903939009 
Luhya         Japanese      Han.Chinese 
0.055903525      0.008026471      0.008772699 

media <- with (tbl2, tapply(tbl2$V6, tbl2$Pop, FUN = mean)) # IBS/TSI (S. Eur.)
media
Finnish     N/W.European          Britain          Iberian           Tuscan 
0.119342778      0.591764333      0.598298330      0.797083336      0.844179486 
Aymara.Quechua         Peruvian          Mexican        Colombian       Portorican 
0.009750233      0.171193471      0.393836781      0.569760000      0.632621692 
Brazilian    Sickle.Brazil        Sickle.US African.American  SW.Afr.American 
0.444543078      0.388590456      0.117883321      0.100494300      0.122556508 
Barbadian          Yoruban             Esan            Mende          Gambian 
0.070273292      0.003733324      0.002031626      0.002635659      0.026502611 
Luhya         Japanese      Han.Chinese 
0.049482475      0.011769923      0.015934573 

 

stdev1 <- with(tbl2, tapply(tbl2$V1, tbl2$Pop, FUN = sd))  
stdev1
Finnish     N/W.European          Britain          Iberian           Tuscan 
0.074347081      0.033402398      0.024320918      0.025164759      0.020065211 
Aymara.Quechua         Peruvian          Mexican        Colombian       Portorican 
0.009949121      0.018968968      0.038375074      0.026175833      0.031209305 
Brazilian    Sickle.Brazil        Sickle.US African.American  SW.Afr.American 
0.050215388      0.024328222      0.036152865      0.028125609      0.039122632 
Barbadian          Yoruban             Esan            Mende          Gambian 
0.025874219      0.005637843      0.004688572      0.005708386      0.005603932 
Luhya         Japanese      Han.Chinese 
0.005640789      0.010304520      0.011316010 

stdev2 <- with(tbl2, tapply(tbl2$V2, tbl2$Pop, FUN = sd))  #JPT
Aymara Quechua African American        Sickle US        Brazilian    Sickle Brazil          Britain          Finnish       Portorican        Colombian 
0.027960908      0.010236640      0.026361574      0.007684217      0.008104403      0.004671778      0.009588524      0.007843467      0.009867678 
Iberian         Peruvian        Barbadian          Gambian             Esan            Mende     N/W European          Yoruban      Han Chinese 
0.007136242      0.042513364      0.016149888      0.004036764      0.003513315      0.004428274      0.006700963      0.003418542      0.011601876 
Japanese            Luhya  SW Afr American          Mexican           Tuscan 
0.006726881      0.005122909      0.016121509      0.015803260      0.007604607

stdev3 <- with(tbl2, tapply(tbl2$V3, tbl2$Pop, FUN = sd))  
stdev3
Finnish     N/W.European          Britain          Iberian           Tuscan 
0.005454529      0.005115398      0.005648939      0.005169867      0.005368987 
Aymara.Quechua         Peruvian          Mexican        Colombian       Portorican 
0.039588656      0.153590925      0.180202997      0.092010486      0.038107102 
Brazilian    Sickle.Brazil        Sickle.US African.American  SW.Afr.American 
0.055541294      0.039024259      0.009505149      0.006674537      0.094546409 
Barbadian          Yoruban             Esan            Mende          Gambian 
0.005032924      0.003055027      0.003137342      0.003783449      0.003078797 
Luhya         Japanese      Han.Chinese 

stdev4 <- with(tbl2, tapply(tbl2$V4, tbl2$Pop, FUN = sd))  
stdev4
Finnish     N/W.European          Britain          Iberian           Tuscan 
0.005837364      0.005832575      0.005123503      0.011868312      0.010451371 
Aymara.Quechua         Peruvian          Mexican        Colombian       Portorican 
0.004245778      0.034531164      0.015220656      0.040211510      0.054396978 
Brazilian    Sickle.Brazil        Sickle.US African.American  SW.Afr.American 
0.193959262      0.104730140      0.076135404      0.063851193      0.081077913 
Barbadian          Yoruban             Esan            Mende          Gambian 
0.046363529      0.023775812      0.019882048      0.022351517      0.040591845 
Luhya         Japanese      Han.Chinese 
0.028636616      0.008061410      0.007995570 

stdev5 <- with(tbl2, tapply(tbl2$V5, tbl2$Pop, FUN = sd))  
stdev5
Finnish     N/W.European          Britain          Iberian           Tuscan 
0.005087753      0.005712915      0.004612070      0.012956475      0.010203078 
Aymara.Quechua         Peruvian          Mexican        Colombian       Portorican 
0.004921506      0.028152645      0.020057516      0.043086663      0.047546589 
Brazilian    Sickle.Brazil        Sickle.US African.American  SW.Afr.American 
0.098883846      0.058837415      0.076162989      0.072170469      0.079531621 
Barbadian          Yoruban             Esan            Mende          Gambian 
0.041403631      0.023199515      0.019790116      0.021577562      0.043602159 
Luhya         Japanese      Han.Chinese 
0.028304101      0.007744129      0.008104284 

stdev6 <- with(tbl2, tapply(tbl2$V6, tbl2$Pop, FUN = sd))  
stdev6
Finnish     N/W.European          Britain          Iberian           Tuscan 
0.078032186      0.033736967      0.024127563      0.020265744      0.019249674 
Aymara.Quechua         Peruvian          Mexican        Colombian       Portorican 
0.013114174      0.106485474      0.154104328      0.116420446      0.090942373 
Brazilian    Sickle.Brazil        Sickle.US African.American  SW.Afr.American 
0.250088266      0.120924288      0.065188609      0.050173465      0.059288297 
Barbadian          Yoruban             Esan            Mende          Gambian 
0.047953729      0.005350253      0.003502837      0.004762134      0.013611784 
Luhya         Japanese      Han.Chinese 
0.010590994      0.009824241      0.012381038 
