####=== Here I am using normalized iHS values from selscan to do a 
##  ---- circular plot ----------
# First, I will use by RCircos to draw circular plot, then SOFIA
# see https://cran.r-project.org/web/packages/RCircos/
## Pedro 2020.05.01

library(RCircos)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(Homo.sapiens)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(VariantAnnotation)
library(annotate) 

setwd("~/Dropbox/input_chr16/haplotype.scores/")

# Use dplyr/tidyr to build data frame with chr information from filename
# Downside is that you'll have to feed RCircos with data as.data.frame() function, otherwise it crashes
# See https://www.biostars.org/p/251987/
# Check if 'dim(MXL)' matches 'cat *.MXL.ihs.out.ihs.out.100bins.norm | wc -l' in terminal

## Use *ihs.out.100bins.norm normalized iHS files from selscan as input
# Mexicans
Filenames <- dir(pattern = "*.MXL.ihs.out.ihs.out.100bins.norm")
MXL = tibble(File = Filenames) %>%
  extract(File, "Chromosome", "([^.]+)", remove = FALSE) %>%  # Regex to get everything but the period  
  mutate(Data = lapply(File, read.table)) %>%
  unnest(Data) %>%
  dplyr::select(-File)
colnames(MXL) <- c("Chromosome", "chromStart", "chromEnd", "1_freq", "ihh_1", "ihh_0", "ihs.unnorm", "ihs", "signal")

# Peruvians
Filenames <- dir(pattern = "*.PEL.ihs.out.ihs.out.100bins.norm")
PEL = tibble(File = Filenames) %>%
  extract(File, "Chromosome", "([^.]+)", remove = FALSE) %>%
  mutate(Data = lapply(File, read.table)) %>%
  unnest(Data) %>%
  dplyr::select(-File)
colnames(PEL) <- c("Chromosome", "chromStart", "chromEnd", "1_freq", "ihh_1", "ihh_0", "ihs.unnorm", "ihs", "signal")

# Colombians
Filenames <- dir(pattern = "*.CLM.ihs.out.ihs.out.100bins.norm")
CLM = tibble(File = Filenames) %>%
  extract(File, "Chromosome", "([^.]+)", remove = FALSE) %>%
  mutate(Data = lapply(File, read.table)) %>%
  unnest(Data) %>%
  dplyr::select(-File)
colnames(CLM) <- c("Chromosome", "chromStart", "chromEnd", "1_freq", "ihh_1", "ihh_0", "ihs.unnorm", "ihs", "signal")

# Puerto Ricans
Filenames <- dir(pattern = "*.PUR.ihs.out.ihs.out.100bins.norm")
PUR = tibble(File = Filenames) %>%
  extract(File, "Chromosome", "([^.]+)", remove = FALSE) %>%
  mutate(Data = lapply(File, read.table)) %>%
  unnest(Data) %>%
  dplyr::select(-File)
colnames(PUR) <- c("Chromosome", "chromStart", "chromEnd", "1_freq", "ihh_1", "ihh_0", "ihs.unnorm", "ihs", "signal")

# Brazilians
Filenames <- dir(pattern = "*.BR.ihs.out.ihs.out.100bins.norm")
BR = tibble(File = Filenames) %>%
  extract(File, "Chromosome", "([^.]+)", remove = FALSE) %>%
  mutate(Data = lapply(File, read.table)) %>%
  unnest(Data) %>%
  dplyr::select(-File)
colnames(BR) <- c("Chromosome", "chromStart", "chromEnd", "1_freq", "ihh_1", "ihh_0", "ihs.unnorm", "ihs", "signal")

#######============ 1. Further filtering and analyses ===========
## Play around with custom data
Mx <- as.data.frame(MXL[which(MXL$signal == "1"), ]) # 1 means suggestive selection ##9554 variants
# Make a GenomicRanges list 
grMx <- makeGRangesFromDataFrame(Mx, keep.extra.columns=TRUE) 
# # get Entrez gene ID
# genesMx <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), grMx)
# # Then get gene SYMBOL
# gene.sb.Mx <- getSYMBOL(genesMx$gene_id, data='org.Hs.eg')
# # Write gene list to file for further analysis
# gMXL <-as.data.frame(gene.sb.Mx, data = "org.Hs.eg")
# write.table(gMXL, "GenesUnderSelectionMXL.txt",col.names = T, row.names = F, quote = F)

# Same as above for other Latinos: PEL
Pl <- as.data.frame(PEL[which(PEL$signal == "1"), ]) # 1 means suggestive selection ##8606 variants
# Make a GenomicRanges list 
grPl<- makeGRangesFromDataFrame(Pl, keep.extra.columns=TRUE) 
# # get Entrez gene ID
# genesPl <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), grPl)
# # Then get SYMBOL
# gene.sb.Pl <- getSYMBOL(genesPl$gene_id, data='org.Hs.eg')
# # Write gene list to file for further analysis
# gPEL <-as.data.frame(gene.sb.Pl, data = "org.Hs.eg")
# write.table(gPEL, "GenesUnderSelectionPEL.txt",col.names = T, row.names = F, quote = F)

# Same as above for other Latinos: CLM
Cl <- as.data.frame(CLM[which(CLM$signal == "1"), ]) # 1 means suggestive selection ##9804 variants
# Make a GenomicRanges list 
grCl<- makeGRangesFromDataFrame(Cl, keep.extra.columns=TRUE) 
# # get Entrez gene ID
# genesCl <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), grCl)
# # Then get SYMBOL
# gene.sb.Cl <- getSYMBOL(genesCl$gene_id, data='org.Hs.eg')
# # Write gene list to file for further analysis
# gCLM <-as.data.frame(gene.sb.Cl, data = "org.Hs.eg")
# write.table(gCLM, "GenesUnderSelectionCLM.txt",col.names = T, row.names = F, quote = F)

# Same as above for other Latinos: PUR
Pr <- as.data.frame(PUR[which(PUR$signal == "1"), ]) # 1 means suggestive selection ##10186 variants
# Make a GenomicRanges list 
grPr<- makeGRangesFromDataFrame(Pr, keep.extra.columns=TRUE) 
# # get Entrez gene ID
# genesPr <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), grPr)
# # Then get SYMBOL
# gene.sb.Pr <- getSYMBOL(genesPr$gene_id, data='org.Hs.eg')
# # Write gene list to file for further analysis
# gPUR <-as.data.frame(gene.sb.Pr, data = "org.Hs.eg")
# write.table(gPUR, "GenesUnderSelectionPUR.txt",col.names = T, row.names = F, quote = F)

# Same as above for other Latinos: BR
Br <- as.data.frame(BR[which(BR$signal == "1"), ]) # 1 means suggestive selection ##11066 variants
# Make a GenomicRanges list 
grBr<- makeGRangesFromDataFrame(Br, keep.extra.columns=TRUE) 
# # get Entrez gene ID
# genesBr <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), grBr)
# # Then get SYMBOL
# gene.sb.Br <- getSYMBOL(genesBr$gene_id, data='org.Hs.eg')
# # Write gene list to further analysis
# gBR <-as.data.frame(gene.sb.Br, data = "org.Hs.eg")
# write.table(gBR, "GenesUnderSelectionBR.txt",col.names = T, row.names = F, quote = F)

###== A lot of markers and different genes by population, poorly informative ====
##= Let's filter variants in coding region 
# needs library(VariantAnnotation)
vcf <- readVcf("/Volumes/PEDRO 70GB/Genetics/ALLinfo/VCFinfo/ALL.vcf.gz", "hg19")
rd <- rowRanges(vcf)
rm(vcf)   # too large to keep
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# change notation to chr#
library(diffloop)
rd <- addchr(rd)

# Here is the GRanges list with REF/ALT information
mergedMx <- merge(grMx, rd)
mergedPl <- merge(grPl, rd)
mergedCl <- merge(grCl, rd)
mergedPr <- merge(grPr, rd)
mergedBr <- merge(grBr, rd)

# # If interested in the functional location of markers - coding
# loc <- locateVariants(mergedMx, txdb, CodingVariants())
# # Same for all types of region
# allvar <- locateVariants(rd, txdb, AllVariants())

# Inspect for synonymous/nonsynonymous (stored in CONSEQUENCE column)
codingMx <- VariantAnnotation::predictCoding(mergedMx, txdb, seqSource=Hsapiens, varAllele=mergedMx@elementMetadata$REF)
codingPl <- VariantAnnotation::predictCoding(mergedPl, txdb, seqSource=Hsapiens, varAllele=mergedPl@elementMetadata$REF)
codingCl <- VariantAnnotation::predictCoding(mergedCl, txdb, seqSource=Hsapiens, varAllele=mergedCl@elementMetadata$REF)
codingPr <- VariantAnnotation::predictCoding(mergedPr, txdb, seqSource=Hsapiens, varAllele=mergedPr@elementMetadata$REF)
codingBr <- VariantAnnotation::predictCoding(mergedBr, txdb, seqSource=Hsapiens, varAllele=mergedBr@elementMetadata$REF)

##-- Getting genes
# get Entrez gene ID
genesMx <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), codingMx)
# Then get gene SYMBOL
gene.sb.Mx <- getSYMBOL(genesMx$gene_id, data='org.Hs.eg')
# Write gene list to file for further analysis
gMXL <-as.data.frame(gene.sb.Mx, data = "org.Hs.eg")
write.table(gMXL, "GenesUnderSelectionMXL.txt",col.names = T, row.names = F, quote = F)

# get Entrez gene ID
genesPl <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), codingPl)
# Then get SYMBOL
gene.sb.Pl <- getSYMBOL(genesPl$gene_id, data='org.Hs.eg')
# Write gene list to file for further analysis
gPEL <-as.data.frame(gene.sb.Pl, data = "org.Hs.eg")
write.table(gPEL, "GenesUnderSelectionPEL.txt",col.names = T, row.names = F, quote = F)

# get Entrez gene ID
genesCl <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), codingCl)
# Then get SYMBOL
gene.sb.Cl <- getSYMBOL(genesCl$gene_id, data='org.Hs.eg')
# Write gene list to file for further analysis
gCLM <-as.data.frame(gene.sb.Cl, data = "org.Hs.eg")
write.table(gCLM, "GenesUnderSelectionCLM.txt",col.names = T, row.names = F, quote = F)

# get Entrez gene ID
genesPr <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), codingPr)
# Then get SYMBOL
gene.sb.Pr <- getSYMBOL(genesPr$gene_id, data='org.Hs.eg')
# Write gene list to file for further analysis
gPUR <-as.data.frame(gene.sb.Pr, data = "org.Hs.eg")
write.table(gPUR, "GenesUnderSelectionPUR.txt",col.names = T, row.names = F, quote = F)

# get Entrez gene ID
genesBr <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), codingBr)
# Then get SYMBOL
gene.sb.Br <- getSYMBOL(genesBr$gene_id, data='org.Hs.eg')
# Write gene list to further analysis
gBR <-as.data.frame(gene.sb.Br, data = "org.Hs.eg")
write.table(gBR, "GenesUnderSelectionBR.txt",col.names = T, row.names = F, quote = F)

# Save as csv                                                                                                         All coding |25 perc. iHS && coding
write.table(x = data.frame(codingMx), file = "codingMXL.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE ) # 255 rows |95 rows # write.table(x = data.frame(codingMx), file = "~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/codingMXL.xls", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table(x = data.frame(codingPl), file = "codingPEL.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE ) # 235 rows |69 
write.table(x = data.frame(codingCl), file = "codingCLM.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE ) # 203 rows |52 
write.table(x = data.frame(codingPr), file = "codingPUR.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE ) # 265 rows |64 
write.table(x = data.frame(codingBr), file = "codingBR.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE ) #  308 rows |68


# Make a df out of each granges object
codingMx_df <- data.frame(codingMx,  col.names=TRUE, row.names=1:255)
codingPl_df <- data.frame(codingPl,  col.names=TRUE, row.names=1:235)
codingCl_df <- data.frame(codingCl,  col.names=TRUE, row.names=1:203)
codingPr_df <- data.frame(codingPr,  col.names=TRUE, row.names=1:265)
codingBr_df <- data.frame(codingBr,  col.names=TRUE, row.names=1:308)

# Add gene symbol to dataframe
codingMx_df$SYMBOL<- gMXL[match(codingMx_df$GENEID, row.names(gMXL)),]
codingPl_df$SYMBOL<- gPEL[match(codingPl_df$GENEID, row.names(gPEL)),]
codingCl_df$SYMBOL<- gCLM[match(codingCl_df$GENEID, row.names(gCLM)),]
codingPr_df$SYMBOL<- gPUR[match(codingPr_df$GENEID, row.names(gPUR)),]
codingBr_df$SYMBOL<- gBR[match(codingBr_df$GENEID, row.names(gBR)),]

# Keep also discrepant values (percentil)
codingBr_df < 
Mx <- tmp[tmp$ihs <= quantile(tmp$ihs, 0.125) | tmp$ihs >= quantile(tmp$ihs, 0.875), ] 
Pl <- tmp[tmp$ihs <= quantile(tmp$ihs, 0.125) | tmp$ihs >= quantile(tmp$ihs, 0.875), ] 
Cl <- tmp[tmp$ihs <= quantile(tmp$ihs, 0.125) | tmp$ihs >= quantile(tmp$ihs, 0.875), ] 
Pr <- tmp[tmp$ihs <= quantile(tmp$ihs, 0.125) | tmp$ihs >= quantile(tmp$ihs, 0.875), ] 
Br <- tmp[tmp$ihs <= quantile(tmp$ihs, 0.125) | tmp$ihs >= quantile(tmp$ihs, 0.875), ] 

  
# Combine all genes to build our gene track in RCircos circular plot
# Let's filter iHS do build different gene tracks while keeping important columns
tmp  <- rbind(
  na.omit(codingMx_df[codingMx$ihs<= quantile(codingMx$ihs, 0.125) | codingMx$ihs>= quantile(codingMx$ihs, 0.875), c(1,2,3,32,26)]), 
  na.omit(codingPl_df[codingPl$ihs<= quantile(codingPl$ihs, 0.125) | codingPl$ihs>= quantile(codingPl$ihs, 0.875),c(1,2,3,32,26)]), 
  na.omit(codingCl_df[codingCl$ihs<= quantile(codingCl$ihs, 0.125) | codingCl$ihs>= quantile(codingCl$ihs, 0.875),c(1,2,3,32,26)]), 
  na.omit(codingPr_df[codingPr$ihs<= quantile(codingPr$ihs, 0.125) | codingPr$ihs>= quantile(codingPr$ihs, 0.875),c(1,2,3,32,26)]), 
  na.omit(codingBr_df[codingBr$ihs<= quantile(codingBr$ihs, 0.125) | codingBr$ihs>= quantile(codingBr$ihs, 0.875),c(1,2,3,32,26)])) 
gene.track25 <- distinct(tmp) # dplyr to supress duplicates 

# Take a look on nonsynonymous
# gene.track25[gene.track25$CONSEQUENCE=='nonsynonymous',] #26  ### 286 (ALL iHS)

# Now for the remaining percentil
tmp  <- rbind(
  na.omit(codingMx_df[codingMx$ihs>= quantile(codingMx$ihs, 0.125) | codingMx$ihs<= quantile(codingMx$ihs, 0.875), c(1,2,3,32,26)]), 
  na.omit(codingPl_df[codingPl$ihs>= quantile(codingPl$ihs, 0.125) | codingPl$ihs<= quantile(codingPl$ihs, 0.875),c(1,2,3,32,26)]), 
  na.omit(codingCl_df[codingCl$ihs>= quantile(codingCl$ihs, 0.125) | codingCl$ihs<= quantile(codingCl$ihs, 0.875),c(1,2,3,32,26)]), 
  na.omit(codingPr_df[codingPr$ihs>= quantile(codingPr$ihs, 0.125) | codingPr$ihs<= quantile(codingPr$ihs, 0.875),c(1,2,3,32,26)]), 
  na.omit(codingBr_df[codingBr$ihs>= quantile(codingBr$ihs, 0.125) | codingBr$ihs<= quantile(codingBr$ihs, 0.875),c(1,2,3,32,26)])) 
gene.track75 <- distinct(tmp)

# Take a look on nonsynonymous
# gene.track75[gene.track75$CONSEQUENCE=='nonsynonymous',] #97  #### 286

####================== 2. Plot data with RCircos (not so pleasing aestetically) ==============================
# Initialize RCircos core components
data(UCSC.HG19.Human.CytoBandIdeogram)
chr.exclude <- c("chrX", "chrY") # chr.exclude <- c("chrX", "chrY", "chr22", "chr21", "chr20", "chr19", "chr18", "chr17", "chr16", "chr15", "chr14", "chr13", "chr12", "chr11")
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
tracks.inside <- 12
tracks.outside <- 2

# IF modifying core components use this block
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$base.per.unit <- 15000
RCircos.Reset.Plot.Parameters(rcircos.params)
RCircos.List.Plot.Parameters()


RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

# Initialize Graphic Device
out.file <- "RCircosDemoHumanGenome_10.pdf"
pdf(file=out.file, height=8, width=8, compress=T)
RCircos.Set.Plot.Area()
title("Extended haplotype score for Latinos")

# Plot Chromosome Ideogram
RCircos.Chromosome.Ideogram.Plot()
RCircos.Scatter.Plot(codingMx_df, data.col = 8, track.num = 5, side = "in")
RCircos.Scatter.Plot(codingPl_df, data.col = 8, track.num = 6, side = "in")
RCircos.Scatter.Plot(codingCl_df, data.col = 8, track.num = 7, side = "in")
RCircos.Scatter.Plot(codingPr_df, data.col = 8, track.num = 8, side = "in")
RCircos.Scatter.Plot(codingBr_df, data.col = 8, track.num = 9, side = "in")

# Add genes of the iHS 25% percentil
name.col <- 4
side <- "in"
track.num <- 1
RCircos.Gene.Connector.Plot(gene.track25[gene.track25$CONSEQUENCE=='nonsynonymous',], track.num, side);

track.num <- 2;
RCircos.Gene.Name.Plot(gene.track25[gene.track25$CONSEQUENCE=='nonsynonymous',], name.col,track.num, side);

# Add genes the remaining genes
name.col <- 4
side <- "out"
track.num <- 1
RCircos.Gene.Connector.Plot(gene.track75[gene.track75$CONSEQUENCE=='nonsynonymous',], track.num, side);

track.num <- 2;
RCircos.Gene.Name.Plot(gene.track75[gene.track75$CONSEQUENCE=='nonsynonymous',], name.col,track.num, side);

#
dev.off()

# # Toy example
# data(RCircos.Line.Data); 
# data.col <- 5;
# track.num <- 7; 
# side <- "in"; 
# RCircos.Line.Plot(RCircos.Line.Data, data.col = 5, track.num, side);

mark0<-UCSC.HG19.Human.CytoBandIdeogram

data<-data.frame(locus=mark0$Band,pos=mark0$ChromStart/1000000,chr=mark0$Chromosome,map=rep(1,length(mark0$Chromosome)))
data$locus<-as.character(data$locus)
data$chr<-as.character(data$chr)
data$text<-NA
data$text[sample(1:length(data$locus),150)]<-as.character(RCircos.Gene.Label.Data$Gene[sample(1:length(RCircos.Gene.Label.Data$Gene),150)])



data$kar<-data$pos

data$heatmap<-RCircos.Heatmap.Data[sample(1:nrow(RCircos.Heatmap.Data),nrow(mark0)),5]
data$scatter<-RCircos.Scatter.Data[sample(1:nrow(RCircos.Scatter.Data),nrow(mark0)),5]
data$line<-RCircos.Line.Data[sample(1:nrow(RCircos.Line.Data),nrow(mark0)),5]


tmp<-numeric()
while(length(tmp)<nrow(data)){
  tmp<-c(tmp,rep(runif(1),sample(1:20,1)))
}
tmp<-tmp[1:nrow(data)]
data$histogram<-tmp


data2<-data[1:200,]
data2[,5:ncol(data2)]<-NA
data2$map<-2
data3<-rbind(data,data2)


####=============== 3. Plot data with circos through SOFIA ==============================
library(SOFIA)
library(biomaRt)
# From supplementary data from SOFIA: https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/jhered/108/4/10.1093_jhered_esx023/2/esx023_suppl_Supplementary_File_1.docx?Expires=1602384109&Signature=y8dTUGKxAycqAEQMREf0PI9ME~fNq7moQqhn2RRsMnG~5-JObthlCKR~skuFbIP1urqfAxrS-rYp6VVJoQHBkGntgf-NZlipI1a888N80r7PrrSeZYXJPJ7fyEN2j4yUpb6uXz9boGl9Xc6Pt8cicQTqguEeTsy~w~FJX5tUnWicRvxkjWpZg0augfnpOn8HLzc2twqZBayPMcEMc~GKrzYd5qxIJJFHQ9f3jG-YsEgHopMuniw~6LSD9jkz~fuDynG6Zl8ZdEjPlX~M6un5W6Jor9lbTQkoMy3CLInBmTJXtamddHBTc3ia5HWGzNFbspc5eOEBsaamPUf0QjgymQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA
# In this topic I will also try to filter by percentile alone (to undercut to ~88-90 genes).
# These 88 genes are equivalent to 1-percentil of total genes (8893) found to be associated with singals (I have used RGmatch). See cms.sh and iHS.all.xlsx

# Load RGmatch-generated list of associated genes. Genetate this output like the following:
# python ~/Downloads/pfurio-rgmatch-e9289746a5bd/rgmatch.py  \
# -r gene \
# -g /Volumes/PEDRO\ 70GB/Genetics/GTF/hg19/gencode.v29lift37.annotation.gtf.gz \
# -b ~/Dropbox/input_chr16/haplotype.scores/iHS/CLM.xls \
# -o ~/Dropbox/input_chr16/haplotype.scores/iHS/CLM.genomic.regions.gene.xls
Mx <- read.table("./iHS/MXL.genomic.regions.gene.xls", header = T)
  Mx$Gene <- gsub("\\..*","",Mx$Gene)
Pl <- read.table("./iHS/PEL.genomic.regions.gene.xls", header = T)
  Pl$Gene <- gsub("\\..*","",Pl$Gene)
Cl <- read.table("./iHS/CLM.genomic.regions.gene.xls", header = T)
  Cl$Gene <- gsub("\\..*","",Cl$Gene)
Pr <- read.table("./iHS/PUR.genomic.regions.gene.xls", header = T)
  Pr$Gene <- gsub("\\..*","",Pr$Gene)
Br <- read.table("./iHS/BR.genomic.regions.gene.xls", header = T)
  Br$Gene <- gsub("\\..*","",Br$Gene)

# Because I am currently not working with a GRanges object, lets use biomart
# It has a very comprehensive db of biologic data. See in 'listAttributes(mart)'.   
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- c(Mx$Gene, Pl$Gene,Cl$Gene,Pr$Gene,Br$Gene) 
#Mx$id <- NA
Genelist <- getBM(filters= "ensembl_gene_id", 
               attributes= c("ensembl_gene_id", "hgnc_symbol"),
               values=genes,mart= mart)

# Add gene symbol to dataframe
Mx$SYMBOL<- Genelist[match(Mx$Gene, Genelist$ensembl_gene_id),2]
Pl$SYMBOL<- Genelist[match(Pl$Gene, Genelist$ensembl_gene_id),2]
Cl$SYMBOL<- Genelist[match(Cl$Gene, Genelist$ensembl_gene_id),2]
Pr$SYMBOL<- Genelist[match(Pr$Gene, Genelist$ensembl_gene_id),2]
Br$SYMBOL<- Genelist[match(Br$Gene, Genelist$ensembl_gene_id),2]

length(unique(Genelist$hgnc_symbol)) # 6467 (out of 8893 RGmatch generated ENSGs) 

##=== Start actual preparation to plot
# Combine all genes to build our gene track for circular plot
# Let's filter iHS do build different gene tracks while keeping important columns - 1 percentile                                                Count  
extreme.iHS.Mx <- na.omit(Mx[Mx$thickEnd<=quantile(Mx$thickEnd, 0.005) | Mx$thickEnd>=quantile(Mx$thickEnd, 0.995),c(1,2,17,15)]) # 82 write.table(extreme.iHS.Mx, file="~/Dropbox/input_chr16/haplotype.scores/iHS/extreme.iHS.Mx.txt",  row.names = F, col.names = T, quote = F )
extreme.iHS.Pl <- na.omit(Pl[Pl$thickEnd<=quantile(Pl$thickEnd, 0.005) | Pl$thickEnd>=quantile(Pl$thickEnd, 0.995),c(1,2,17,15)]) # 75
extreme.iHS.Cl <- na.omit(Cl[Cl$thickEnd<=quantile(Cl$thickEnd, 0.005) | Cl$thickEnd>=quantile(Cl$thickEnd, 0.995),c(1,2,17,15)]) # 83
extreme.iHS.Pr <- na.omit(Pr[Pr$thickEnd<=quantile(Pr$thickEnd, 0.005) | Pr$thickEnd>=quantile(Pr$thickEnd, 0.995),c(1,2,17,15)]) # 85
extreme.iHS.Br <- na.omit(Br[Br$thickEnd<=quantile(Br$thickEnd, 0.005) | Br$thickEnd>=quantile(Br$thickEnd, 0.995),c(1,2,17,15)]) # 92
#write.table(rbind(extreme.iHS.Mx,extreme.iHS.Pl,extreme.iHS.Cl,extreme.iHS.Pr,extreme.iHS.Br), file = "./iHS/1percentil.extreme.all.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE ) #  308 rows |68
# in bash: cat /Users/Pedro/Dropbox/input_chr16/haplotype.scores/iHS/1percentil.extreme.all.txt | awk '{print $2}' | grep "[A-Z][A-Z]*" | wc -l #287 (actually 204 genes)

iHS.to.plot<- rbind(extreme.iHS.Mx, extreme.iHS.Pl,extreme.iHS.Cl, extreme.iHS.Pr, extreme.iHS.Br)

# To get chr only: gsub("[_][0-9]*", "", extreme.iHS.Br$Region)
iHS.gen<-data.frame(map=1,chr=gsub("[_][0-9]*", "", iHS.to.plot$Region),pos=iHS.to.plot$Midpoint,locus=iHS.to.plot$SYMBOL, text=iHS.to.plot$SYMBOL)

# Set iHS columns for each pop
iHS.gen$MXL<-NA 
iHS.gen$PEL<-NA
iHS.gen$CLM<-NA
iHS.gen$PUR<-NA
iHS.gen$BR<-NA

# Populate with signals from Latinos
iHS.gen$MXL[1:82] <- extreme.iHS.Mx$thickEnd
iHS.gen$PEL[83:157] <- extreme.iHS.Pl$thickEnd
iHS.gen$CLM[158:240] <- extreme.iHS.Cl$thickEnd
iHS.gen$PUR[241:325] <- extreme.iHS.Pr$thickEnd
iHS.gen$BR[326:417] <- extreme.iHS.Br$thickEnd

# NAs in empty symbols (actually in whole df, in this particular case only $symbols has empty lines) and "duplicated" genes (i.e, occurring in > 1 pop)
iHS.gen[iHS.gen==""]<-NA
iHS.gen$text[duplicated(iHS.gen$text)] <- NA

#order by chr and pos
chrOrder <-paste(rep("chr",22),1:22, sep="")
iHS.gen$chr <- factor(iHS.gen$chr, chrOrder, ordered=TRUE)
iHS.gen<-iHS.gen[do.call(order, iHS.gen[, c("map", "chr", "pos")]), ]

# obtuse conversion to gen. distance
iHS.gen$pos <- iHS.gen$pos/1000000

# keep only first occurrence of a gene symbol
iHS.gen$text <- ave(
  iHS.gen$text, 
  iHS.gen$locus, 
  FUN = function(a) replace(a, duplicated(a), NA_integer_)
)

###=====================
### 2022 update, now with markers that have overcome simulation thresholds 
#####=============
iHS.to.plot<- read.table("~/Dropbox/paperultima/iHS\ validation/iHS.surpassed.genomic.regions.genic.uniq.xls", header = T, sep = "\t")

iHS.gen<-data.frame(map=1,chr=gsub("[_][0-9]*", "", iHS.to.plot$Region),pos=iHS.to.plot$Midpoint,locus=iHS.to.plot$Gene, text=iHS.to.plot$Gene, iHS=iHS.to.plot$thickEnd, pop=iHS.to.plot$blockCount)

iHS.gen$MXL<-NA 
iHS.gen$PEL<-NA
iHS.gen$CLM<-NA
iHS.gen$PUR<-NA
iHS.gen$BR<-NA

# Populate with signals from Latinos (this is a dumb way to do it, refine that)
iHS.gen$BR[1:2] <- iHS.gen$iHS[iHS.gen$pop=="BR"]
iHS.gen$CLM[3:22] <- iHS.gen$iHS[iHS.gen$pop=="CLM"]
iHS.gen$MXL[23:80] <- iHS.gen$iHS[iHS.gen$pop=="MXL"]
iHS.gen$PEL[81:163] <- iHS.gen$iHS[iHS.gen$pop=="PEL"]
iHS.gen$PUR[164:173] <- iHS.gen$iHS[iHS.gen$pop=="PUR"]

iHS.gen$pop<- NULL
iHS.gen$iHS<-NULL
# NAs in empty symbols (actually in whole df, in this particular case only $symbols has empty lines) and "duplicated" genes (i.e, occurring in > 1 pop)
# iHS.gen[iHS.gen==""]<-NA
# iHS.gen$text[duplicated(iHS.gen$text)] <- NA

#order by chr and pos
chrOrder <-paste(rep("chr",22),1:22, sep="")
iHS.gen$chr <- factor(iHS.gen$chr, chrOrder, ordered=TRUE)
iHS.gen<-iHS.gen[do.call(order, iHS.gen[, c("map", "chr", "pos")]), ]

# obtuse conversion to gen. distance
iHS.gen$pos <- iHS.gen$pos/1000000

# keep only first occurrence of a gene symbol
iHS.gen$text <- ave(
  iHS.gen$text, 
  iHS.gen$locus, 
  FUN = function(a) replace(a, duplicated(a), NA_integer_)
)

# To test if the above manipulations meet SOFIA format requirements to generate plot do: (thus rule out input format problems if the cmds below fail to generate a plot). This is because SOFIA run is poorly informative and error msg can be cryptical
#SOFIA(data=iHS.gen, chromoConfiguration=NULL, circosLocation='~/Dropbox/circos-0.69-9')

# Set plot aesthetics
plotLocation<-data.frame(r0=c(0.90,.81,.71,.61,.51,.41),r1=c(.99,.9,.8,.7,.6,.5))
plotBackground<-data.frame(backgroundShow=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE),
                           backgroundColor=rep('vvlgrey',6),axisShow=rep(TRUE,6),axisSep=rep(4,6))                          
chromoConfiguration<-data.frame(paste(rep("chr",22),1:22, sep=""), map=rep(1,22),
                                rev=c(rep(FALSE,22)),color=c(rep(c('vvdblue_a3','lblue_a3'),11)),
                                radius=rep(1,22))
plotType<-c('text',rep('scatter',5)) 

# to explore colors: library(RColorBrewer) then display.brewer.all(colorblindFriendly=TRUE) and brewer.pal.info
# original colors were piyg-11-div for all 5 scatter plots
plotColor<-c('black','green','green','green','green','green')

# Defining marker size. 
# For the text plot, the marker size defines the font size. 
# For scatter, it defines the circle size while for line plots, it defines the line thickness. 
# For heatmaps, any random number can be included since it is not used at all. 
markerSize<-c(16,16,16,16,16,16)

# to plot background without data open /Users/Pedro/Dropbox/circos-0.69-9/bin/circos.conf and delete
# everything inside 2 "<backgrounds>" (5x because there are 5 scatterplots) and run circos on terminal:
# ~/Dropbox/circos-0.69-9/bin/circos -conf circos.conf
SOFIA(data=iHS.gen,
      linkColor='chr',
      linksFlag=F,
      chromoConfiguration=chromoConfiguration,
      plotBackground=plotBackground,
      chrPrefixFont='lower',
      plotLocation=plotLocation,
      plotType=plotType,
      plotColor=plotColor
      ,markerSize=markerSize,
      circosLocation='~/Dropbox/circos-0.69-9',
      tickSuffix='cM',
      returnConf=TRUE,
      circosDisplay=TRUE)

#######======== Stuff for making supplementary tables to pulish =================#########

##== 1) To write iHS table of 1percentile markers============
iHS.gen$pos <- iHS.gen$pos*1000000
write.table(iHS.gen, "1percentilebyPop.txt", col.names = T, row.names = F, quote = F, sep="\t")

##== 2) To write iHS table of 1percentile markers============
# Print out to run RGmatch
iHS.codingMXL <- read.csv("/Users/Pedro/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/codingMXL.csv")
iHS.codingPEL <- read.csv("/Users/Pedro/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/codingPEL.csv")
iHS.codingCLM <- read.csv("/Users/Pedro/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/codingCLM.csv")
iHS.codingPUR <- read.csv("/Users/Pedro/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/codingPUR.csv")
iHS.codingBR <- read.csv("/Users/Pedro/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/codingBR.csv")

# RGmatch cmd: python  ~/Downloads/pfurio-rgmatch-e9289746a5bd/rgmatch.py -r gene  -g /Volumes/PEDRO\ 70GB/Genetics/GTF/hg19/gencode.v29lift37.annotation.gtf.gz  -b ~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/codingPEL.xls -o  ~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/PEL.rgmatch.txt 
# Then
write.table(x =iHS.codingMXL , file = "~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/codingMXL.xls", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table(x =iHS.codingPEL , file = "~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/codingPEL.xls", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table(x =iHS.codingCLM , file = "~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/codingCLM.xls", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table(x =iHS.codingPUR , file = "~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/codingPUR.xls", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table(x =iHS.codingBR , file = "~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/codingBR.xls", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )

# The rest is basically repeating the commands to produce a SOFIA plot dataframe (copied from above)
Mx <- read.table("~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/MXL.rgmatch.txt", header = T)
Mx$Gene <- gsub("\\..*","",Mx$Gene)
Mx<- Mx[which(Mx$Distance == 0),]  #286 rows
Pl <- read.table("~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/PEL.rgmatch.txt", header = T)
Pl$Gene <- gsub("\\..*","",Pl$Gene)
Pl<- Pl[which(Pl$Distance == 0),]   #272 rows
Cl <- read.table("~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/CLM.rgmatch.txt", header = T)
Cl$Gene <- gsub("\\..*","",Cl$Gene)
Cl<- Cl[which(Cl$Distance == 0),]   #243 rows
Pr <- read.table("~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/PUR.rgmatch.txt", header = T)
Pr$Gene <- gsub("\\..*","",Pr$Gene)
Pr<- Pr[which(Pr$Distance == 0),]   #274
Br <- read.table("~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/BR.rgmatch.txt", header = T)
Br$Gene <- gsub("\\..*","",Br$Gene)
Br<- Br[which(Br$Distance == 0),]   #308 rows

# Because I am currently not working with a GRanges object, lets use biomart
# It has a very comprehensive db of biologic data. See in 'listAttributes(mart)'.   
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- c(Mx$Gene, Pl$Gene,Cl$Gene,Pr$Gene,Br$Gene) 
#Mx$id <- NA
Genelist <- getBM(filters= "ensembl_gene_id", 
                  attributes= c("ensembl_gene_id", "hgnc_symbol"),
                  values=genes,mart= mart)

# Add gene symbol to dataframe
Mx$SYMBOL<- Genelist[match(Mx$Gene, Genelist$ensembl_gene_id),2]
Pl$SYMBOL<- Genelist[match(Pl$Gene, Genelist$ensembl_gene_id),2]
Cl$SYMBOL<- Genelist[match(Cl$Gene, Genelist$ensembl_gene_id),2]
Pr$SYMBOL<- Genelist[match(Pr$Gene, Genelist$ensembl_gene_id),2]
Br$SYMBOL<- Genelist[match(Br$Gene, Genelist$ensembl_gene_id),2]

length(unique(Genelist$hgnc_symbol)) # 241 (out of 8893 RGmatch generated ENSGs) 

iHS.to.plot <- rbind(Mx, Pl, Cl, Pr, Br)

# To get chr only: gsub("[_][0-9]*", "", extreme.iHS.Br$Region)
iHS.gen<-data.frame(map=1,chr=gsub("[_][0-9]*", "", iHS.to.plot$Region),pos=iHS.to.plot$Midpoint,locus=iHS.to.plot$SYMBOL, text=iHS.to.plot$SYMBOL)

# Set iHS columns for each pop
iHS.gen$MXL<-NA 
iHS.gen$PEL<-NA
iHS.gen$CLM<-NA
iHS.gen$PUR<-NA
iHS.gen$BR<-NA

# Populate with signals from Latinos
iHS.gen$MXL[1:286] <- Mx$blockCount
iHS.gen$PEL[287:558] <- Pl$blockCount
iHS.gen$CLM[559:801] <- Cl$blockCount
iHS.gen$PUR[802:1075] <- Pr$blockCount
iHS.gen$BR[1076:1383] <- Br$blockCount

# NAs in empty symbols (actually in whole df, in this particular case only $symbols has empty lines) and "duplicated" genes (i.e, occurring in > 1 pop)
iHS.gen[iHS.gen==""]<-NA
iHS.gen$text[duplicated(iHS.gen$text)] <- NA

#order by chr and pos
chrOrder <-paste(rep("chr",22),1:22, sep="")
iHS.gen$chr <- factor(iHS.gen$chr, chrOrder, ordered=TRUE)
iHS.gen<-iHS.gen[do.call(order, iHS.gen[, c("map", "chr", "pos")]), ]

# obtuse conversion to gen. distance - uncomment to plot
# iHS.gen$pos <- iHS.gen$pos/1000000

# keep only first occurrence of a gene symbol
iHS.gen$text <- ave(
  iHS.gen$text, 
  iHS.gen$locus, 
  FUN = function(a) replace(a, duplicated(a), NA_integer_)
)

# Finally, remove duplicated rows:
iHS.gen <- unique(iHS.gen)

# Write the table
write.table(iHS.gen, "AllinCodingRegions.txt", col.names = T, row.names = F, quote = F, sep="\t")

