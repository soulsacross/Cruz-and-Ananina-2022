#----------------- Haplotype diversity from vcf files using pegas pkg ------------------------
library(pegas) 
library(vcfR)

clm <- read.vcfR( "~/Dropbox/paperultima/Diversity/HaplotypeDiversity/CLM.vcf", verbose = FALSE )

# get vcf into DNAbin format required for the Haplotype Diversity function in pegas
clm.dnabin <- vcfR2DNAbin(clm, extract.haps =T)

# data directly to hap.div: it considers all genome as a halotype, thus the chance of 
# getting 2 different haplotypes is always 1.0 (Hd = 1.0)
#hap.div(clm.dnabin, TRUE) # [1] 1.000000e+00 1.504965e-07

# thankfully, there is a sliding window function that accepts DNAbin format 
# This will make 5-steps walking of a window of size 10 seg. sites and calculate Haplotype Diversity inside that window
div.values <- sw(clm.dnabin, width = 10, step = 5, FUN = hap.div, rowAverage = TRUE, quiet = FALSE)

## if 50 marker windows were used, almost all "haplotypes" unique (value near 1.0,  Hd=0.9880358 in CLM)
mean(div.values)  # whole genome: 0.7991145 (chr22 alone: 0.8257916)


rm(clm,clm.dnabin,div.values)

###=== same thing for the remaining populations ======
br <- read.vcfR( "~/Dropbox/paperultima/Diversity/HaplotypeDiversity/BR.vcf", verbose = FALSE )
pur <- read.vcfR( "~/Dropbox/paperultima/Diversity/HaplotypeDiversity/PUR.vcf", verbose = FALSE )
mxl <- read.vcfR( "~/Dropbox/paperultima/Diversity/HaplotypeDiversity/MXL.vcf", verbose = FALSE )
pel <- read.vcfR( "~/Dropbox/paperultima/Diversity/HaplotypeDiversity/PEL.vcf", verbose = FALSE )

br.dnabin <- vcfR2DNAbin(br, extract.haps =T)
pur.dnabin <- vcfR2DNAbin(pur, extract.haps =T)
mxl.dnabin <- vcfR2DNAbin(mxl, extract.haps =T)
pel.dnabin <- vcfR2DNAbin(pel, extract.haps =T)

div.values.br <- sw(br.dnabin, width = 10, step = 5, FUN = hap.div, rowAverage = TRUE, quiet = FALSE) 
div.values.pur <- sw(pur.dnabin, width = 10, step = 5, FUN = hap.div, rowAverage = TRUE, quiet = FALSE)
div.values.mxl <- sw(mxl.dnabin, width = 10, step = 5, FUN = hap.div, rowAverage = TRUE, quiet = FALSE)
div.values.pel <- sw(pel.dnabin, width = 10, step = 5, FUN = hap.div, rowAverage = TRUE, quiet = FALSE)

mean(div.values.br) # BR: 0.8371473     CLM: 0.7991145 (see above)
mean(div.values.pur)# PUR:0.8126395 
mean(div.values.mxl)# MXL:0.7769805 
mean(div.values.pel)# MXL:0.7256145 (DivStat after running ~7h: 0.737965 [not finished])

rm(br,br.dnabin,div.values.br)
rm(pur,pur.dnabin,div.values.pur)
rm(mxl,mxl.dnabin,div.values.mxl)
rm(pel,pel.dnabin,div.values.pel)

