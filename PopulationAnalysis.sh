###===================== Calculate effective population size (Ne) ============================================
# for admixed populations: phased VCF in GT fomat -> Refined IBD -> IBDNe
# NAM: unphased GT -> IBDseq -> IbdNe # for homogeneous (as NAM or any other genetically uniform pop) 
# IBDSeq: https://faculty.washington.edu/browning/ibdseq.html#download
# IBDNe: https://faculty.washington.edu/browning/ibdne.html#download

# Nam is homogeneous population, the pipeline is a bit more simple
# Prepare file for IBDseq (we will use this for NAM)
# Get a nice population tags file for plink
awk 'FNR==NR{a[NR]=$2; next}{$6=a[FNR]}1'  Populations.tags.txt   BRn171-nam-eur-afr-amr.ordered.fam > Pops.to.plink.txt # simply putting populations' names in fam file

# This command will siplit large file into several different by population name (col 6, thus $6)
awk '{print >> $6; close($6)}' ./Pops/Pops.to.plink.txt # List of ind to all populations, quite nice

# Get NAM vcf file
plink --bfile BRn171-nam-eur-afr-amr.ordered --keep ./Pops/NAM --recode-vcf --out ./Ne/NAM.to.IBDSeq

### IBDseq is a software program for detecting segments of identity-by-descent (IBD) and homozygosity-by-descent (HBD) 
## in unphased genetic sequence data.
# IBDseq has two required arguments: the gt argument to specify an input VCF file, and the out argument to specify the output file prefix.
# The nthreads argument enables parallel computation. Here I am also using 7000Mb, that is, 7GB
 java -Xmx7000m -jar ibdseq.r1206.jar gt=./Ne/NAM.to.IBDSeq.vcf out=./Ne/NAM.to.IbdNe nthreads=8
 
# Get effective population size to NAM
# get plink map here: http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/
cat ./GRCh37/*.GRCh37.map > ./GRCh37/grch37.map # exclude X chromosome first 

# Run IBDNe. I think it might need to actually read chr by chr
# There were an error ERROR: Cannot bootstrap sample from a single chromosomes. Rerun the program with "nboots=0". So I did:
cat ./Ne/NAM.to.IbdNe.ibd | java -jar ibdne.19Sep19.268.jar map=./GRCh37/grch37.map out=./Ne/NAM.Ne nthreads=8 nboots=0

# This gives 'NAM.Ne.ne' file as output with Ne in different generations

#====== For admixed pipeline is different ====
# phased VCF in GT fomat -> Refined IBD -> IBDNe
# The difference is that we'll use phased input and the 'Refined IBD' instead of 'IBDseq'
# Refined IBD is a software program for detecting identity-by-descent (IBD) segments within and between individuals. 
# B L Browning and S R Browning (2013). Improving the accuracy and efficiency of identity by descent detection in population data. Genetics 194(2):459-71. doi:10.1534/genetics.113.150029

# Getting phased chr16 from different Latino pops
bcftools view -Ov -M2 -v snps -S ./Pops/BR.txt ALL.chr16.PHASED.vcf.gz > ./Ne/BR.to.RefinedIBD.vcf
bcftools view -Ov -M2 -v snps -S ./Pops/PUR.txt ALL.chr16.PHASED.vcf.gz > ./Ne/PUR.to.RefinedIBD.vcf
bcftools view -Ov -M2 -v snps -S ./Pops/CLM.txt ALL.chr16.PHASED.vcf.gz > ./Ne/CLM.to.RefinedIBD.vcf 
bcftools view -Ov -M2 -v snps -S ./Pops/MXL.txt ALL.chr16.PHASED.vcf.gz > ./Ne/MXL.to.RefinedIBD.vcf 
bcftools view -Ov -M2 -v snps -S ./Pops/PEL.txt ALL.chr16.PHASED.vcf.gz > ./Ne/PEL.to.RefinedIBD.vcf 

# Run Refined IBD
#There are only two required command line arguments: the gt argument to specify the input VCF file and the out argument to specify the output file prefix. A third parameter, the map parameter, is recommended, but not required. Other parameters have sensible default values
 java -Xss5m -Xmx8g -jar refined-ibd.17Jan20.102.jar gt=./Ne/BR.to.RefinedIBD.vcf map=./GRCh37/plink.chr16.GRCh37.map out=./Ne/BR.Ne nthreads=8 # uncompress the ouput (.ibd.gz)
 
# Run IBDNe. Same message:  ERROR: Cannot bootstrap sample from a single chromosomes. Rerun the program with "nboots=0". So I did:
cat ./Ne/BR.toIbdNe.ibd | java -jar ibdne.19Sep19.268.jar map=./GRCh37/plink.chr16.GRCh37.map out=./Ne/BR.Ne nthreads=8 nboots=0


###===================== 11. Calculate EHH alone & Explore iHS among populations further==================================
##== Theoretical stuff
#Extended haplotype homozygosity and relative EHH at a distance x from the core region is defined as the probability that two randomly chosen chromosomes carrying a tested core haplotype are homozygous at all SNPs8 for the entire interval from the core region to the distance x. EHH is on a scale of 0 (no homozygosity, all extended haplotypes are different) to 1 (complete homozygosity, all extended haplotypes are the same). Relative EHH is the ratio ofthe EHHon the tested core haplotype compared with the EHH ofthe grouped set ofcore haplotypes at the region not including the core haplotype tested. Relative EHH is therefore on a scale of 0 to infinity

#=== Sabeti, 2002
##===== But see Szpiech, Z. A. & Hernandez, R. D. selscan: An Efficient Multithreaded [...] (2014) for mathematical formulas.


# Although haplotype "selection" calculations (iHS and XP-EHH) were perfomed in topics 2 and 3 on this guide,
# right now I want to calculate EHH directly. 
# There are 2 options for this purpose: selscan --EHH (which seems to behave awkward, see Issues in git) and PopLDdecay, that have a nice output 

#PopLDdecay used in topic 9 does it better.
cd /Users/Pedro/Dropbox/input_chr16/

# Top values in chr16
# BRZ peak (iHS: -3.56) rs10871351
./PopLDdecay/bin/PopLDdecay -InVCF ALL.chr16.INFO.vcf.gz -OutStat ./LDfigures/Lddecay.BRZ.EHH.gz -SubPop ./Pops/BR.txt -EHH 16:78522906 

# NAM peak (iHS: 4.50) rs2966244
./PopLDdecay/bin/PopLDdecay -InVCF ALL.chr16.INFO.vcf.gz -OutStat ./LDfigures/Lddecay.NAM.EHH.gz -SubPop ./Pops/NAM.txt -EHH 16:82114741

# CLM peak (iHS: 4.22)
./PopLDdecay/bin/PopLDdecay -InVCF ALL.chr16.INFO.vcf.gz -OutStat ./LDfigures/Lddecay.CLM.EHH.gz -SubPop ./Pops/CLM.txt -EHH 16:78524099

# PEL peak (iHS: 3.84)
./PopLDdecay/bin/PopLDdecay -InVCF ALL.chr16.INFO.vcf.gz -OutStat ./LDfigures/Lddecay.PEL.EHH.gz -SubPop ./Pops/PEL.txt -EHH 16:11193930

# PUR peak (iHS: 3.70)
./PopLDdecay/bin/PopLDdecay -InVCF ALL.chr16.INFO.vcf.gz -OutStat ./LDfigures/Lddecay.PUR.EHH.gz -SubPop ./Pops/PUR.txt -EHH 16:78466758

# Alternatively, selscan calculates EHH alone with the flag --ehh, or in our CMS suite, like this
# scans.py selscan_ehh inputTped outFile locusID (--optionals). Avoid using this method, selscan --ehh is not working
source activate cms-env3

../chr16/cms/scans.py selscan_ehh chr16.scans.input.tped.gz /Users/Pedro/Dropbox/input_chr16/haplotype.scores/ locusID

# Get markers with suggestive selection ($8==1, get the modulo values of iHS and avg them). For chr16 only:
awk '$7 !~ /^-/{print $0}' ./haplotype.scores/chr16.BRZ.ihs.out.100bins.norm | awk '$8 == "1" {print $0}' | wc -l #140
awk '$7 !~ /^-/{print $0}' ./haplotype.scores/chr16.BRZ.ihs.out.100bins.norm | awk '$8 == "1"  {print $0}' | awk  '{sum+=$7} END { print "Average = ",sum/NR}' # Average =  2.41831

#NAM
awk '$7 !~ /^-/{print $0}' ./haplotype.scores/chr16.NAM.ihs.out.100bins.norm | awk '$8 == "1" {print $0}' | wc -l  #114
awk '$7 !~ /^-/{print $0}' ./haplotype.scores/chr16.NAM.ihs.out.100bins.norm | awk '$8 == "1"  {print $0}' | awk  '{sum+=$7} END { print "Average = ",sum/NR}' # Average =  2.41262

# CLM
awk '$7 !~ /^-/{print $0}' ./haplotype.scores/chr16.CLM.ihs.out.100bins.norm | awk '$8 == "1" {print $0}' | wc -l  #130
awk '$7 !~ /^-/{print $0}' ./haplotype.scores/chr16.CLM.ihs.out.100bins.norm | awk '$8 == "1"  {print $0}' | awk  '{sum+=$7} END { print "Average = ",sum/NR}' # Average =  2.49427

# PEL
awk '$7 !~ /^-/{print $0}' ./haplotype.scores/chr16.PEL.ihs.out.100bins.norm | awk '$8 == "1" {print $0}' | wc -l  #139
awk '$7 !~ /^-/{print $0}' ./haplotype.scores/chr16.PEL.ihs.out.100bins.norm | awk '$8 == "1"  {print $0}' | awk  '{sum+=$7} END { print "Average = ",sum/NR}' # Average =  2.37859

# PUR 
awk '$7 !~ /^-/{print $0}' ./haplotype.scores/chr16.PUR.ihs.out.100bins.norm | awk '$8 == "1" {print $0}' | wc -l  #125
awk '$7 !~ /^-/{print $0}' ./haplotype.scores/chr16.PUR.ihs.out.100bins.norm | awk '$8 == "1"  {print $0}' | awk   '{sum+=$7} END { print "Average = ",sum/NR}' # Average =  2.39703

##========================= Plot signals for different Latino populations ==============================
##== I have decided to extend analysis beyond chr16 using VariantAnnotation (R, to get synonymous or not), RCircos/SOFIA for circular plot and RGmatch (python: https://bitbucket.org/pfurio/rgmatch/src/master/) to classify variants as intronic, intergenic (upstream/downstream), promoter, exonic (1st_EXON and GENE_BODY), TSS and TTS (see https://bitbucket.org/pfurio/rgmatch/src/master). see also my script circularPlot.R
## I've installed Circos than RCircos pkg to plot, might also use SOFIA (Circos is needed)
# See script for SOFIA in https://academic.oup.com/jhered/article/108/4/443/3192401#supplementary-data
# Install modules as in http://www.circos.ca/documentation/tutorials/configuration/perl_and_modules/
# For Circos: basically  'brew install cpanminus' in my Mac or apt-get in a Linux and then take a look in missing modules
# Take a look at modules by doing 'cd bin' and './list.modules' also 'circos -modules'
#./circos -modules
# Install missing with cpan minus by
#sudo cpanm Clone Config::General Font::TTF::Font GD GD::Polyline Math::Bezier Math::Round Math::VecStat Params::Validate Readonly Regexp::Common SVG Set::IntSpan Statistics::Basic Text::Format
## Install missing
sudo cpanm -f GD   # This does the trick 
#sudo cpanm Module::Build::Tiny --force
#sudo cpanm Readonly
## Test it
#cd ../example/
#./run 
# Beggin by taking a look in quick start tutorial http://www.circos.ca/tutorials/lessons/quick_start/ or more broadly http://www.circos.ca/documentation/tutorials/

# In R do setwd("~/Dropbox/input_chr16/haplotype.scores/") then load("circos.RData") and draw the circos directly with (see CircularPlot.R)
SOFIA(data=iHS.gen,linkColor='chr',linksFlag=F,
      chromoConfiguration=chromoConfiguration,plotBackground=plotBackground,
      plotLocation=plotLocation,plotType=plotType,plotColor=plotColor,markerSize=markerSize,
      circosLocation='~/Dropbox/circos-0.69-9',tickSuffix='cM',returnConf=TRUE,circosDisplay=TRUE)

# circosDisplay=TRUE will create the circos.svg and png in R working directory, and regardless of this flag, conf files (3) will be created in ~/Dropbox/circos-0.69-9/bin/. Also, data will be written to /Dropbox/circos-0.69-9/data/. To check # of genes do 
wc -l  /Users/Pedro/Dropbox/circos-0.69-9/data/tr_5_1.txt #204 OK!

# In data directory, edit the ideo1.txt file (remove "1x" from chr names)
# To get backgrounds when there are no datapoints to be displayed, remove the "show=data" lines in circos.conf. Then, do
# (circos.png and circos.svg will be written in current directory - svg can be editted in GIMP)
~/Dropbox/circos-0.69-9/bin/circos -conf circos.conf

####======================== Overlapping of signals ===================================
##===== See Sabeti and Tourbenize. "selectionSabeti.csv" and "iHS.all.xlsx" spreadsheets have iHS signals summarized
# Gene lists are like /Users/Pedro/Dropbox/input_chr16/haplotype.scores/GenesUnderSelectionCLM.txt. I have also reunited selection signals and their overlapping in [use to overlap: http://bioinformatics.psb.ugent.be/webtools/Venn/]
# /Users/Pedro/Dropbox/input_chr16/haplotype.scores/selectionSabeti.csv
## See also Tournebize et al. McSwan: A joint site frequency spectrum method to detect and date selective sweeps across multiple population genomes. Mol. Ecol. Resour. 19, 283295 (2019).
# Copy from the supplementary material (men12957-sup-0004-tables3.xlsx) the genes under selection in CEU and LWK then
## (turning spaces into newlines)   | sed does not work for some reason with a similar regex (maybe gnu sed will, though)
tr ' ' '\n' < Tournebize > SelectionTournebize.txt # 848 genes ## See selectionSabeti.csv for gene-wise overview of literature

##========================================

##====================== RGmatch to get genes nearby/wherein markers lie ========================================
##===== There is still the problem of how to communicate the volume of iHS data and extract sense out of it
## Let's try the following rationale: selection events will, teoretically, more often affects genes than other genomic regions. Let's get sets of genes which iHS signals are nearby/in.
# Use RGmatch to this end https://bitbucket.org/pfurio/rgmatch/src/master/
# For input from populations, do in R when working on CircularPlot.R scipt: write.table(Pl, "~/Dropbox/input_chr16/haplotype.scores/iHS/PEL.xls",col.names = T, row.names = F, quote = F, sep="\t")

# ALL iHS - it can also output the gene SYMBOL use "--gene gene_name" (I discovered this later)
python ~/Downloads/pfurio-rgmatch-e9289746a5bd/rgmatch.py  \
    -g /Volumes/PEDRO\ 70GB/Genetics/GTF/hg19/gencode.v29lift37.annotation.gtf.gz \
    -b ~/Dropbox/input_chr16/haplotype.scores/iHS/CLM.xls \
    -o ~/Dropbox/input_chr16/haplotype.scores/iHS/CLM.genomic.regions.xls
        

# The above cmd will populate the output with many many overlaps for a same marker (for some reason TSS 240kb away from the query marker will show up as an associated marker), lets try to restrict associations to TSS
python ~/Downloads/pfurio-rgmatch-e9289746a5bd/rgmatch.py  \
    -r gene \
    -g /Volumes/PEDRO\ 70GB/Genetics/GTF/hg19/gencode.v29lift37.annotation.gtf.gz \
    -b ~/Dropbox/input_chr16/haplotype.scores/iHS/CLM.xls \
    -o ~/Dropbox/input_chr16/haplotype.scores/iHS/CLM.genomic.regions.gene.xls
 
 # Also for coding variants. Do in R: 
 #write.table(x =iHS.coding , file = "~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/codingMXL.xls", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
 python  ~/Downloads/pfurio-rgmatch-e9289746a5bd/rgmatch.py 
    -r gene
    -g /Volumes/PEDRO\ 70GB/Genetics/GTF/hg19/gencode.v29lift37.annotation.gtf.gz  
    -b ~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/codingMXL.xls 
    -o  ~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/MXL.rgmatch.txt 

# Everyone else
python  ~/Downloads/pfurio-rgmatch-e9289746a5bd/rgmatch.py -r gene  -g /Volumes/PEDRO\ 70GB/Genetics/GTF/hg19/gencode.v29lift37.annotation.gtf.gz  -b ~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/codingPEL.xls -o  ~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/PEL.rgmatch.txt  ;  ~/Downloads/pfurio-rgmatch-e9289746a5bd/rgmatch.py -r gene  -g /Volumes/PEDRO\ 70GB/Genetics/GTF/hg19/gencode.v29lift37.annotation.gtf.gz  -b ~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/codingCLM.xls -o  ~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/CLM.rgmatch.txt ;  ~/Downloads/pfurio-rgmatch-e9289746a5bd/rgmatch.py -r gene  -g /Volumes/PEDRO\ 70GB/Genetics/GTF/hg19/gencode.v29lift37.annotation.gtf.gz  -b ~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/codingPUR.xls -o  ~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/PUR.rgmatch.txt ;  ~/Downloads/pfurio-rgmatch-e9289746a5bd/rgmatch.py -r gene  -g /Volumes/PEDRO\ 70GB/Genetics/GTF/hg19/gencode.v29lift37.annotation.gtf.gz  -b ~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/codingBR.xls -o  ~/Dropbox/input_chr16/haplotype.scores/coding_VariantAnnotation_pkg/BR.rgmatch.txt 

####=============================== Estimate positive selection dates (tMRCA) ===============================
##See Voight 2006 and Kelley 2012
#========================================== Theoretical stuff ======================================================== 
############################################## Voight 2006 #################################################
#To obtain a crude estimate of the ages of sweeps, we assumed (1) a star phylogeny for the selected haplotypes, (2) that the two chromosomes which descended from a common ancestor and carry a selected allele are identical between the core SNP and the point ofthe nearest recombination event on either lineage, (3) that the two chromosomes start to differ immediately beyond the nearest recombination event. Hence, Pr[Homoz] = e(-2rg) where
# where Pr[Homoz] is the probability that two chromosomes are homozygous at a recombination distance r from the selected site, given a common ancestor g generations before the present. Taking the generation time to be 25y, the ancestor time in years becomes t = 25g (and, g  = t/25).
# The average total distance between the first point to the left, and to the right of the core SNP at which EHH on the selected haplotypes drops below 0.25 is 0.52 cM in both East Asians and Europeans, and 0.32 cM in Yoruba. Hence, for example in Africans, when Pr[Homoz] = 0.25, we observed that 2r = 0.32% (the distances in the text are the sum of r in both directions).


############################################## Kelley 2012 #################################################

#=======================================================================================================
# As PopLDdecay, selscan can compute EHH alone with (downside is no graph is generated)
# ./selscan --ehh <locusID> --vcf <vcf> --map <mapfile> --out <outfile>
# On the other hand, valuable genetic distance (cM) is generated, output details:
# Output: <physical dist> <genetic dist> <'1' EHH> <'0' EHH> #Default: __NO_LOCUS__
# With the scans.py wrapper from cms software:
#selscan_ehh [-h] [--gapScale GAPSCALE] [--maf MAF]
#                                       [--threads THREADS] [--window WINDOW]
#                                       [--cutoff CUTOFF]
#                                       [--maxExtend MAXEXTEND]
#                                       [--loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL,EXCEPTION}]
#                                       [--version] [--tmpDir TMPDIR]
#                                       [--tmpDirKeep]
##                                       inputTped outFile locusID
#conda activate cms-env2 <<<<- completly wrong, use PopLDdecay instead
#~/Dropbox/chr16/cms/scans.py selscan_ehh ~/Dropbox/input_chr16/TPED/chr16.CLM.scans.input.tped.gz ~/Dropbox/input_chr16/haplotype.scores/CLM.test 57503213 --cutoff 0.001 --maxExtend 10000000

#====== Run in-house script tMRCA.sh
# See the script tMRCA.sh, outputs are in the format BR.GeneticDistances_tMRCA.txt

# For 1 percentil most extreme iHS values in each Latino pop do
cat /Users/Pedro/Dropbox/input_chr16/haplotype.scores/iHS/extreme.iHS.Mx.txt | sed 's/  */,/g' | sed 's/_/,/g' > ~/Dropbox/input_chr16/haplotype.scores/iHS/extreme.iHS.Mx.csv
# Run tMRCA.sh
 ~/Documents/GitHub/tMRCA/tMRCA.sh -p MXL -f ~/Dropbox/input_chr16/haplotype.scores/iHS/extreme.iHS.Mx.csv -o ~/Dropbox/input_chr16/haplotype.scores/tMRCA/1percentil/

# Get averages
cat /Users/Pedro/Dropbox/input_chr16/haplotype.scores/BR.GeneticDistances_tMRCA.txt | grep tMRCA | awk '{if( $2 >1000) {print $2}}' |awk '{ total += $1; count++ } END { print total/count }' 
# 191047
cat /Users/Pedro/Dropbox/input_chr16/haplotype.scores/PUR.GeneticDistances_tMRCA.txt | grep tMRCA | awk '{if( $2 >1000) {print $2}}' |awk '{ total += $1; count++ } END { print total/count }'
# 128652
cat /Users/Pedro/Dropbox/input_chr16/haplotype.scores/CLM.GeneticDistances_tMRCA.txt | grep tMRCA | awk '{if( $2 >1000) {print $2}}' |awk '{ total += $1; count++ } END { print total/count }'
# 108086
cat /Users/Pedro/Dropbox/input_chr16/haplotype.scores/MXL.GeneticDistances_tMRCA.txt | grep tMRCA | awk '{if( $2 >1000) {print $2}}' |awk '{ total += $1; count++ } END { print total/count }'
# 93668
cat /Users/Pedro/Dropbox/input_chr16/haplotype.scores/PEL.GeneticDistances_tMRCA.txt | grep tMRCA | awk '{if( $2 >1000) {print $2}}' |awk '{ total += $1; count++ } END { print total/count }'
# 64221.8

#=================== Demography simulation with SLiM  =============##=======================
##=================== ##=================== ##=================== ##=========================
##=================== ##=================== ##=================== ##====##=================== 
##=========================================================================================== 

# SLiM has a nice documentation, take a look in manual (Mendeley) and workshops: http://benhaller.com/workshops/workshops.html
# Also in the discuss group: https://groups.google.com/g/slim-discuss
# For Latinos, here are nice scripts for simulating (isolate populations) population dynamics: https://github.com/LohmuellerLab/ROH_Latin_American_Isolates - from paper from Mooney et al. 2018
# My own script is in ~/Dropbox/paperultima/ 

# in SLiMGui, load previous populations in mutation-drift equilibrium by doing (script in ~/Dropbox/paperultima/LoadPopulation.slim)
1 { sim.readFromPopulationFile("~/Desktop/SLiM_models/slim_1676180791791.txt");
           }

# Pass the following command to generate a vcf file for downstream analyses. (It has to be copied from output section, modify header for vcf=4.3 instead of 4.2)
71667 late() { p1.outputVCFSample(80); }

bgzip -c  Latinos_out.vcf > Latinos_out.vcf.gz 
tabix -p vcf Latinos_out.vcf.gz 

# To quickly generate LD decay:
~/Dropbox/input_chr16/PopLDdecay/bin/PopLDdecay -InVCF Latinos_out.vcf.gz  -OutStat Lddecay.stat.gz

# Now plot
perl ../../input_chr16/PopLDdecay/bin/Plot_OnePop.pl -inFile LDDecay_B125_i50.stat.gz -output Latino.simulation.B125i50

# Or get r2 with plink
plink --vcf Latinos_outBottleck1000_inbreed75.vcf --maf 0.05 --r2

# For the empirical populations
 plink --vcf /Users/Pedro/Dropbox/input_chr16/ALL.chr16.INFO.vcf --double-id  --maf 0.05 --keep /Users/Pedro/Dropbox/input_chr16/Pops/BR --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 90200000 --out BR
# get average
awk  '{sum+=$7} END { print "Average = ",sum/NR}' ./Testing\ LD/plink.ld


# Uncompress Lddecay.stat.gz and open in excel, use the following to calculate avg r^2 (>2kb apart)
=AVERAGEIF(A:A, ">2000",B:B)
# This gets the following values
#### Empirical (chr16) plink   #### SLiM                                              LDDecay       plink - MAF > 0.05**
## BR:  0.034614054   0.525587  # Bottleneck Ne=1000 only:                           0.025204705   Average =  0.658736
## PUR: 0.043170909   0.546617  # Bottleneck Ne=100:                                 0.066492911   Average =  0.672229
## CLM: 0.04762064    0.550755  # Bottleneck Ne=100 and inbreeding (50% of the time):0.064234373   Average =  0.675566
## MXL: 0.051648419   0.56722   # Bottleneck Ne=50 and inbreeding (75% of the time): 0.111001131   Average =  0.677671
## PEL: 0.055237637   0.58506   # Bottleneck Ne=1000 and inbreeding (40%):           0.026735256   Average =  0.654722
                                # Bottleneck Ne=100 and inbreeding (60%)             0.063754087   Average =  0.667909 
                                # Bottleneck Ne=75 and inbreeding (50%)              0.077405142   Average =  0.67247 
                                # Bottleneck Ne=1000 and 75% inbreed                               Average =  0.656759
        #This is ~BR!! >>>>>>>  # Bottleneck Ne=500 and 30% inbreed                  0.032724057   Average =  0.647615             
                                # Bottleneck Ne=300 and 40% inbreed                  0.039164192   Average =  0.654451 
      #This is ~MXL!! >>>>>>>   # Bottleneck Ne=175 and 50 inbreed                   0.051076464   Average =  0.659617
                                # Bottleneck Ne=100 and inbreeding (60%)             0.063754087   Average =  0.667909 
      #This is ~MXL!! >>>>>>>   # Bottleneck Ne=150 and 50 inbreed                   0.052340127   
                                # Bottleneck Ne=125 and 70 inbreed                   0.058200559   
      #This is ~PEL!! >>>>>>>   # Bottleneck Ne=125 and 50 inbreed                   0.0576711   


#     ** plink --r2 command automatically removes r2 < 0.2 (this is why averages are always higher)

# Now for LD plotting to check if simulation ~ to empirical
cd /Users/Pedro/Dropbox/input_chr16/LDfigures 
perl ../PopLDdecay/bin/Plot_MultiPop.pl -inList Pop.ReslutPath.list  -output ALL.sim #just include simulation LDDecay result in Pop.ReslutPath.list 

# See biostars.org/p/300381/
# See https://www.biostars.org/p/347796/ also 


##=================== LD decay analysis and plotting  =============##=======================
##=================== ##=================== ##=================== ##=========================
##=================== ##=================== ##=================== ##====##=================== 
##=========================================================================================== 

############################################################################################################################
####### For empirical populations (BRn171-nam-eur-afr-amr.ordered has 302,369 markers because merges 5.0 and 6.0 and no QC)
# Based on https://www.biostars.org/p/347796/ also 
#################################################################################################################### 
plink --bfile ../input_chr16/BRn171-nam-eur-afr-amr.ordered --maf 0.01 --mind 0.1 --geno 0.1 --hwe 0.000000001 --keep ../input_chr16/Pops/PEL --recode-vcf --out ./Populations/PEL.all.markers
plink --bfile ../input_chr16/BRn171-nam-eur-afr-amr.ordered --maf 0.01 --mind 0.1 --geno 0.1 --hwe 0.000000001 --keep ../input_chr16/Pops/MXL --recode-vcf --out ./Populations/MXL.all.markers 
plink --bfile ../input_chr16/BRn171-nam-eur-afr-amr.ordered --maf 0.01 --mind 0.1 --geno 0.1 --hwe 0.000000001 --keep ../input_chr16/Pops/CLM --recode-vcf --out ./Populations/CLM.all.markers 
plink --bfile ../input_chr16/BRn171-nam-eur-afr-amr.ordered --maf 0.01 --mind 0.1 --geno 0.1 --hwe 0.000000001 --keep ../input_chr16/Pops/PUR --recode-vcf --out ./Populations/PUR.all.markers 
plink --bfile ../input_chr16/BRn171-nam-eur-afr-amr.ordered --maf 0.01 --mind 0.1 --geno 0.1 --hwe 0.000000001 --keep ../input_chr16/Pops/BR --recode-vcf --out ./Populations/BR.all.markers

# After filtering, kept 41598 out of a possible 302369 Sites (for all populations)
for i in $(ls ~/Dropbox/paperultima/Populations/ | cut -d . -f 1);    
    do vcftools --vcf ./Populations/$i.all.markers.vcf --recode --recode-INFO-all --thin 50000 --out ./LD\ decay/$i.thin;   
done

for i in $(ls ~/Dropbox/paperultima/LD\ decay/| cut -d . -f 1);    
    do mv ./LD\ decay/$i.thin.recode.vcf ./LD\ decay/$i.thin.vcf;   
done

for i in $(ls ~/Dropbox/paperultima/LD\ decay/| cut -d v -f 1);        
    do bcftools view ./LD\ decay/$i"vcf" -O b -o ./LD\ decay/$i"bcf"; 
done


# Instructions on using tomahawk is detailed here. Convert BCF to the tomahawk format. 
## note: BR fails to import to tomahawk, idk why. I comented /lib/core.ccp lines 116 and 117 and complie a different instance of tomahawk just to import BR
for i in $(ls *bcf | cut -d b -f 1);    # cd to ./LD decay/
    do ../tomahawk/tomahawk import -i $i"bcf" -o $i ; 
done

# Arguments: -u (unphased), -i (input file), -o (output file name prefix), -r (min R2 cut-off). See tomahawk help for more descriptions.
 for i in $(ls *twk | cut -d . -f 1,2);       # this takes quite some time
    do ../tomahawk/tomahawk calc -ui $i.twk -o $i  -r 0.1 -P 0.1 -C 1 -t 1;  
done

# The .two output is binary and needs to be exported to a text file. The -I arguments restricts the output to a specified chromosome. So, -I chr1 (UCSC) of -I 1 (Ensembl) exports only chromosome 1. It is probably better to have separate files for each chromosome as the files are a more managable size to read into R. The LD metrics are computed between SNPs within chromosomes anyway.
../tomahawk/tomahawk view -i PEL.thin.two -I 1 > pel.chr1.ld

## Use script LD_decay.R to plot (Roteiros/)

####################################################################################################################
#################### For SLiM simulation-generated population ###########
####################################################################################################################
# So, iHS is calculated for each population, my idea is to a qq-plot of p-values on simulated and each Latino population. 
# Note that each population has between 64 (MXL) and 171 (BR) individuals, so the calculation will take this 171 individuals and
# multiply to the # of autosomal chr (22) and the average size in MB (total sequenced size: 2684.5, divided by 22: 122Mb):
# 171 * 22 * 122 = 458964 ind*chr*Mb 
# Since I am simulating 10Mb chr, I'll need 45897 individuals from SLiM (each ind is a diploid for a single chr). I will simplfy
# by generating 5 samples of 10,000, getting 50,000 diploid individuals with 10Mb chromosomes (then I will change the chr #
# in each sample to merge all, thus creating 10,000 samples of individuals diploid for 5 chromosomes) 

## code sample for qq-plot: https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
# This post is useful on the interpretaion of a qqplot

# To test LD decay (use LDDecay with multi_pop option, look above)
# So I tried a really aggressive LD approach in my simulated populations (bottleneck of 100, 70% inbreeding and only recovering 
# from bottleneck 3 generations ago)
plink --vcf  ~/Dropbox/paperultima/Simulation/Latinos_Bottleck100_inbreed70.vcf.gz --maf 0.01 --mind 0.1 --geno 0.1 --hwe 0.000000001  --recode-vcf --out Bottleck100_inbreed70 # QC

bcftools view Bottleck100_inbreed70.vcf -O b -o Bottleck100_inbreed70.bcf 

# change chromosome numbers
echo "1 2" >> chr_name_conv.txt
echo "1 3" >> chr_name_conv.txt
echo "1 4" >> chr_name_conv.txt
echo "1 5" >> chr_name_conv.txt
bcftools annotate --rename-chrs chr_name_conv.txt original.vcf.gz  | bgzip > rename.vcf.gz

# concat since they are the "same" individuals (merge for different individuals), open, change to VCFv=4.3 and bgzip/tabix
bcftools concat --file-list bcfs.to.merge.txt | bcftools norm -d all -o Latinos_Bottleck100_inbreed70.all.vcf
# 108185 variants and 10000 people pass filters and QC.
plink --vcf /Users/Pedro/Dropbox/paperultima/Simulation/Latinos_Bottleck100_inbreed70.vcf.gz  --make-just-bim --out checking 

# Here I am getting the 1st position of each chr (variable $first), importing it to awk (as 'var') than subtracting it from every 
# other physical position (i.e. distance), multiplying the result by 1e-8 (recomb rate I set in SLiM by generation by bp) 
# and adjusting to be more realistic. See https://www.biostars.org/p/178465/#9503825

awk '{print>$1}' checking.bim # split the .bim by chr is easier to calculate genetic distances
for file in {1..5};
    do first=$(awk 'NR==1{print $4}' $file); 
#    {print $1 "\t" $2 "\t" $4-$first "\t" $4}' $file;
    awk -v var="$first" '{print $1 "\t" $2 "\t" ($4-var)*0.00000001*1000 "\t" $4}' $file > $file.map;
done

# Make sure physical positions are sorted if needed (pipe to this: | sort -k1,1 -k4,4n)

# selscan works with one chr at a time :/
 bcftools filter Latinos_Bottleck100_inbreed70.vcf.gz -r 1 > Latinos_Bottleck100_inbreed70.chr1.vcf.gz 
  bcftools filter Latinos_Bottleck100_inbreed70.vcf.gz -r 2 > Latinos_Bottleck100_inbreed70.chr2.vcf.gz 
   bcftools filter Latinos_Bottleck100_inbreed70.vcf.gz -r 3 > Latinos_Bottleck100_inbreed70.chr3.vcf.gz 
    bcftools filter Latinos_Bottleck100_inbreed70.vcf.gz -r 4 > Latinos_Bottleck100_inbreed70.chr4.vcf.gz 
     bcftools filter Latinos_Bottleck100_inbreed70.vcf.gz -r 5 > Latinos_Bottleck100_inbreed70.chr5.vcf.gz 

# Finally, calculate iHS in simulated dataset. This takes a lot of time - Ive aboned this approach
#~/Documents/Sickle\ Paper/Original\ Files/updated/Haplotypes/selscan/bin/osx/selscan --ihs  \
#                        --vcf Latinos_Bottleck100_inbreed70.chr1.vcf.gz \
#                        --map 1.map \
#                        --maf 0.05 

# It is a lot of individuals and seems to me that: 1 - selscan can't handle and 2 - will hamper false positives
# Possible solution: break in 80 ind lists and run it one at a time, so sampling might get false iHS-positives 
# "chr1" has 21562 markers (grep -v "^#" ~/Dropbox/paperultima/Simulation/Latinos_Bottleck100_inbreed70_chr1.vcf | wc -l)
query -l Latinos_Bottleck100_inbreed70.chr1.vcf.gz > samples/all.samples.txt # get list of samples
cd ./samples/
for f in *.txt; do gsplit -d -a3 -l80 --additional-suffix=.txt "$f" "${f%.txt}-"; done # break into 125 lists
rm all.samples.txt

# Lets automate the rest: generate vcf files from the 125 lists of ind
for sample in *txt; do     
    bcftools view -Oz -S $sample  -o $sample.vcf.gz ../Latinos_Bottleck100_inbreed70.chr1.vcf.gz; 
    
done

# iHS for all 80-individuals list (125 in total). It took ~72h to complete
for file in *.vcf.gz; do 
     ~/Documents/Sickle\ Paper/Original\ Files/updated/Haplotypes/selscan/bin/osx/selscan --ihs  \
                        --vcf $file \
                        --map ../1.map  \
                        --ihs-detail \
                        --threads 8 \
                        --maf 0.05 \
                        --out $file ;             
done

# normalize 'em - incredibly quick
for file in *.ihs.out; do 
     ~/Documents/Sickle\ Paper/Original\ Files/updated/Haplotypes/selscan/bin/osx/norm --ihs \
     --files $file;                      
done

for file in *.norm; do     # yields 32721 signals (Ive got 29175 unique in all empirical samples) ranging from -5.46596 to 4.68154
    awk '$8=="1"{print $0}' $file >> all.simulated.signals.txt;
done    

for file in ~/Dropbox/input_chr16/haplotype.scores/*.norm; do  # yields 50832, some unique, some shared between populations
    awk '$8=="1"{print $0}' $file >> all.empirical.signals.txt;      #ranging from -4.51503 to 5.35662
done    

# get iHS for empirical
for chr in {1..22} ; do awk '$8=="1"{print $0}' chr$chr.MXL.ihs.out.ihs.out.100bins.norm >> MXL.signals.txt ; done && wc -l MXL.signals.txt #9554 signals (97 surpassed simulations [look below]: 1.015 percentile)

for chr in {1..22} ; do awk '$8=="1"{print $0}' chr$chr.PEL.ihs.out.ihs.out.100bins.norm >> PEL.signals.txt ; done && wc -l PEL.signals.txt # 8606 signals (surpassed: 134: 1.56 percentile)

for chr in {1..22} ; do awk '$8=="1"{print $0}' chr$chr.PUR.ihs.out.ihs.out.100bins.norm >> PUR.signals.txt ; done && wc -l PUR.signals.txt # 10186 signals (surpassed: 12, 0.12 percentile)

for chr in {1..22} ; do awk '$8=="1"{print $0}' chr$chr.CLM.ihs.out.ihs.out.100bins.norm >> CLM.signals.txt ; done && wc -l CLM.signals.txt # 9804 signals (surpassed: 28, 0.28559 percentile)

for chr in {1..22} ; do awk '$8=="1"{print $0}' chr$chr.BR.ihs.out.ihs.out.100bins.norm >> BR.signals.txt ; done && wc -l BR.signals.txt  #11066 signals (surpassed: 2, 0.01807 percentile)


### I now reckon I have to match number of tests ran in each, simulated and empirical populations, so I need total markers 
# in each iHS evaluater-population and than match to iHS from the same number of markers in simulation
for file in $(ls *PEL.ihs.out.ihs.log); do grep Removed $file | cut -d" "  -f2  ; done | awk '{ total += $1 } END { print total }' 
# 75314 markers removed from MXL, so there are 167390-75314 = 92076 tests in PEL 
for file in $(ls *MXL.ihs.out.ihs.log); do grep Removed $file | cut -d" "  -f2  ; done | awk '{ total += $1 } END { print total }' 
# 55829 markers removed from MXL, so there are 167390-55829 = 111561 tests in MXL
for file in $(ls *CLM.ihs.out.ihs.log); do grep Removed $file | cut -d" "  -f2  ; done | awk '{ total += $1 } END { print total }'
# 48610 removed, thus 118780 tests
for file in $(ls *PUR.ihs.out.ihs.log); do grep Removed $file | cut -d" "  -f2  ; done | awk '{ total += $1 } END { print total }' 
# 43449 removed, thus 123941 tests
for file in $(ls *BR.ihs.out.ihs.log); do grep Removed $file | cut -d" "  -f2  ; done | awk '{ total += $1 } END { print total }' 
# 34071 removed, thus 133319 tests

# Now lets see how many samples we need to match the above empirical data
# wc -l /Users/Pedro/Dropbox/paperultima/Simulation/1.map gives 21529 markers in total, so do the following to see how many 
# removed markers each sample got than work with numbers in Ratio.tests.simulation.emprirical.xlsx
for file in $(ls *.ihs.log); do grep Removed $file | cut -d" "  -f2  ; done  

# I figured what simulations I need to ~ match tests in each population (I actually selected row in R, see qqplot.R)
for n in {0..7} ; do awk '$8=="1"{print $0}' all.samples-00$n.txt.vcf.gz.ihs.out.100bins.norm >> simulated.iHS.PEL.txt  ; done # PEL 2235 signals
for n in {8..9} ; do awk '$8=="1"{print $0}' all.samples-00$n.txt.vcf.gz.ihs.out.100bins.norm >> simulated.iHS.MXL.txt  ; done
for n in {10..16} ; do awk '$8=="1"{print $0}' all.samples-0$n.txt.vcf.gz.ihs.out.100bins.norm >> simulated.iHS.MXL.txt  ; done #  MXL 2269
for n in {17..26} ; do awk '$8=="1"{print $0}' all.samples-0$n.txt.vcf.gz.ihs.out.100bins.norm >> simulated.iHS.CLM.txt  ; done # CLM 2492
for n in {38..48} ; do awk '$8=="1"{print $0}' all.samples-0$n.txt.vcf.gz.ihs.out.100bins.norm >> simulated.iHS.PUR.txt  ; done # PUR 3018
for n in {27..37} ; do awk '$8=="1"{print $0}' all.samples-0$n.txt.vcf.gz.ihs.out.100bins.norm >> simulated.iHS.BR.txt  ; done # BR 3137

# Now let's trim the resulting iHS values fot simulation-derived cutoffs
# make sure chr and pos are sorted (I did it manually) BUT, here's a nice tip, use gnu sort with V option
# gsort -k1,1V file (alpha-numeric sort, see https://www.biostars.org/p/64687/)
sort -rk7,7 /Users/Pedro/Dropbox/input_chr16/haplotype.scores/iHS/BR.xls  | awk '$8 > 4.21354 {print $0"\tBR"  }' | gsort -k1,1V -k2,2n > iHS.surpassing.simulation.xls
sort -rk7,7 /Users/Pedro/Dropbox/input_chr16/haplotype.scores/iHS/PUR.xls  | awk '$8 > 4.12386 {print  $0"\tPUR"  }' | gsort -k1,1V -k2,2n  | awk 'NR>1{print $0}'  >> iHS.surpassing.simulation.xls
sort -rk7,7 /Users/Pedro/Dropbox/input_chr16/haplotype.scores/iHS/CLM.xls  | awk '$8 > 3.82635 {print  $0"\tCLM"  }' | gsort -k1,1V -k2,2n  | awk 'NR>1{print $0}'  >> iHS.surpassing.simulation.xls
sort -rk7,7 /Users/Pedro/Dropbox/input_chr16/haplotype.scores/iHS/MXL.xls  | awk '$8 > 3.38214 {print  $0"\tMXL"  }' | gsort -k1,1V -k2,2n  | awk 'NR>1{print $0}'  >> iHS.surpassing.simulation.xls
sort -rk7,7 /Users/Pedro/Dropbox/input_chr16/haplotype.scores/iHS/PEL.xls  | awk '$8 > 3.29419 {print  $0"\tPEL"  }' | gsort -k1,1V -k2,2n  | awk 'NR>1{print $0}'  >> iHS.surpassing.simulation.xls


# Get the genic ones (173 are strictly intronic/exonic)
# -q 0: map to the exact position the SNP lays 
# --gene gene_name gets gene SYMBOL right on (no need to biomart)
python2 ~/Downloads/pfurio-rgmatch-e9289746a5bd/rgmatch.py -g /Volumes/PEDRO\ 70GB/Genetics/GTF/hg19/gencode.v29lift37.annotation.gtf.gz \
        --gene gene_name \
        -q 0 \
        -b iHS.surpassing.simulation.xls \
        -o iHS.surpassed.genomic.regions.0kb.xls

# Just conforming the genes we got (way easier than using UCSC genome browser to inspect) - one can get symbols directly from RGmatch. This is also good for a short functional view on each gene
# This gets us a lot of ensembl-notation gene IDs (ENSG000001), getting HGNC symbols
## UPDATE: I think SNPEff will be better suited for this purpose, see below
awk ' { print $3 } ' iHS.surpassed.genomic.regions.xls | cut -d . -f1 | uniq # copy and go to https://www.ensembl.org/biomart/martview/42d4cef6d1b039820dadf5c7d0ba9163
# select Ensembl Genes 105 and GRCh38, then in Filters go to Gene and then "Input external references ID" as paste IDs
# Attributes panel is what we expect to retrieve from biomart, go External > HGNC symbol
# Attributes > Phenotype description is also interesting

###### Now get tMRCA on those 
#remember tab is 'control + v' than hit tab in Mac OS
# ADD HEADER TO THEM, LDDecay (and thus, tMRCA script) EXPECTS HEADER (see BR output for a header)
sed 's/      /,/g'  iHS.surpassing.simulation.xls | grep -w  PEL > PEL.to.tMRCA.txt
sed 's/      /,/g'  iHS.surpassing.simulation.xls | grep -w  MXL > MXL.to.tMRCA.txt
sed 's/      /,/g'  iHS.surpassing.simulation.xls | grep -w  CLM > CLM.to.tMRCA.txt
sed 's/      /,/g'  iHS.surpassing.simulation.xls | grep -w  PUR > PUR.to.tMRCA.txt
sed 's/      /,/g'  iHS.surpassing.simulation.xls | grep -w  BR > BR.to.tMRCA.txt


# Run tMRCA.sh
# get avg like this: cat ./tMRCA/PEL.GeneticDistances_tMRCA.txt | grep tMRCA | awk '{if( $2 >1000) {print $2}}' |awk '{ total += $1; count++ } END { print total/count }' 
~/Documents/GitHub/tMRCA/tMRCA.sh -p PUR -f PUR.to.tMRCA.txt -o ./tMRCA/ # avg 17608.2
~/Documents/GitHub/tMRCA/tMRCA.sh -p BR -f BR.to.tMRCA.txt -o ./tMRCA/   # 15484 
~/Documents/GitHub/tMRCA/tMRCA.sh -p CLM -f CLM.to.tMRCA.txt -o ./tMRCA/ # avg 77102.5
~/Documents/GitHub/tMRCA/tMRCA.sh -p MXL -f MXL.to.tMRCA.txt -o ./tMRCA/ # avg 44322.8
~/Documents/GitHub/tMRCA/tMRCA.sh -p PEL -f PEL.to.tMRCA.txt -o ./tMRCA/ # avg  18310.4

# Since tMRCA output comes sorted by coordinates (within each population), one can simply grep "tMRCA" and paste into iHS.surpassing.simulation.xls. Like so: cat ./tMRCA/PEL.GeneticDistances_tMRCA.txt | grep tMRCA  

# Now get unique genic markers in each population (remember, here are only GENIC variants, either intron or exon) 
# 173 total markers (161 intronic [86 1st intron] and 12 exons) spanning 120 unique genes
# 39 overlapping markers in all populations (exact same position)
# thus, 173-(120+39) = 14 genes overlapping to different positions in same and/or different populations
sort -k17,17 iHS.surpassed.genomic.regions.0kb.xls | awk '!seen[$1,$2,$17]++' > iHS.surpassed.genomic.regions.genic.uniq.xls # this is input to circular plot

# With little modification, this will become a supplementary table (maybe complement with phenotypical info from biomart)
gsort -k1,1V -k2,2n iHS.surpassed.genomic.regions.genic.uniq.xls | sed 's/_[0-9]*//g'  > Supplementary.tables.xls 

## We still need the remaining 100 markers (intergenic ones), let's get first a comprehensive list for all markers
python2 ~/Downloads/pfurio-rgmatch-e9289746a5bd/rgmatch.py -g /Volumes/PEDRO\ 70GB/Genetics/GTF/hg19/gencode.v29lift37.annotation.gtf.gz         --gene gene_name         -q 700         -b iHS.surpassing.simulation.xls         -o iHS.surpassed.genomic.regions.700kb.xls

# Then exclude the genic markers to get intergenic alone - just adjust a little to get a supplementary table
awk 'NR==FNR{a[$1];next} !($1 in a)'  iHS.surpassed.genomic.regions.genic.uniq.xls iHS.surpassed.genomic.regions.700kb.xls > iHS.surpassed.genomic.regions.intergenic.xls 

# Finally export to add to supplementary tables
gsort -k1,1V -k2,2n iHS.surpassed.genomic.regions.intergenic.xls | sed 's/_[0-9]*//g'  > Supplementary.tables2.xls

############## ideas to futher analyze selection signals #################################################################
# - get better tMRCA estimates (~starTMRCA, see Joel Smith's 2008 paper)
# - and/or explore human genome dating database: https://human.genome.dating/snp/rs182549
# - inspect variants function impact, SNPEff to the rescue: http://pcingola.github.io/SnpEff/examples/ 
# - (or the https://pophumanvar.uab.cat/ pipeline, already integrated to SNPEff)
# - Different tMRCA within populations for the same variant: populations size and selection coefficient are responsible?
####################=========================================================#############

## Some quick SNPEff evaluation 
cat ../iHS.surpassing.simulation.xls | gsort -k1,1V -k2,2n | awk 'NR>1{print $1"\t"$2}' | sed 's/chr//g'   > markers.to.effect.txt 

bcftools filter /Volumes/PEDRO\ 70GB/Genetics/ALLinfo/VCFinfo/ALL.vcf.gz -R markers.to.effect.txt | bgzip > markers.to.effect.vcf.gz

# get general stats output
java -Xmx8g -jar ~/snpEff/snpEff.jar -v -stats effect.html  GRCh37.75 markers.to.effect.vcf > markers.effect.ann.vcf

# get rsID for supplementary table
awk 'NR==FNR {vcf[$1 FS $2] = $3} {ihs = $1 FS $2} ihs in vcf {print ihs "\t"  vcf[ihs]} ' markers.to.effect.vcf  ../iHS.surpassing.simulation.xls

#####======================== Get marker location information (intron, exon, etc) ================================
# get a BED-like file for rgmatch
#####=============================================================================================================
awk '{ print "chr"$1"\t"$4"\t"$4}' /Users/Pedro/Dropbox/input_chr16/GALINAthinned.map > array.markers.to.gdmatch.xls

# important: output is unsorted (chr order is messy), but all chr are included
# -q 700 up to 700kb associations, in order to comprehend all markers
# -r gene outputs only genic SNPs (itron, 1st_EXON and GENE_BODY: exons other than 1st one)
python2 ~/Downloads/pfurio-rgmatch/rgmatch.py -g /Volumes/PEDRO\ 70GB/Genetics/GTF/hg19/gencode.v29lift37.annotation.gtf.gz \   
        -q 700 \
        -b array.markers.to.gdmatch.xls \
        -o array.markers.info.exonic700kb.txt   


# Check if all markers were included (ie, 176390 plus header). uniq because RGmatch reports 
# intron/up/downstream assignments for different transcripts, thus creating ambiguosity
awk '{print $1}' array.markers.info.exonic.txt | sort | uniq | wc -l


# 4th column in RGmatch output is transcript, and a SNP may have different distance to TSS, intron position or intron/exon
# assignment depending on the transcript, yielding duplicate entries for some of the markers. Lets keep only first occurrence
# based on the coordinates column (1st). awk creates an array called 'seen' and increments with values in 1st column, the !
# sign will remove the next occurrences of these values, thus keeping only 1st lines for each value in column 1.   
awk '!seen[$1]++' array.markers.info.exonic700kb.txt > array.markers.info.all.txt

# Finally, get numbers
grep INTRON   array.markers.info.all.txt | wc -l #  7264
grep STREAM   array.markers.info.all.txt | wc -l # 159649
grep 1st_EXON   array.markers.info.all.txt | wc -l #  2578
grep GENE_BODY   array.markers.info.all.txt | wc -l #298
grep TSS  array.markers.info.all.txt | wc -l  #1033-1 = 1032 (header)
grep PROMOTER  array.markers.info.all.txt | wc -l #5569

# Get distances between markers and standard deviation
# for chr 1 ($1=="1" - I did manually by chrs)
awk '{if( $1 == "1") {print $4}}' ../ROH/BR_n171.thin6.bim |  #get markers position
    awk 'NR > 1 { print $0 - prev } { prev = $0 }' |          # subtract each position by previous and proceed to mean/std
    awk '/^#/{next}{for(i=1;i<=NF;i++) {sum[i] += $1; sumsq[i] += ($1)^2}}    
          END {for (i=1;i<=NF;i++) {
          printf "%f %f \n", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR)}
         }'; done


####=== Now for estimates os nucleotide diversity, extract each population of interest and use TASSEL Analysis>Diversity panel
## Use "Open as" and than set plink file
for pop in $(echo ACB ASW CEU CLM ESN FIN GBR GWD IBS LWK MSL MXL NAM PEL PUR TSI YRI) ; do 
    plink --bfile ../BRn171-nam-eur-afr-amr.ordered \
    --extract ../ROH/BR_n171.thin6.bim \
    --keep ~/Dropbox/input_chr16/Pops/$pop \
    --recode --tab --out $pop; done


###===== IMPORTANT: nucleotide divresity is inflated in array SNP panel, so let's proceed to haplotype diversity in vcflib
bcftools query -l /Volumes/PEDRO\ 70GB/Genetics/ALLinfo/VCFinfo/ALL.vcf.gz > indlist.txt # for order they occour in vcf file

for pop in $(echo ACB ASW BR CEU CLM ESN FIN GBR GWD IBS LWK MSL MXL NAM PEL PUR TSI YRI) ; do cat ~/Dropbox/input_chr16/Pops/$pop | sed 's/ /      /g'  >> poplist.txt; done  # for populations membership (indlist.xls for combined population and position in vcf)


# use seq to get the list of ind to vcflib
# BR
for i in $(seq -s, 894 1063); do /Users/Pedro/Documents/Sickle\ Paper/Original\ Files/updated/Haplotypes/vcflib.source/sequenceDiversity --file  /Volumes/PEDRO\ 70GB/Genetics/ALLinfo/VCFinfo/ALL.vcf.gz --type GT  --target $i  ; done > HaplotypeDiversity/BR.hap.txt

# PUR
for i in $(seq -s, 789 892); do /Users/Pedro/Documents/Sickle\ Paper/Original\ Files/updated/Haplotypes/vcflib.source/sequenceDiversity --file  /Volumes/PEDRO\ 70GB/Genetics/ALLinfo/VCFinfo/ALL.vcf.gz --type GT  --target $i  ; done > HaplotypeDiversity/PUR.hap.txt

# CLM
for i in $(seq -s, 695 788); do /Users/Pedro/Documents/Sickle\ Paper/Original\ Files/updated/Haplotypes/vcflib.source/sequenceDiversity --file  /Volumes/PEDRO\ 70GB/Genetics/ALLinfo/VCFinfo/ALL.vcf.gz --type GT  --target $i  ; done > HaplotypeDiversity/CLM.hap.txt

# MXL
for i in $(seq -s, 631 694); do /Users/Pedro/Documents/Sickle\ Paper/Original\ Files/updated/Haplotypes/vcflib.source/sequenceDiversity --file  /Volumes/PEDRO\ 70GB/Genetics/ALLinfo/VCFinfo/ALL.vcf.gz --type GT  --target $i  ; done > HaplotypeDiversity/MXL.hap.txt

# PEL 
for i in $(seq -s, 546 630); do /Users/Pedro/Documents/Sickle\ Paper/Original\ Files/updated/Haplotypes/vcflib.source/sequenceDiversity --file  /Volumes/PEDRO\ 70GB/Genetics/ALLinfo/VCFinfo/ALL.vcf.gz --type GT  --target $i  ; done > HaplotypeDiversity/PEL.hap.txt


# awk '{ total += $4; count++} END{print total/count}' BR.hap.txt 
awk '{ total += $4; count++} END{print total/count}' BR.hap.txt #2.32128e-05
awk '{ total += $4; count++} END{print total/count}' PUR.hap.txt #2.25745e-05
awk '{ total += $4; count++} END{print total/count}' CLM.hap.txt #2.22268e-05
awk '{ total += $4; count++} END{print total/count}' MXL.hap.txt #2.15263e-05
awk '{ total += $4; count++} END{print total/count}' PEL.hap.txt #1.98863e-05



# vcflib in /Users/Pedro/Documents/Sickle\ Paper/Original\ Files/updated/Haplotypes/vcflib.source/          

#finestructure markers used by Galina
 wc -l /Volumes/PEDRO\ 70GB/Genetics/final/EIG-master/EIGENSTRAT/galina.snp # 302369
 
# finding out number of SNPs used without LD pruing 
# Yields 237982, minus the number os flipping SNPs: 237982-44123 = 193,859 (close to 196,390, probably used different flip list)
plink --bfile ../../BRn171-nam-eur-afr-amr.ordered --maf 0.05 --geno 0.05  --hwe 0.000000001 --write-snplist --out total.snps.txt

# For LD decay: 196749 (From BR.LD-decay.log - see gDrive) - so no LD filtering 

#####======================== ROH analysis ================================
##  Based on Galina's "ROHanalysis" scipt
##=======================================================
##=======================================================
#####=============================================================================================================
## See Ceballos 2018, Howrigan 2011 (see table 3) and Meyermans 2020 on the chose of parameters

## Brazilian
# 192987 variants and 171 people pass filters and QC.
#--homozyg: Scan complete, found 2082 ROH.
plink --bfile BR_n171.thin6 \
        --homozyg group \
        --homozyg-snp 50 \
        --homozyg-kb 1000 \
        --homozyg-gap 1000 \
        --homozyg-window-snp 50 \
        --homozyg-window-missing 5 \
        --homozyg-window-threshold 0.05 \
        --homozyg-het 0 \
        --out BR_n171.thin6.1Mb.ROHs

## All Others
# for pop in $(echo ACB ASW CEU CLM ESN FIN GBR GWD IBS LWK MSL MXL NAM PEL PUR TSI YRI) ; do
for pop in $(ls ~/Dropbox/input_chr16/Pops/ | cut -d . -f1 | uniq); do \  
plink --bfile ../BRn171-nam-eur-afr-amr.ordered \
        --extract BR_n171.thin6.bim \
        --keep ~/Dropbox/input_chr16/Pops/$pop \
        --homozyg group \
        --homozyg-snp 50 \
        --homozyg-kb 1000 \
        --homozyg-gap 1000 \
        --homozyg-window-snp 50 \
        --homozyg-window-missing 5 \
        --homozyg-window-threshold 0.05 \
        --homozyg-het 3 \
        --out $pop.thin6.1Mb.ROHs; 
done

# check general counts:
for i in $(ls *.log); do grep "Scan" $i && echo $i  ; done
# BR_n171.thin6.1Mb.ROHs.log
--homozyg: Scan complete, found 2082 ROH.
#BR_n171.thin6.300kb.ROHs.2.log
#--homozyg: Scan complete, found 10406 ROH.
#BR_n171.thin6.ROHs.log
#--homozyg: Scan complete, found 10303 ROH.
#CEU.thin6.1Mb.ROHs.log
#--homozyg: Scan complete, found 2718 ROH.
#CLM.thin6.1Mb.ROHs.log
--homozyg: Scan complete, found 2555 ROH.
#ESN.thin6.1Mb.ROHs.log
#--homozyg: Scan complete, found 1184 ROH.
#FIN.thin6.1Mb.ROHs.log
#--homozyg: Scan complete, found 3512 ROH.
#GBR.thin6.1Mb.ROHs.log
#--homozyg: Scan complete, found 2727 ROH.
#GWD.thin6.1Mb.ROHs.log
#--homozyg: Scan complete, found 1096 ROH.
#IBS.thin6.1Mb.ROHs.log
#--homozyg: Scan complete, found 2825 ROH.
#LWK.thin6.1Mb.ROHs.log
#--homozyg: Scan complete, found 1137 ROH.
#MSL.thin6.1Mb.ROHs.log
#--homozyg: Scan complete, found 726 ROH.
#MXL.thin6.1Mb.ROHs.log
--homozyg: Scan complete, found 2342 ROH.
#NAM.thin6.1Mb.ROHs.log
--homozyg: Scan complete, found 4575 ROH.
#PEL.thin6.1Mb.ROHs.log
--homozyg: Scan complete, found 5551 ROH.
#PUR.thin6.1Mb.ROHs.log
--homozyg: Scan complete, found 2190 ROH.
#TSI.thin6.1Mb.ROHs.log
#--homozyg: Scan complete, found 2689 ROH.
#YRI.thin6.1Mb.ROHs.log
#--homozyg: Scan complete, found 1167 ROH.



# Plot with the roh.R scipt


##=======================================================

#####======================== Diversity ================================
##  See https://www.cog-genomics.org/plink/1.9/basic_stats
##=======================================================
##=======================================================
#####========================================================================

for pop in $(echo ACB ASW BR CEU CLM ESN FIN GBR GWD IBS LWK MSL MXL NAM PEL PUR TSI YRI) ; do 
    plink --bfile ../BRn171-nam-eur-afr-amr.ordered --extract ../ROH/BR_n171.thin6.bim --keep ~/Dropbox/input_chr16/Pops/$pop --hardy  --out $pop  ; done

# filter p-value ($9) and average expected heterozygosity
awk 'NR>1{ total += $8; count++  }END{ print total/count}' BR.hwe # 0.383972
awk 'NR>1{ total += $8; count++  }END{ print total/count}' PEL.hwe  #   0.3505
awk 'NR>1{ total += $8; count++  }END{ print total/count}' PUR.hwe # 0.378918
awk 'NR>1{ total += $8; count++  }END{ print total/count}' MXL.hwe  # 0.372061
awk 'NR>1{ total += $8; count++  }END{ print total/count}' CLM.hwe  # 0.377863

# for Haplotype diversity:
cat ../input_chr16/Pops/BR.txt ../input_chr16/Pops/PUR.txt ../input_chr16/Pops/CLM.txt ../input_chr16/Pops/MXL.txt ../input_chr16/Pops/PEL.txt > ../input_chr16/Pops/latinos.txt

bcftools view -Ov -M2 -v snps -S ../input_chr16/Pops/latinos.txt ALL.vcf > Diversity/HaplotypeDiversity/DnaSP/latinos.vcf

# extract haplotypes
bcftools convert --hapsample --vcf-ids latinos.vcf -o latinos.haplotypes.vcf

# Use DnaSP (Windows only): File > Multi-MSA Data file analysis (SNP positions, *.vcf)...





##=======================================================
