##======================== 17. Finestructure/Chromosome "Painting" ==============================# 
############################# In progress #############################
######## Important note: nowadays finestructure pipeline (fs + cp) can be run in a single wrapper function, see https://people.maths.bris.ac.uk/~madjl/finestructure/manualse2.html#x5-40002
## Install finestructure:
## This is a really helpfull package to get inputs to ChromoPainter (and globetrotter):
## https://github.com/maarjalepamets/human-admixture - especially the "create_population_list_infile_and_idfile.py"
## script for creating the *optional* "-f switch" input file ("donor_list_infile" - see Instruction Manual for “ChromoPainter : a copying model for exploring admixture in population data ”)
## ˆ Isn't really necessary, as -a flag can scan target individuals against everyone else
## I performed these steps by taking a look in these major manuals: 
##      Instruction Manual for “ ChromoPainter : a copying model for exploring admixture in population data ” (Mendeley)
##      Practical Guide to the Analytical Pipeline (Mendeley)
##      https://people.maths.bris.ac.uk/~madjl/finestructure-old/data_example.html (first section)
##      ~/Documents/Sickle Paper/Original Files/updated/Haplotypes/badMIXTUREexample-master/convert.sh (from the human-admixture pkg above mentioned) 
##      /Users/Pedro/Documents/Sickle\ Paper/Original\ Files/updated/Haplotypes/human-admixture/Lepamets_admixture_poster.pdf ## See the "chromopainter.sh" for reference on 
## this "human-admixture" pipeline (https://github.com/maarjalepamets/human-admixture)
##=======================================================================================# 
####################################  Prelude ##############################################
# Get interesting HGDP samples that might help us discuss the data > http://hgdp.uchicago.edu/Phased_data/
grep -v "^\d" /Users/Pedro/Downloads/chrom11_hapguess_switch.out # labels > save in HGDP.labels.txt

# Duplicate each row
awk -v num=1 'BEGIN {OFS=FS=" "} {tmp=$2; print; for (i=1;i<=num;i++) {$2=tmp"_"i; print}} ' ~/Downloads/HGDP.labels.txt | awk '{print $1} '  > ~/Downloads/HGDP.labels.to.ped.txt 

# Assign ind name to first column of each haplotype
grep  "^\d" /Users/Pedro/Downloads/chrom11_hapguess_switch.out > chrom11_hapguess_switch.to.ped
paste -d" "  HGDP.labels.to.ped.txt  chrom11_hapguess_switch.to.ped > chrom11_hapguess_switch.labeled.to.ped 
 
# Unite both haplotypes by individual. The resulting file is very similar to plink PED (lacking cols 2-6)
awk  '                                                                                                                     function p(n,A) {
                     s = n
                     for (i=2; i<=NF; i++)  {s = s FS  A[i]
                                             A[i] = $i
                                            }
                     if (n) print s
                    }


    $1==n           {
    for (i=2; i<=NF; i++)  A[i] = A[i] "|" $i
                     next
                    }

                    {p(n,A)
                     n = $1
                    }

    END             {p(n,A)
                    }
    ' ~/Downloads/chrom11_hapguess_switch.labeled.to.ped  > ~/Downloads/chrom11_hapguess_switch.united.to.ped 
    
# Transpose (rows > columns). Resulting file ~ to vcf except it lacks #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT 
# cols and also a header, as the following one:
##fileformat=VCFv4.2
##fileDate=20181130
##source=PLINKv1.90
##contig=<ID=11,length=134942627>
##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}'  chrom11_hapguess_switch.united.to.ped |  sed 's/  */     /g'  > chrom11_hapguess_switch.united.to.vcf # this sed is ctrl+v and hit tab key (substitute any numeber of spaces: 2 spaces and *)


# Add vcf fields
awk 'BEGIN{print "#CHROM" "\t" "POS" "\t" "ID" "\t" "REF" "\t" "ALT" "\t" "QUAL" "\t" "FILTER" "\t" "INFO" "\t" "FORMAT"};
NR>1{print "11" "\t" $4 "\t" $1"\t" $5"\t" $6"\t" "." "\t" "." "\t" "." "\t" "GT" }; ' chrom11.final | paste - chrom11_hapguess_switch.united.to.vcf > HGDP.chr11.vcf

# getting a pgen/pvar/psam set from vcf preserving phase
plink2 --vcf HGDP.chr11.vcf.gz --export vcf --out HGDP.plink.chr11 

# HGDP has fucked up marker positions, change to match my vcf (when the markers match, otherwise will keep the positions)
awk  ' 
    BEGIN {OFS=FS="\t"} 
    NR==FNR{a[$3]=$2; next} NF>6{if($3 in a) {
        $2=a[$3];}else{$2=$2}}1
    '  ALL.NOT.IMPUTED.BUT.PHASED.vcf HGDP.plink.chr11.vcf > HGDP.final.chr11.vcf

# This makes the resulting vcf to be unsorted, so do 
cat HGDP.final.chr11.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > HGDP.final.sorted.vcf 

bgzip -c HGDP.final.sorted.vcf > HGDP.final.sorted.vcf.gz
tabix -p vcf  HGDP.final.chr11.vcf.gz 


#isec_output/0000.vcf.gz variants unique to 1.vcf.gz
#isec_output/0001.vcf.gz variants unique to 2.vcf.gz
#isec_output/0002.vcf.gz variants shared by 1.vcf.gz and 2.vcf.gz as represented in 1.vcf.gz
#isec_output/0003.vcf.gz variants shared by 1.vcf.gz and 2.vcf.gz as represented in 2.vcf.gz
bcftools isec -p isec_output -Ozvcf HGDP.final.sorted.vcf.gz ALL.NOT.IMPUTED.BUT.PHASED.vcf.gz 

# Go to line 988 to make a haplotype network out of them
bcftools merge --merge all ./isec_output/0002.vcf.gz ./isec_output/0003.vcf.gz > ALL.HGDP.vcf
    
###################### (A) Download page gives compiled binaries (install_fs.sh just chose fs and set paths)##############
# It might be interesting to compile from source and try to make openMP work, as it is required to multithreading/
# In fs-2.1.3/
###############################################################################################################################
./configure CC=/usr/local/opt/llvm/bin/clang CXX=/usr/local/opt/llvm/bin/clang++  
make install
# Use like this: ~/Documents/Sickle Paper/Original Files/updated/CPainting/fs-2.1.3/fs -h

####################### (B) Preparing (original) inputs ###############################################################
plink --double-id  --vcf ALL.NOT.IMPUTED.BUT.PHASED.vcf  --recode12 --out ALL.to.ChromoPainter

plink --file /Users/Pedro/Documents/Sickle\ Paper/Original\ Files/updated/Haplotypes/ALL.to.ChromoPainter --remove BRZ.txt --recode12 --out /Users/Pedro/Documents/Sickle\ Paper/Original\ Files/updated/Haplotypes/ALL.to.ChromoPainter

perl ./fs-2.1.3/scripts/plink2chromopainter.pl -p ALL.to.ChromoPainter.ped -m ALL.to.ChromoPainter.map -o ALL.to.ChromoPainter -d IDS.to.ChromoPainter

# Transfer to broadsword server (only one I have available right now)
scp -P 2222 /Users/Pedro/Documents/Sickle\ Paper/Original\ Files/updated/Haplotypes/IDS.to.ChromoPainter pcruz@143.106.4.111:/home/pcruz/CPainting 
scp -P 2222 /Users/Pedro/Documents/Sickle\ Paper/Original\ Files/updated/Haplotypes/ALL.to.ChromoPainter pcruz@143.106.4.111:/home/pcruz/CPainting &
scp -P 2222 /Users/Pedro/Dropbox/Documentos/OneDrive/Projeto_Doutorado/Análise/2013.04.03_plink/Manhattan\ Plot/genetic_map_b37/genetic_map_chr11_combined_b37.txt   pcruz@143.106.4.111:/home/pcruz/CPainting 

# ====== Now in server ========
# previously: mkdir CPainting
mv ALL.to.ChromoPainter ALL.to.ChromoPainter.phase 

#### OPTION 1 Getting a recombination file with the same recombination rate all over the chromosome (I have been using this, but the next one I guess is better)
perl ../fs-2.1.3/scripts/makeuniformrecfile.pl ALL.to.ChromoPainter.phase filep.recomb
###

#### OPTION 2 Using the hapmap genetic map available do
sudo apt-get install libswitch-perl # for the 'cant locate switch.pm' issue
while read A B C ; do echo chr11 $A $B $C ; done < genetic_map_chr11_combined_b37.txt > genetic_map_chr11_b37.txt # change the header in vi (just add name to 1st col)

# convertrecfile.pl script requires same markers in files so do
for i in $(awk ' NR==4 ' ALL.to.ChromoPainter.phase.a.option); # reads markers pos in input file
    do grep -w  $i genetic_map_chr11_b37.txt ;      # look for each one in the genetic map and print the line
done > genetic_map_chr11.txt

# Finally get the final chr11.ALL.recombfile with actual recombination rates
./convertrecfile.pl -M hapmap ALL.to.ChromoPainter.phase.a.option genetic_map_chr11.txt chr11.ALL.recombfile

# OPTIONS: -I snp -U centimorgans -M cdf -T norm -t 2,4 -s
# Effective unit= 0.01
# FILES: PHASE FILE: ALL.to.ChromoPainter.phase.a.option, REC FILE: genetic_map_chr11.txt -> REC FILE: chr11.ALL.recombfile
# Detected Chromopainter v1 format
# Detected 1737 individuals
# And 7324 SNPs

#### 
./convertrecfile.pl -M hapmap ALL.to.ChromoPainter.phase.a.option genetic_map_chr11_combined_b37.txt chr11.ALL.recombfile # gives an error
 ## ERROR: No valid lines found in the location of SNPs. Something is wrong
 
  plink2chromopainter.pl is a nice script, but first lines are not in agreement to CP:
# • The The first line of the file contains the number of donor haplotypes. (NOTE: if the ’-a’ switch is specified, this first line should be 0.)
# • The second line of the file contains the total number of donor and recipient individuals. If the ’-j’ switch is used, indicating all individuals are hap- loid, this will be the total number of donor and recipient haplotypes in the file. If the ’-j’ switch is not specified, indicating individuals are diploid, this number should be half the total number of donor and recipient haplo- types in the file. (In this latter caes, note that unless the ’-a’ switch is used, only recipient individuals are assumed to be diploid. Thus each donor pop- ulation can have any number of haplotypes representing it – the program does not ever assume that any given pair of donor haplotypes forms an individual. For this reason, the number entered here can be a fraction; for example if there are 7 total donor haplotypes and 3 diploid recipient individuals, the number entered here should be “6.5”(=[7+3*2]/2).)
# • The third line contains the number of SNPs. 3
# • The fourth line contains the letter “P” followed by a vector of the basepair positions of each SNP, in monotonically increasing order. The basepair positions do not strictly need to be in monotonically increasing order if you are including genetic information from multiple chromosomes, as specified in recom rate infile below. However, within each chromosome, basepairs must be in order.
# • The fifth line contains a vector of “S”s, with one “S” per SNP (this line is ignored in the current implementation, but some value must still be specified in line 5).
# • The remaining lines of the file contain the genetic variation information of each donor and recipient haplotype, with each row a new haplotype and each column the allelic type at each biallelic SNP, in the same order as the “positions” line. The accepted allelic type values are “0”, “1”, “A”, “G”, “C”, or “T”. There should be NO missing values!! There should be no spaces between columns for the genotype rows. Donor haplotypes are listed in the initial rows of haplotypes, followed by recipient haplotypes in the final rows of the input file. If individuals are diploid, each pair of 2 contiguous rows should be the two haplotypes for a single individual.
## Let's solve this issue by doing
# also, sed command in Mac os X behaves differently and so let's do the following in Linux (GNU sed), see https://stackoverflow.com/questions/6537490/insert-a-line-at-specific-line-number-with-sed-or-awk
sed '1i0' ALL.to.ChromoPainter | sed -e '2s/3474/1737/'  | sed "5i$(printf 'S%.0s' {1..7324})" > ALL.to.ChromoPainter.phase.a.option

#                                             ˆ couldn't find a way to divide value in second line (3474) by 2 

# Should be # chromosomes (ind*2) + 5 (header lines)
wc -l ALL.to.ChromoPainter.phase.a.option # 3479, OK

###########====================== (C) Running cp alone =========================
## In the "FinestructureExample.R" line 177: "chromopainter -a 1 1 -b -in -iM -i 10 -g EastAsiaSimple.chrom1.phase -r EastAsiaSimple.chrom1.trecombfile -o EastAsiaSimple.chrom1.linked.hap1"
## In the "/Users/Pedro/Documents/Sickle\ Paper/Original\ Files/updated/Haplotypes/human-admixture/pipeline/chromopainter.sh":
## ./ChromoPainterv2 -a 0 0 -i 10 -in -iM -s 0 -g ${1}${2}.haplotypes -r ${1}${2}.recomrates -t ${1}${2}.idfile -f ${1}${2}.poplist $n $n -o $1EMest/${2}.$n > $1EMest/log.$n &
## ./ChromoPainterv2 -s 10 $necmd -g ${1}${2}.haplotypes -r ${1}${2}.recomrates -t ${1}${2}.idfile -f ${1}${2}.poplist 0 0 -o ${1}${2} > ${1}${2}.log
# My first ever run: chromopainter -g ALL.to.ChromoPainter.phase.a.option -r filep.recomb -in -i 10 -a 134 226
#################################################################################

# estimating parameters
chromopainter -a 0 0 -i 10 -in -iM -s 0 -g ALL.to.ChromoPainter.phase.a.option -r filep.recomb -t IDS.to.ChromoPainter -o ALL.to.ChromoPainter.phase.a.option >ALL.to.ChromoPainter.phase.a.option.log

# calculating final parameters for the actual run (neavarage is part of tools in page https://people.maths.bris.ac.uk/~madjl/finestructure-old/toolsummary.html and human-admixrure pkg, but it does simply the average of the estimates, I'd be able to do in bash)
./neaverage.pl -o ALL.neaverage.txt -v  ALL.to.ChromoPainter.phase.a.option.EMprobs.out

necmd=$(cat ALL.neaverage.txt)

# actual ChromoPainter run <<<--- this one is in "Haplotypes" folder in Mac
chromopainter -a 134 226 -s 10 $necmd -g ALL.to.ChromoPainter.phase.a.option -r filep.recomb -t IDS.to.ChromoPainter -b -o SBR > SBR.log  #around 4 hours in broadsword

# Will condition each individual on every other individual...
#  Number of EM-runs = 0
#  Number of samples = 10
#  N_e value = 31.948488
#  Region size = 100.000000
#  Global mutation value = 0.005311
#  Number of donor haplotypes = 0
#  Number of recipient haplotypes = 3474
#  num donor pops = 1736

# Or with actual recombination rates from HapMap <<<--- this one is in "CPainting" folder in Mac
chromopainter -a 134 226 -b -g ALL.to.ChromoPainter.phase.a.option -r chr11.ALL.recombfile -t IDS.to.ChromoPainter -o SBR > SBR.log

# Will condition each individual on every other individual...
#  Number of EM-runs = 0
#  Number of samples = 10
#  N_e value = 115.141048
#  Region size = 100.000000
#  Global mutation value = -9.000000
#  Number of donor haplotypes = 0
#  Number of recipient haplotypes = 3474
#  num donor pops = 1736                    

# Final on SBR
chromopainter  -a 134 226 -s 10 $necmd -g ALL.to.ChromoPainter.phase.a.option -r chr11.ALL.recombfile -t IDS.to.ChromoPainter -b -o SBR2 > SBR2.log

# For all subjects
chromopainter -a 0 0 -s 10 $necmd  -g  ALL.to.ChromoPainter.phase.a.option -r chr11.ALL.recombfile -t IDS.to.ChromoPainter -b -o ALL
#  Region size = 100.000000
#  Global mutation value = 0.005311
#  Number of donor haplotypes = 0
#  Number of recipient haplotypes = 3474
#  num donor pops = 1736

# rownames (recipients) do not match SBR samples (lines 134-226), do
#cat SBRtags.chunkcounts.out | awk 'NR==2 {print $0}' | tr ' ' '\n' > right.rownames.cp.txt > SBR.righttags.regionchunkcounts.out
#awk 'NR==FNR{a[NR]=$1;next} {$1=a[FNR]}1' right.rownames.cp.txt SBRtags.chunkcounts.out > SBR.righttags.chunkcounts.out
#awk 'NR==FNR{a[NR]=$1;next} {$1=a[FNR]}1' right.rownames.cp.txt SBRtags.chunklengths.out > SBR.righttags.chunklengths.out
#awk 'NR==FNR{a[NR]=$1;next} {$1=a[FNR]}1' right.rownames.cp.txt SBRtags.mutationprobs.out > SBR.righttags.mutationprobs.out
#awk 'NR==FNR{a[NR]=$1;next} {$1=a[FNR]}1' right.rownames.cp.txt SBRtags.regionchunkcounts.out > 
#awk 'NR==FNR{a[NR]=$1;next} {$1=a[FNR]}1' right.rownames.cp.txt SBRtags.regionsquaredchunkcounts.out > SBR.righttags.regionsquaredchunkcounts.out

#Not sorting due to name differences between rows and columns.
#This could mean that your data has been read incorrectly!
#Some suggestions are:
#1. Did you split the data into many files? Check that all the files you expected to be created by chromopainter exist, and have the correct number of lines. You may wish to try rerunning stating the input files explicitly, in case autodetection of files has failed.
#2. Are there file line-ending issues? (did you run chromocombine
# on a Windows/UNIX/Mac machine and chromopainter elsewhere?)
#3. Alternatively, you may be using chromopainter in donor population mode, in which case this warning is expected.
#Continuing, in case you intended this. Expect problems!
#Successfully summed 1 file root(s) containing 93 individuals, to new file root SBRtags.combined
#Successfully written chunkcounts file SBRtags.combined.chunkcounts.out with c value 0.238063
#When using these files with fineSTRUCTURE there is now no need to specify "c".
#ChromoCombine completed successfully
~/Documents/Sickle\ Paper/Original\ Files/updated/CPainting/fs_4.1.1/fs combine -i indivnamelist.lineendings.txt  -o SBRtags.combined SBRtags


# Name tags for individuals (<POP><number>)
perl chromopainterindivrename.pl indivnamelist.txt SBR2 SBRtags


export PATH="/Users/Pedro/Documents/Sickle\ Paper/Original\ Files/updated/CPainting/fs_4.0.1:$PATH"


#./fs_mac finestructure -x 100000 -y 100000 -z 100 SBRtags.combined chunkcounts.out SBR.mcmc.xml
ˆ#ERROR: row and column name 4 are incompatible!
#row name=NativeAmerican2, col name=NativeAmerican5
#Writing data to SBRtags.chunkcounts.out.testdata check it manually!
#This is the symptom if you have the wrong line endings.
#Try converting them to that of your platform (or UNIX are usually handled correctly everywhere)
#Alternatively, you may not be using a square matrix.
#Data error

# For CP alone (after running section D below)
awk '{ print $4}' IDS.to.pops.lineendings.txt > IDS.to.ChromoPainter.2020.txt

#chromopainter -a134 226 -s 10 -n 31.9484876524292 -M 0.0053106644375 -g ALL.to.ChromoPainter.phase.a.option -r chr11.ALL.recombfile -tIDS.to.ChromoPainter -b -o SBR2
# It took ~1h30 in my macbook
./fs-2.1.3/fs cp -a 134 226 -s 10 -n 31.9484876524292 -M 0.0053106644375 -g ALL.to.ChromoPainter -r filep.recomb -t IDS.to.ChromoPainter.2020.txt -b -o SBR2020 

############################### (D) Trying to run whole pipeline through fs wrapper ########################
# By looking to SBR Ive end up with a not squared matrix, finesrtucture alone can't work with it. 
# Lets try to run the wrapper to the whole pipeline and then test the mcmc and tree in 'R library'
############################################################################################################
# Because the downloaded latest version (fs_4.1.1) is not working with multithreading, lets compile from source (https://people.maths.bris.ac.uk/~madjl/finestructure/fs-2.1.3.tar.gz)
# Follow (just once)
# See https://ryanhomer.github.io/posts/build-openmp-macos-catalina and
# https://unix.stackexchange.com/questions/149359/what-is-the-correct-syntax-to-add-cflags-and-ldflags-to-configure
cd ./fs-2.1.3/
./configure CC=/usr/local/opt/llvm/bin/clang CXX=/usr/local/opt/llvm/bin/clang++  LDFLAGS="-L/usr/local/opt/llvm/lib"
make

#What do I do if I'm interested in some individuals more than others? The best thing to do is still to use fs2 to paint ALL INDIVIDUALS AGAINST ALL OTHERS. You should then construct a `finestructure force file' which will allow you to manually define populations within the `less interesting' indivduals. The individuals of interest can then be clustered by FineSTRUCTURE using the information about differential rates of ancestry sharing with these populations, usually leading to increased power relative to either running the whole sample simultaneuously or leaving out some individuals. R users can check out the (documented) R scripts I used to generate the HGDP results.
cd .. # i.e. ~/Documents/Sickle Paper/Original Files/updated/CPainting
./fs-2.1.3/fs  ALL2.cp -idfile IDS.to.ChromoPainter -phasefiles ALL.to.ChromoPainter -recombfiles filep.recomb -go

# Estimated params:
# Successfully run ChromoCombine stage! Inferred a 'c' value of 0.116812
#N_e value = 56.242100 
# Global mutation value = 0.003431
# Inferred a 'c' value of 0.116812

# Stopped in the beggining of stage3 (after combine) because of individual names starting with number solve this by
awk 'NR==FNR{a[NR+2]=$1;next} NR>2{$1=a[FNR]}1' indivnamelist.lineendings.txt ALL2_linked.chunkcounts.out.old >  ALL2_linked.chunkcounts.out  # fix the 2 first fields of column 1 by hand 
awk 'NR==FNR{a[NR+1]=$1;next} NR>2{$1=a[FNR]}1' indivnamelist.lineendings.txt ALL2_linked.mutationprobs.out.old >  ALL2_linked.mutationprobs.out # fix first fields of column 1 by hand 
awk 'NR==FNR{a[NR+1]=$1;next} NR>2{$1=a[FNR]}1' indivnamelist.lineendings.txt ALL2_linked.regionchunkcounts.out.old >  ALL2_linked.regionchunkcounts.out # fix first fields of column 1 by hand 
awk 'NR==FNR{a[NR+1]=$1;next} NR>2{$1=a[FNR]}1' indivnamelist.lineendings.txt ALL2_linked.chunklengths.out.old >  ALL2_linked.chunklengths.out # fix first fields of column 1 by hand 
awk 'NR==FNR{a[NR+1]=$1;next} NR>2{$1=a[FNR]}1' indivnamelist.lineendings.txt ALL2_linked.regionsquaredchunkcounts.out.old >  ALL2_linked.regionsquaredchunkcounts.out # fix first fields of column 1 by hand

# Resume
 ./fs-2.1.3/fs ALL2.cp -go
 
################################# (E)  Plotting #########################################
####################################
setwd("/Users/Pedro/Documents/Sickle Paper/Original Files/updated/CPainting/FinestructureRcode/")
source("/Users/Pedro/Documents/Sickle Paper/Original Files/updated/CPainting/FinestructureRcode/FinestructureLibrary.R") # read in the R functions, which also calls the needed packages
setwd("/Users/Pedro/Documents/Sickle Paper/Original Files/updated/CPainting/")
chunkfile<-"SBR2.chunkcounts.out"
dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=1))
copyprobsfile<-"SBR2.copyprobsperlocus.out.gz"

# 267 is the first haplotype in the sample 134 (134*2 - 1), see chromopainter run in line 1339
myhap<-getHap(267,copyprobsfile,verbose=TRUE)

# Plotting first 1000 SNPs
cpdensityplot(myhap$snps[1:1000],myhap$probs[1:1000,],simplecollist)

# Plotting from 4,500,031 to 5,777,053 (open ALL.to.ChromoPainter.map in excel to look at coordinates)
simplecollist<-MakeColorYRP(0.1) 
cpdensityplot(myhap$snps[197:321],myhap$probs[197:321,],simplecollist)s