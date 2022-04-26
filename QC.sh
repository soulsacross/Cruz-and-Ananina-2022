# --- Pedro Cruz --- #
#    Script for general population analysis from PED/BED files


# --- 1. Getting VCF from 1000 Genomes and converting to plink format --- #
#    Start here if you are fresh starting
#
# 


#======================================
## Download files from the site:
## ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/bcf_files/
## Or get it through DataSlicer: http://grch37.ensembl.org/Homo_sapiens/Tools/DataSlicer
#======================================

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/bcf_files/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf

# For Native Americans
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130711_native_american_admix_train/native_amr_train_20130711.*

#======================================
## Convert to vcf format applying some filters
## bcftools options
## -Ov,,output file in vcf format
## -M2,,keep only markers with at most 2 alleles
## -v snps,extract snps variants only
## -S [file],,use file containing the samples list


#======================================

## Extract samples ID for the desired populations:
grep -e EUR -e AMR -e AFR 1000genoms-release-20130502-integrated_call_samples.txt | cut -f 1 > Eur.Amr.Afr.txt

## Filter and convert to vcfs (I had som trouble to convert from bcf to plink format however in theory it is possible)
bcftools view -Ov -M2 -v snps -S Eur.Amr.Afr.txt ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf > filtered/eur-afr-amr.chr22.M2.snps.only.vcf

#======================================
## Convert vcf to plink format
#======================================

plink --vcf eur-afr-amr.chr22.M2.snps.only.vcf --double-id --keep-allele-order --vcf-half-call m --out plink_format/chr22.eur-afr-amr

#======================================
## extract affy6.0 markers only
#======================================

plink --bfile plink_format/chr22.eur-afr-amr --extract Affy6.snps.list --make-bed --out plink_format/chr22.eur-afr-amr.aff6



# --- 2. Primary Quality Control --- #
#    Start from here if your input are plink PED/BED files raw data (Quality Control check needed)
#
#

### Preliminary stuff

# Take a look at the directory tree
#ls -R | grep ":$" | sed -e 's/:$//' -e 's/[^-][^\/]*\//--/g' -e 's/^/   /' -e 's/-/|/' 
#    |-para_Pedro > Old files
#    |-updated    > New files 
#    |---Poplists

# Get the individuals list (by population) in the proper plink input format (two columns)
for i in $( ls ./Poplists/)  ; 
   do awk 'BEGIN { FS="\t"; OFS="\t" } { $1=$1 "\t" $1 }1' $i  > $i; 
done


# Assuming a path named ./Poplists/ with the ind lists (and only them) is present. Note that HWE is 0.001 by 859106 markers 
for i in $( ls ./Poplists/) ;
 do  plink --bfile ALL2017 --keep ./Poplists/$i  --maf 0.01 --mind 0.1 --geno 0.1 --hwe 0.000000001 --out ./snplist/$i --write-snplist; 
 done 
# IBS strictly plink --bfile ../ALL2017 --keep ../Poplists/IBS.txt  --maf 0.01 --mind 0.1 --geno 0.1 --hwe 0.000000001 --chr 11 --make-bed --out IBS.chr11 --write-snplist
# BRZ (all) strictly plink --file ../ALLBrazilian --keep ../Poplists/ALL.BRAZILIAN.txt  --maf 0.01 --mind 0.1 --geno 0.1 --hwe 0.000000001 --chr 11 --make-bed --out  BRZ.chr11 --write-snplist
# For the old Populations (here I'm just getting a list of .ped files names [chopping ".ped" out] and then I'm calling plink iteratively for each ped/map set
for i in $( ls | sed 's/.\{4\}$//' | uniq ) ;  
   do  plink --file $i  --maf 0.01 --mind 0.1 --geno 0.1 --hwe 0.000000001 --out ../snplist$i --write-snplist;   
done 


# Sort all *.snplist files in directory (it may take a few minutes)
for file in *.snplist; 
   do sort -o $file $file ; 
done


# Get common list in all snplist file (477981 markers), then clean the directory from the bunch of snplists
awk '{$1=$1} ++A[$0]>=ARGC-1' *.snplist > ALL2017.snplist
rm -r snplist*

echo "=== end of QC step ==="


# --- 2. Merging files --- #
#    Start from here if your input is already filtered for population-specific parameters but you still need them to get merged
#
#

# In the first merge try it gave an error, so lets exclude discordant SNPs
grep -vf 2017ALL-merge.missnp ../ALL2017.snplist > to_mergeSNP.txt
plink --bfile ALL2017 --extract to_mergeSNP.txt --make-bed --out ALL2017.b 

# Finally merge
 plink --file AGCT.ALLCAG --allow-no-sex  --extract to_mergeSNP.txt --merge-list to_merge.txt --make-bed --out 2017ALL


echo "=== end of merging step ==="

# --- 3. Further quality filtering --- #
#    Start from here if you already have all populations you're working with on a single file and are starting
#    dealing with AT/GC SNPs and LD
#

# Getting a list of A-T and G-C SNPs
plink --bfile 2017ALL --freq --out 2017ALL  
awk '{print $1 "\t" $2 "\t" $3 $4 "\t" $5 "\t" $6}' 2017ALL.frq > 2017ALL.b.frq
grep -e AT -e TA -e GC -e CG 2017ALL.b.frq | awk '{print $2}'  > 2017ALL.AT-GC.snplist

echo "=== end of AT-GC flip removal step ==="


# LD handling
plink  --bfile 2017ALL_noflip --allow-no-sex --indep-pairwise 50 5 0.5 --out ALL.to && plink --bfile 2017ALL_noflip --allow-no-sex --extract ALL.to.prune.in --make-bed --out ALL.localLD\ 
&& plink --bfile ALL.localLD --allow-no-sex --indep-pairwise 1500 150 0.5 --out ALL.long.to && plink --bfile  ALL.localLD --allow-no-sex  --extract ALL.long.to.prune.in --make-bed --out ALLfinal  


# Remove sites of long LD in the human genome (hg19)
plink --bfile ALLfinal --allow-no-sex --chr 20 --from-mb 32 --to-mb 34 --write-snplist --out hr20.mb.32-34; plink --bfile ALLfinal --allow-no-sex --chr 12 --from-mb 109.5 --to-mb 112 --write-snplist --out chr12.mb.109.5-112; plink --bfile ALLfinal --allow-no-sex  --chr 12 --from-mb 33 --to-mb 40 --write-snplist --out  chr12.mb.33-40; plink --bfile ALLfinal --allow-no-sex --chr 11 --from-mb 87.5 --to-mb 90.5 --write-snplist --out chr11.mb.87.5-90.5; plink --bfile ALLfinal --allow-no-sex --chr 11 --from-mb 46 --to-mb 57 --write-snplist --out chr11.mb.46-57; plink --bfile ALLfinal --allow-no-sex --chr 10 --from-mb 37 --to-mb 43 --write-snplist --out chr10.mb.37-43; plink --bfile ALLfinal --allow-no-sex --chr 8 --from-mb 112 --to-mb 115 --write-snplist --out chr8.mb.112-115; plink --bfile ALLfinal --allow-no-sex --chr 8 --from-mb 43 --to-mb 50 --write-snplist --out chr8.mb.43-50; plink --bfile ALLfinal --allow-no-sex --chr 8 --from-mb 8 --to-mb 12 --write-snplist --out chr8.mb.8-12; plink --bfile ALLfinal --allow-no-sex --chr 7 --from-mb 55 --to-mb 66 --write-snplist --out chr7.mb.55-66; plink --bfile ALLfinal --allow-no-sex --chr 6 --from-mb 140 --to-mb 142.5 --write-snplist --out chr6.mb.140-142.5 ; plink --bfile ALLfinal --allow-no-sex --chr 6 --from-mb 25.5 --to-mb 35.5 --write-snplist --out chr6.mb.25.5-35 ; plink --bfile ALLfinal --allow-no-sex --chr 6 --from-mb 57 --to-mb 64 --write-snplist --out chr6.mb.57-64 ; plink --bfile ALLfinal --allow-no-sex --chr 5 --from-mb 135.5 --to-mb 138.5 --write-snplist --out chr5.mb.135.5-138.5; plink --bfile ALLfinal --allow-no-sex --chr 5 --from-mb 129 --to-mb 132 --write-snplist --out chr5.mb.129-132 ; plink --bfile ALLfinal --allow-no-sex --chr 5 --from-mb 98 --to-mb 100.5 --write-snplist --out chr5.mb.98-100.5 ; plink --bfile ALLfinal --allow-no-sex --chr 5 --from-mb 44.5 --to-mb 50.5 --write-snplist --out chr5.mb.44.5-50.5 ; plink --bfile ALLfinal --allow-no-sex --chr 3 --from-mb 89 --to-mb 97.5 --write-snplist --out chr3.mb.89-97.5; plink --bfile ALLfinal --allow-no-sex --chr 3 --from-mb 83.5 --to-mb 87 --write-snplist --out chr3.mb.83.5-87; plink --bfile ALLfinal --allow-no-sex --chr 3 --from-mb 47.5 --to-mb 50 --write-snplist --out chr3.mb.47.5-50;plink --bfile ALLfinal --allow-no-sex --chr 2 --from-mb 183 --to-mb 190 --write-snplist --out chr2.mb.183-190 ;plink --bfile ALLfinal --allow-no-sex --chr 2 --from-mb 134.5 --to-mb 138 --write-snplist --out chr2.mb.134.5-138 ;plink --bfile ALLfinal --allow-no-sex ;plink --bfile ALLfinal --allow-no-sex --chr 2 --from-mb 86 --to-mb 100.5 --write-snplist --out chr2.mb.86-100.5;plink --bfile ALLfinal --allow-no-sex --chr 1 --from-mb 48 --to-mb 52 --write-snplist --out chr1.mb.48-52

# concatenate snplists
cat *snplist > ALL.Long.LD.snplist  



echo "=== end of LD filtering: your dataset is now genotyping erros/LD free ==="


# --- 5. Plotting PCA and Ancestry (admixture) --- #
#    Start from here if you already have all populations filtered, LD thinned, AT/GC flip handleded and   
#    need basic population structure plotting
#

# Individuals order were messed up on the merging step, let's get the proper order
awk ' { print $1 " " $2 } ' ALL.indlist > OLD.indlist

cat FIN.txt GBR.txt  TSI.txt IBS.txt NAT.txt CEU.txt MXL.txt CLM.txt PEL.txt PUR.txt BRZ.txt SCB.txt SCA.txt  AAM.txt ASW.txt ESN.txt GWD.txt ACB.txt LWK.txt MSL.txt YRI.txt > ALL.txt

# Reorder individuals 
plink --allow-no-sex  --bfile ../ALL2017final  --indiv-sort f ALL.txt --make-bed --out ALLfinal

awk ' { print $1 "\t" $2 "\t" $3"\t"$4"\t"$5"\t"$6 } ' ../ALLfinal.fam > ALLindlist.txt # Put pop names in pheno column

# Run R to output ALL.indlist we are going to use in the next step
Rscript Population.Names.R

# Finally set final files for PCA and Admxiture plotting
plink --bfile ALLfinal --allow-no-sex --exclude  ALL.Long.LD.snplist --recode --tab --out ALL.Final && plink --file ALL.Final --allow-no-sex --recode12 --out ALL.final12coded

### Little fix: exclude people from Joinville-SC and test with and without Recife-PB
# 
#

plink --file ALL.Final  --remove ./Poplists/excludeJoinville.txt --allow-no-sex --recode12 --out ALL.final12coded # With Recife
plink --file ALL.Final  --remove ./Poplists/excludeJoinville_Recife.txt --allow-no-sex --recode12 --out ALL.final12coded2 # I decided for this one
# '->Total genotyping rate in remaining samples is 0.99331

perl plink-to-eigen-ped-converter_spacedelimited.pl
perl plink-to-eigen-ped-converter.pl

# __ALL.Final12coded.ped is ready for admixture analysis, so do for instance (test 20 K's)
for K in {1..20} ; 
    do /Users/Pedro/Dropbox/Documentos/OneDrive/Projeto_Doutorado/Análise/2013.04.03_plink/ADMIXTURE/admixture --j4 --cv __ALL.final12coded.ped $K | tee ALL.log${K}.out; 
done

# Pick the lowest error level
grep -h CV ALL.log*.out

# CV error (K=1): 0.64274
# CV error (K=2): 0.59542
# CV error (K=3): 0.58253
# CV error (K=4): 0.58200
# CV error (K=5): 0.58166
# CV error (K=6): 0.58153
# CV error (K=7): 0.58160
# CV error (K=8): 0.58173
# CV error (K=9): 0.58246

# __ALL.Final.ped and .map need go through CONVERTF subroutines
scp -P 22 __ALL.Final.* pedro@143.106.4.75:/home/public/Eigensoft_6.1.4/CONVERTF

# After editting "par.PED.EIGENSTRAT" (already done for this dataset) do in /home/public/Eigensoft_6.1.4/CONVERTF
perl example.perl # Edited to fit input files

# In the folder /EIGENSTRAT/ simply do
./ALL.Final.perl 

# For Fst alone, parameter file (ALL.Final.pca.par) should be edited, adding the following line:
# Try to remember to comment this line after running
fstonly: YES

# Then do 
./smartpca -p ALL.Final.pca.par > Fst.txt

# And also
 ./twstats -t twtable -i ALL.Final.eval > ALL.Final.stats.txt


#======================================
#
## 6. Problem solving
## Galina has found native americans (NAM from now on) to be messed on the SNPs order, lets remove, fix and merge again

#======================================

plink --bfile ALL2017final --remove NAT.fam --make-bed --out AlmostALL

# Run 20170723_bim-convert_forNatives.R

# Problem: I found NAs in chr (column 1, i.e., $1) after above script, get right values:

awk 'NR==FNR{a[$2]=$1} NR>FNR{ print a[$2] "\t" $2 "\t" $3 "\t" $4 "\t"  $5 "\t"  $6}'  AlmostALL.bim nam.bim > NAM2.bim # rename it to nam.bim

# Order is still fucked up, lets set the right markers' order
plink --bfile nam --make-bed --out nam2

plink --bfile nam --allow-no-sex  --extract to_mergeSNP.snplist --bmerge AlmostALL  --make-bed --out ALL2017.b 

plink --bfile ALL2017.b --allow-no-sex --recode12 --out ALL.final12coded

awk ' { print $1 "\t" $2 "\t" $3"\t"$4"\t"$5"\t"$6 } ' ../ALL2017.b.fam > ALLindlist.txt # Here you manually open on excel for instance and add labels

Rscript Population.Names.R

perl plink-to-eigen-ped-converter_spacedelimited.pl

# Finally get the admixture running

for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14; do  ../../public/admixture_1.3.0/admixture -s time  --cv __ALL.final12coded.ped $K -j22  | tee log${K}.out  ; done

grep -h CV  log*.out 

