## --- 12.2018 Pedro Cruz --- #
##   Script for general population analysis from PED/BED files
##   Intended for using in a "manual" manner by now, copying and running a line or a chunk at time.
##
##   This is not so computing-intensive, all steps were performed on a 8GB, 4 cores processor Mac OS X
##   I would nonetheless recommend running on a higher system, in 8GB linux should run smoothly. Don't try to do it in a Virtual
##   Machine with low disk and memory allocated        
##
##   Software needed:
##    bcftools
##    vcftools  
##    tabix  
##    plink 1.9
##    admixture
##    eigensoft
##    shapeit
##    ### Also in house scripts:###
##    plink-to-eigen-ped-converter.pl
##    plotPCA.R  
##    plotAdmixture.R


###########################################################
##
#==========================================================
##=== 1. Give population names and final conversions ======
#==========================================================
## 
###########################################################

## Plink cannot hold population names, we will have to set phenotype column (6th) to population names (ACB, NAM, YRI, etc) manually

## Include NAM and BDS in Poplists directory
cp BDS.to.merge.fam ./Poplists/BDS.txt
cp nam.to.merge.fam ./Poplists/NAM.txt

## This piece of code will get us the population each subject is derived from
awk 'NR==FNR { a[$2]=$1; next } $2 in a { print a[$2] "\t" $2 "\t" FILENAME; next }' ./LDflip/ALL.populations.LD.pruned.fam ./Poplists/*.txt | sed -e 's/\.\/Poplists\///g' | sed -e 's/\.txt//g' > PopulationTags 

## Now make a .fam file named ALL.indlist proper for the perl script in the next step 
awk 'BEGIN { FS=" "; OFS="\t" } NR==FNR{a[$1,$2]=$3;next}{$6=a[$1,$2]?a[$1,$2]:"."}1' PopulationTags ./LDflip/ALL.populations.LD.pruned.fam > ALL.indlist

##
## --- 3.1 Final format conversions to admixture/PCA  --- 
##     
##     

mkdir final && cd ./final

### Copy the scripts plink-to-eigen-ped-converter_spacedelimited.pl and plink-to-eigen-ped-converter.pl in the current directory.
## Copy also ALL.indlist containing population tags 
cp ../plink-to-eigen-ped-converter* .
cp ../ALL.indlist .

## Getting in format suited to admixture software 
## If running a 2nd time:  
##plink --bfile 1kg.to.merge --allow-no-sex --merge-list to_merge.txt --exclude ../LDflip/ALL.LD.flip.to.prune.out.snplist --recode12  --out ALL.Final12coded
plink --bfile ../LDflip/ALL.populations.LD.pruned --allow-no-sex --recode12 --out ALL.Final12coded

## This will run the script to get a .ped file with populations' names in it
perl plink-to-eigen-ped-converter_spacedelimited.pl

## Set map file to match the name of our new ped
mv ALL.Final12coded.map __ALL.Final12coded.map

## Remove the old .ped for saving disk memory
rm ALL.Final12coded.ped

##
##
## ----- __ALL.Final12coded.ped and __ALL.Final12coded.map are ready ----
## At this point you can skip the next fews lines and go on to admixture running (session 4)
## If not so much in a hurry, the following steps will analogously get EIGENSOFT inputs for our PCA (session 5)

## Getting in format suited to EIGENSOFGT 
plink --bfile ../LDflip/ALL.populations.LD.pruned --allow-no-sex --recode --tab --out ALL.Final 


## This will run the script to get a .ped file with populations' names in it
## plink-to-eigen-ped-converter.pl must be in current path, note that this is a perl script
## Although this works just fine, trying set population tags in bash is a good exercise, awk array should get the same result with one line of code 
perl plink-to-eigen-ped-converter.pl

## Set map file to match the name of our new ped
mv ALL.Final.map __ALL.Final.map

## Remove the old .ped for saving disk memory
rm ALL.Final.ped

##
## ----- __ALL.Final.ped and __ALL.Final.ped are ready ----
## This pair is the input to EIGENSOFT
##

###########################################################
##
#==========================================================
##=== 2. Running admixture to find out ancestries in dataset 
#==========================================================
## From now one sessions are basically independent, session
## 4 is not necessary to session 5 (PCA), one could run
## admixture and, since it takes a while, go on running 
## eigenstrat, which is quicker 
###########################################################

## Run admixture 
## Here I am testing 18 K to search for the best model fit (see admixture manual), i. e., the lowest error (cross-validation) 
## -j4 specifies 4 processing cores (all that I have in my computer), having more processors will shorten the process

mkdir "admixture Results"

## This may take some days in a standard machine
for K in {1..18} ;    
    do ./admixture -j4 --cv __ALL.final12coded.ped $K | tee ./admixture\ Results/ALL.log${K}.out;  
done

# After a while, take a look in error values generated in log files:
grep -h CV ./admixture\ Results/ALL**.out
# CV error (K=6): 0.59562
# CV error (K=7): 0.59542
# CV error (K=8): 0.59506
# CV error (K=9): 0.59497
# CV error (K=9): 0.58896
# CV error (K=10): 0.59518

# ˆ Look these values, the lowest CV error should be considered as the more likely true ancestral components number,
# in our case, K = 9. Read admixture documentation 

###########################################################
##
#==========================================================
##=== 3. Running EIGENSOFT to plot the PCA ================ 
#==========================================================
## 
###########################################################

## This chunk is for my own personal annotation purposes - ignore it
#### I finally could install eigensoft locally on my machine
## After checking with brew if openblas and gsl are installed, updated and linked, do
export LDFLAGS="-L/usr/local/opt/openblas/lib"
export CPPFLAGS="-I/usr/local/opt/openblas/include"
export PKG_CONFIG_PATH="/usr/local/opt/openblas/lib/pkgconfig"
make CFLAGS="-I/usr/local/opt/openblas/include -I/usr/local/opt/gsl/include" LDFLAGS="-L/usr/local/opt/openblas/lib -L/usr/local/opt/gsl/lib" # Probably running into an error
cd ../src/
make
make install

## Here you get back to paying attention 
## Start converting further, this time with EIGEN own routines
## To so such, start by going into CONVERTF directory
cd ./EIG-master/CONVERTF

## Open "par.PED.EIGENSTRAT" in text editor and insert 
## our "__ALL.Final.ped" instead of "example.ped" and so on
## Also rename output files to "ALL.Final" or something like this
## Also open "example.perl", comment "par.PED.EIGENSTRAT" and uncomment "par.PED.EIGENSTRAT"
## Then (it takes some 20min on a standard computer - a server do this in no time)
perl example.perl

## This gives us (or whatever prefix you gave):
ALL.Final.eigenstratgeno ALL.Final.ind and ALL.Final.snp 

mv ALL.Final* ../EIGENSTRAT/

## Let's now begin to actually run EIGENSOFT
cd ../EIGENSTRAT/

## There are actually many different ways to run eingenstrat, it is nice to take a look in the documentation
## 
## This step also require manual editting. Here in EIGENSTRSTRAT folder you'll find a file named example.pca.par
## Open it and replace "example.geno" by "ALL.Final.eigenstratgeno"
## Also replace "example.snp" "ALL.Final.snp"
## Moreover, replace "example.ind" "ALL.Final.ind"
## Optionally, give also names to output files, here I am choosing the ALL.Final prefix again
## Might, previously do
export PATH="/Users/Pedro/final/EIG-master/bin:$PATH"

## should work with
../bin/smartpca –p example.pca.par > ALL.Final.Sout

## I had no success with the above command, so let's do it in a more manually way
../bin/smartpca.perl -i ALL.Final.eigenstratgeno -a ALL.Final.snp  -b ALL.Final.ind -k 10 -o ALL.Final.pca -p ALL.Final.plot  -e All.Final.eval -l ALL.FInal.log -m 5 -t 10 -s 6.0

## I guess this took around 40 minutes to run
## There are some outputs (see documentation for details, explore the outputs beyond this script), 
## *.pca.evec file is the file from we get the PCA










