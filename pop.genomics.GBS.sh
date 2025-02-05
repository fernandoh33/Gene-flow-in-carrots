#estimate genetic diversity and divergence across the genome using pixy, the program should be installed locally in the server or personal laptop (good luck with that), here is the website: https://github.com/ksamuk/pixy
#once installed, the program is very simple to use and surprisingly fast, this is the code to estimate pi, dxy, and fst in 150 kb windows across the genome for populations indicated in the pops.file.txt (two columns file with SampleID PopID, without header)
#please note that pixy uses both invariant and variant (biallelic) sites, however, the data should be filtered for quality, missing rate, minimum depth, etc. For filtering while retaining invariant sites, you can use vcftools as follows:

vcftools --gzvcf raw.vcf.vcf.gz --minQ 30 --min-meanDP 4 --keep final.individuals.txt --max-missing 0.7 --max-alleles 2 --remove-indels --recode --stdout|bgzip -c > vcf.with.invariant.sites.vcf.gz

#then run pixy

pixy --stats pi,dxy,fst --vcf vcf.with.invariant.sites.vcf.gz --populations pops.file.txt --window_size 150000 --n_cores 16 --output_folder out_pixy/ --output_prefix out_pixy_150kb  

#Basic filtering for population genomics include minimum allele frequency and low linkage disequilibrium
#I use plink for filtering, the script is designed for compute canada but it should be similar for other servers

#!/bin/bash
#SBATCH --account=your_account
#SBATCH --time=0-4:00
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --cpus-per-task=30

module load StdEnv/2020
module load plink/1.9b_6.21-x86_64
module load vcftools/0.1.16
module load admixture/1.3.0

VCF=input.vcf.gz
OUTVCF=filtered.lowLD.50.10.0.2.maf0.05.vcf.gz
INPUTADM=lowLD.maf0.05

plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.2 --threads 30 --out lowLD.50.10.0.2
tr : $'\t' < lowLD.50.10.0.2.prune.in > lowLD.50.10.0.2.snps.txt
vcftools --gzvcf $VCF --positions lowLD.50.10.0.2.snps.txt --maf 0.05 --recode --stdout|bgzip -c > $OUTVCF

#run admixture
plink --vcf $OUTVCF --make-bed --double-id --out $INPUTADM --allow-extra-chr --set-missing-var-ids @:#
#the next two commands are needed to replace chromosome names when they have non-numeric names (often the case), which is not accepted by admixture
#to check the name of your chromosomes, you can use bcftools
bcftools query -f '%CHROM\n' $VCF | head
sed -i 's/Ha412HOChr//g' $INPUTADM.bim
awk '{$1=0;print $0}' $INPUTADM.bim > $INPUTADM.bim.tmp
for i in `seq 1 10`;
do  admixture --cv=10 $INPUTADM.bed $i -j32 > out.admixture.${i}.out;
do paste pops.file.txt lowLD.maf0.05.$i.Q > ancestry.K$i.txt;
done

#in R
#average the ancestry proportions by populations
ancestry.K2=read.table("ancestry.K2.txt", header=FALSE)
mean_pops_ancestry_K2=aggregate(V3 ~ V2, data=ancestry.k2, FUN=mean)

#run pca using snprelate
#you only need to change the name of the vcf with yours

library(SNPRelate)
vcf="filtered.lowLD.50.10.0.2.maf0.05.vcf.gz"
snpgdsVCF2GDS(vcf, out.fn = "yourfile.gds", method="copy.num.of.ref")
genofile=snpgdsOpen("yourfile.gds")
#the next two lines can be used to filter by ld and/or missing rate, in order to use the same set of snps in all analyses, I prefer to filter the vcfs before and use all the snps here
#snpset=snpgdsLDpruning(genofile, autosome.only=FALSE, ld.threshold=0.2, missing.rate=0.7)
#snpset.id=unlist(snpset)
#it allows multithread, though it is so fast that you won't need to use more than 8 threads for most datasets
pca=snpgdsPCA(genofile, autosome.only=FALSE,snp.id=NULL, num.thread=8)
snpgdsClose(my.genofile)

#plot pca
plot(pca)

#run abba-baba test to identify populations with evidence of gene flow with cultivated carrot
#I use Dsuite to perform abba-baba tests (https://github.com/millanek/Dsuite), you need three files for this analysis: 1) a vcf with all the samples, including at least one outgroup sample; 2) a TREE.txt file indicating the trio of populations you want to evaluate: (Outgroup,(CROP,(NAT,INT)));and 3) a SETS.txt file with two columns: the first with the name of the samples in the vcf file and the second with the population/group information (e.g., Outgroup, CROP, NAT, or INV)
#IMPORTANT: instead of subsetting the vcfs, you can ignore samples using the xxx keyword

./Dsuite Dtrios -t TREE.txt -c -o out_pop.txt INPUT_FILE.vcf SETS.txt
