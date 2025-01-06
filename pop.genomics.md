#estimate genetic diversity and divergence across the genome using pixy
pixy --stats pi,dxy,fst --vcf vcf.with.invariant.sites.vcf.gz --populations pops.file.txt --window_size 150000 --n_cores 16 --output_folder out_pixy/ --output_prefix out_pixy_150kb  

#Basic filtering for population genomics include minimum allele frequency and low linkage disequilibrium
#I use plink for filtering, the script for compute canada but it is similar to other servers

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
#bcftools query -f '%CHROM\n' $VCF | head
sed -i 's/Ha412HOChr//g' $INPUTADM.bim
awk '{$1=0;print $0}' $INPUTADM.bim > $INPUTADM.bim.tmp
for i in `seq 1 10`;do  admixture --cv=10 $INPUTADM.bed $i -j32 > out.admixture.${i}.out;done

#in R
#average the proportion of admixture for individual populations

#run pca using snprelate
#you only need to change the name of the vcf with yours
library(SNPRelate)
vcf="yourvcf.vcf.gz"
snpgdsVCF2GDS(vcf, out.fn = "yourfile.gds", method="copy.num.of.ref")
genofile=snpgdsOpen("yourfile.gds")
#the next two lines can be used to filter by ld and/or missing rate, in order to use the same set of snps in all analyses, I prefer to filter the vcfs before and use all the snps here
#snpset=snpgdsLDpruning(genofile, autosome.only=FALSE, ld.threshold=0.2, missing.rate=0.7)
#snpset.id=unlist(snpset)
#it allows multithread, though it is so fast that you won't need to use more than 8 threads for most datasets
pca=snpgdsPCA(genofile, autosome.only=FALSE,snp.id=NULL, num.thread=8)
snpgdsClose(my.genofile)

#plot pca differentiating wild and cultivated samples

#run abba-baba test to identify populations with evidence of gene flow with cultivated carrot

#back in r, run gradient forest models to identify the environmental variables that better explain introgression frequencies
