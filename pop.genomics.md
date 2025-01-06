#estimate genetic diversity and divergence across the genome using pixy

#Basic filtering for population genomics include minimum allele frequency and low linkage disequilibrium
#I use plink for this

#run admixture


#in R
#average the proportion of admixture for individual populations

#run pca using snprelate
#you only need to change the name of the vcf with yours
vcf="yourvcf.vcf.gz"
snpgdsVCF2GDS(vcf, out.fn = "yourfile.gds", method="copy.num.of.ref")
genofile=snpgdsOpen("yourfile.gds")
#the next two lines can be used to filter by ld and/or missing rate, I prefer to filter the vcfs before and use all the snps here
#snpset=snpgdsLDpruning(genofile, autosome.only=FALSE, ld.threshold=0.2, missing.rate=0.7)
#snpset.id=unlist(snpset)
#it allows multithread, though it is so fast that you won't need to use more than 8 threads for most datasets
pca=snpgdsPCA(genofile, autosome.only=FALSE,snp.id=NULL, num.thread=8)
#plot pca differentiating wild and cultivated samples

#run abba-baba test to identify populations with evidence of gene flow with cultivated carrot

#back in r, run gradient forest models to identify the environmental variables that better explain introgression frequencies
