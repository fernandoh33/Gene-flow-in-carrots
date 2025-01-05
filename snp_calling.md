#!/bin/bash
#SBATCH --account=your_account
#SBATCH --time=4-00:00
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --cpus-per-task=48

#set variables, I previously created three directories: trimmed_reads/ bam_files/ out_fastqc/ and the demultiplexed samples are in a directory named fastq.samples/PE.samples
FASTQ=fastq.samples/PE.samples
TRIMDIR=trimmed_reads
GENOME=reference/DCARv3.definitions.fna
BAM=bam_files
Fw=.1.fq.gz
Rv=.2.fq.gz
VCF=raw.vcf.minQ30.vcf.gz

module load StdEnv/2023
module load gcc/12.3
module load fastp/0.23.4
module load samtools/1.20
module load fastqc/0.12.1
module load bamtools/2.5.2
module load bwa/0.7.18
module load bcftools/1.19
module load stacks/2.67
module load vcftools/0.1.16

#demultiplex the libraries using process_radtags, libraries are in a folder named libraries/ and demultiplex samples are output to fastq.samples or fastq.samples/PE.samples/
#note all libraries but one were single-end sequenced, for SE and PE libraries I mapped the reads separately
#the barcode file include the sample name and the unique barcode for that sample, with no header, here is an example:
#        TGCT	JOBRU00005_001
#        CAGCT	JOBRU00005_004
#        GACTC	JOBRU00005_005
#        ACGTC	JOBRU00005_006
#        CACGTT	JOBRU00005_009
#        GTTAGC	JOBRU00005_010
#        AGCATT	JOBRU00005_015
#        CTCCGA	JOBRU00005_018
#        TTGGCA	JOBRU00005_020
#        GATGTC	JOBRU00007_001

#for demultiplexing, you need to provide the restriction enzyme used, if ddRADseq was used, use --renz-1 enzyme1 --renz-2 enzyme2 instead of -e. Here is the manual of process_radtags with detailed information and the name of supported restriction enzymes (https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php)
#SE libraries
process_radtags -f libraries/s_CCH17ANXX_s_7_fastq.txt.gz -b barcodes_CCH17ANXX_7.fixed -o fastq.samples/ -c -q -r --threads 48 --inline-null -e apeKI 
#PE libraries
process_radtags -1 libraries/Brunet-Simon-6-30-2023-Daucus_S2_L007_R1_001.fastq.gz -2 libraries/Brunet-Simon-6-30-2023-Daucus_S2_L007_R2_001.fastq.gz -b barcodes_Daucus_S2_L007.fixed -o fastq.samples/ -c -q -r --threads 48 --inline-null -e apeKI

#Create Sequence List
ls $FASTQ | grep $Fw | sed 's/.1.fq.gz\|.2.fq.gz//' | sort -u > samples_list.txt

#Trim Sequence Reads, it's RADseq data where pcr is used to amplify the target fragments, then you expect many duplicates, so I do not evaluate duplication (using the command --dont_eval_duplication)
for i in $(cat samples_list.txt);
do fastp -i $FASTQ/$i$Fw -I $FASTQ/$i$Rv -o $TRIMDIR/trimmed_$i$Fw -O $TRIMDIR/trimmed_$i$Rv --cut_right --dont_eval_duplication -h $TRIMDIR/trimmed_$i'.html' -g -w 16;
done

java -XX:InitialHeapSize=256m -XX:MaxHeapSize=4g -version

#Run fastqc on clean reads
fastqc -o out_fastqc/ $TRIMDIR/*

#Run Alignment using bwa and samtools
for FILE in $(cat samples_list.txt);
do bwa mem -t 48 $GENOME $TRIMDIR/trimmed_$FILE$Fw $TRIMDIR/trimmed_$FILE$Rv|samtools view -bh|samtools sort -T tmp -@ 48 -o $BAM/$FILE'.bam';
samtools index -@ 48 $BAM/$FILE'.bam';
bamtools stats -in $BAM/$FILE'.bam' > $BAM/$FILE'_bamstats.txt';
done

#snp calliing per se
bcftools mpileup --threads 48 -a AD,DP,SP -Ou -f $GENOME -Q 30 -q 20 $BAM/*.bam | bcftools call --threads 48 -f GQ,GP -m -Oz -o $VCF

#filter raw vcf to produce a vcf with only snps (used in most of the analyses) and a vcf including invariant sites, used in some analyses (e.g., in pixy to estimate pi or dxy)
vcftools --gzvcf $VCF --min-meanDP 4 --mac 2 --max-missing 0.7 --min-alleles 2 --max-alleles 2 --remove-indels --recode --stdout|bgzip -c > carrot.v3ref.only.snps.vcf.gz
vcftools --gzvcf $VCF --min-meanDP 4 --mac 2 --max-missing 0.7 --max-alleles 2 --remove-indels --recode --stdout|bgzip -c > carrot.v3ref.with.invariant.sites.vcf.gz

