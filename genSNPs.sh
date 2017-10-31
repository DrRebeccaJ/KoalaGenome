
SPP_NAME="koala"
BASE_DIR="/home/jgb/koala/"

MAP_SPP_LIST=${BASE_DIR}/lists/LN02.libs

SPP_REF_FASTA=${BASE_DIR}/gen/phaCin_unsw_v4.1.fa
SPP_REF_DICT=${BASE_DIR}/gen/phaCin_unsw_v4.1.dict

SPP_MAP_DIR=${BASE_DIR}/gmap/
CLN_READ_DIR=/home/jgb/LN/

### Index reference file for bowtie aligner
 bowtie2-build ${SPP_REF_FASTA} ${SPP_REF_FASTA}
 samtools faidx ${SPP_REF_FASTA}
 SEQDICT=/home/jgb/software/picard-tools-1.88/CreateSequenceDictionary.jar
 java -jar $SEQDICT R=$SPP_REF_FASTA O=$SPP_REF_DICT

### Map to common set of references 
 cat ${MAP_SPP_LIST} | parallel --gnu --max-procs=8 perl /home/jgb/koala/scripts/pl/map2refs.pl -l {} -q ${CLN_READ_DIR} -m ${SPP_MAP_DIR} -r ${SPP_REF_FASTA}

### Do some GATK processing of the BAM files
ADDRG=/home/jgb/software/picard-tools-1.88/AddOrReplaceReadGroups.jar
GATK=/home/jgb/software/gatk-330/GenomeAnalysisTK.jar

cat ${MAP_SPP_LIST} | parallel --gnu --max-procs=8 perl /home/jgb/koala/scripts/pl/processBAM.pl -l {} -m ${SPP_MAP_DIR} -r ${SPP_REF_FASTA} -p $SEQDICT -g $GATK -a $ADDRG


### Merge BAMs and index merged bam
MERGESAM=/home/jgb/software/picard-tools-1.88/MergeSamFiles.jar
MERGEDBAMFILE=${SPP_MAP_DIR}samplesMerged.bam
LIBARRAY=$(cat ${MAP_SPP_LIST})
BAMLIST=""
for B in $LIBARRAY
do
   BAMLIST="$BAMLIST I=${SPP_MAP_DIR}${B}.proc.bam"
done
echo "$BAMLIST"

java -jar $MERGESAM $BAMLIST SO=coordinate AS=true VALIDATION_STRINGENCY=SILENT O=$MERGEDBAMFILE
samtools index $MERGEDBAMFILE


cd ${SPP_MAP_DIR}
### Run GATK
HAPCALL="HaplotypeCaller"
RAWVCF=samplesMerged.vcf

java -Xmx12g -jar ${GATK} -R ${SPP_REF_FASTA} -T ${HAPCALL} -I ${MERGEDBAMFILE} -o ${RAWVCF}

#vcftools --vcf samplesMerged.vcf --remove-indels --hwe 0.00001 --minDP 12 --minGQ 30 --min-alleles 2 --max-alleles 2 --recode --out samplesMerged.filt

#vcftools --vcf samplesMerged.filt.recode.vcf --012 --out samplesMerged.filt




