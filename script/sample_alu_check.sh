#! /bin/bash

SAMPLE=$1

bcftools view -s ${SAMPLE} ../db/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz | bcftools filter -i 'INFO/SVLEN > 250 && INFO/SVTYPE == "ALU"' > ${SAMPLE}.alu.vcf

python proc_1000genomes.py ${SAMPLE}.alu.vcf > ${SAMPLE}.alu.hg19.bed

liftOver -bedPlus=3 ${SAMPLE}.alu.hg19.bed ../db/hg19ToHg38.over.chain.gz ${SAMPLE}.alu.hg38.bed.tmp ${SAMPLE}.alu.unmapped
 
sort -k1,1 -k2,2n -k3,3n ${SAMPLE}.alu.hg38.bed.tmp > ${SAMPLE}.alu.hg38.bed

rm -rf ${SAMPLE}.alu.vcf
rm -rf ${SAMPLE}.alu.hg19.bed
rm -rf ${SAMPLE}.alu.hg38.bed.tmp
rm -rf ${SAMPLE}.alu.unmapped

