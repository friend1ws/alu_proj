#! /usr/bin/env bash

##########
# chain file
if [ ! -f hg38ToHg19.over.chain.gz ]
then
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
fi

if [ ! -f hg19ToHg38.over.chain.gz ]
then
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
fi
##########

# 1000 geonme SV file
if [ ! -f ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz ]
then
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz
fi

if [ ! -f ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.tbi ]
then
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz.tbi
fi

bcftools filter -i 'INFO/SVLEN > 250 && INFO/SVTYPE == "ALU"' ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz | cut -f 1-8 > 1000genomes.alu.hg19.vcf

python proc_1000genomes.py 1000genomes.alu.hg19.vcf | sort -k1,1 -k2,2n -k3,3n > 1000genomes.alu.hg19.bed

liftOver -bedPlus=3 1000genomes.alu.hg19.bed hg19ToHg38.over.chain.gz 1000genomes.alu.hg38.bed.tmp 1000genomes.alu.unmapped

sort -k1,1 -k2,2n -k3,3n 1000genomes.alu.hg38.bed.tmp > 1000genomes.alu.hg38.bed
##########

bgzip -c 1000genomes.alu.hg19.bed > 1000genomes.alu.hg19.bed.gz
tabix -p bed 1000genomes.alu.hg19.bed.gz

bgzip -c 1000genomes.alu.hg38.bed > 1000genomes.alu.hg38.bed.gz
tabix -p bed 1000genomes.alu.hg38.bed.gz

rm -rf 1000genomes.alu.hg19.vcf
rm -rf 1000genomes.alu.hg19.bed
rm -rf 1000genomes.alu.hg38.bed
rm -rf 1000genomes.alu.hg38.bed.tmp 
rm -rf 1000genomes.alu.unmapped

