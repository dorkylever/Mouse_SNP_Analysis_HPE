# load python 3.9 and install CrossMap
module load python3/3.9.2

pip3 install CrossMap==0.7.0

# download chain and reference
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
gunzip mm39.fa.gz

CrossMap vcf --compress mm10ToMm39.over.chain.gz C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf.gz mm39.fa C57BL_6NJ.mgp.v5.snps.dbSNP142.GRCm39.vcf.gz
