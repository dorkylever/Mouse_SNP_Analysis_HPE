
module load bcftools/1.12

bcftools view -s C3H_HeH,C57BL_6NJ -O z -o mgp_REL2021_snps_C3H_BL6.vcf.gz mgp_REL2021_snps.vcf.gz
bcftools view -s C3H_HeH --private mgp_REL2021_snps_C3H_BL6.vcf.gz -Ou | bcftools view -e 'FORMAT/FI[*]="." || FORMAT/FI[*]!=1' -Oz > unique_to_C3H.vcf.gz 
