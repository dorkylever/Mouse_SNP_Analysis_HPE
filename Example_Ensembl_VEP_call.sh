#!/bin/bash
#PBS -P ri21
#PBS -q express
#PBS -l walltime=24:00:00
#PBS -l mem=190GB
##PBS -l jobfs=200GB
#PBS -lstorage=gdata/if89+gdata/nm24
#PBS -l ncpus=48
##PBS -l ngpus=4
##PBS -V
## For licensed software, you have to specify it to get the job running. For unlicensed software, you should also specify it to help us analyse the software usage on our system.
##PBS -l software=lama_get_walkthrough_data
## The job will be executed from current working directory instead of home.


#### Some preliminary stuff to get it working ####
# Get the current working directory of this script


#module load nci-parallel/1.0.0a
#export ncores_per_task=4
#export ncores_per_numanode=12

cd /g/data/nm24/kd0793/mouse_data
#pwd

module load singularity



singularity exec /g/data/nm24/kd0793/ables-software-installations/ensembl-vep/111.0/vep.sif \
  vep --dir /g/data/nm24/kd0793/ables-software-installations/ensembl-vep/111.0/mouse_data/ \
      --cache --offline --format vcf --vcf --force_overwrite --compress_output bgzip \
      --sift b --pick --assembly GRCm39 --species mus_musculus --verbose  --fork 48 \
      --regulatory --gene_phenotype --mirna --numbers --max_af --canonical --domains --biotype\
      --input_file C57BL_6NJ.mgp.v5.snps.dbSNP142.GRCm39.vcf.gz \
      --output_file pred_C57BL_6NJ.mgp.v5.snps.dbSNP142.GRCm39.vcf.gz \
      --plugin Conservation,/g/data/nm24/kd0793/mouse_data/gerp_conservation_scores.mus_musculus.GRCm39.bw --plugin NMD \
      --plugin Blosum62
