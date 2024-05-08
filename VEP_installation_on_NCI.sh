module load singularity

singularity pull --name vep.sif docker://ensemblorg/ensembl-vep:release_111.0

singularity exec vep.sif INSTALL.pl -c mouse_data -a cf -s mus_musculus -y GRCm39
