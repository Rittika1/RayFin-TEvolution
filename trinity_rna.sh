#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --job-name=trinity
#SBATCH --time=600:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=100GB
#SBATCH --error=/scratch/rmallik1/PhD_EVILab//GenomeAnnouncement/tirnity_gar3.out


module load trinity/2.14.0

cd /scratch/rmallik1/PhD_EVILab/GenomeAnnouncement/RNAs-seq

Trinity --seqType fq --left gar-rnaseq/NS035_1.fq --right gar-rnaseq/NS035_2.fq --CPU 12 --max_memory 100G