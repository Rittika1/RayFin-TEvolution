#!/bin/bash
#SBATCH --partition=Nebula
#SBATCH --job-name=buscorunbichir
#SBATCH --time=46:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=100GB
#SBATCH --error=/scratch/rmallik1/PhD_EVILab/GenomeAnnouncement/longnsoegar_RNA_actinop.out

module load busco


filepath=$1
inputfile=$2
lineage=$3 ##-- can be vertebrata_odb10, or actinopterygii_odb10
outputfolder=$4
mode=$5 ##-- mode can be genome, proteins, transcriptome

# cd /scratch/rmallik1/PhD_EVILab/GenomeAnnouncement/longnosegar
echo "cd $1"
cd $1

# busco -i dataset/GCF_016835505.1_ASM1683550v1_genomic.fna -l bichir/vertebrata_odb10 -o busco_vertebrate_bichir --out_path /scratch/rmallik1/SIRPs/bichir -m genome >> busco_vertebrate_bichir.txt -f -c 32
# busco -i lepisosteus_osseus_noadapter_nodups_noemptylines.fasta -l vertebrata_odb10 -o busco_longnosegar -m genome -f 

echo "busco -i $2 -l $3 -o $4 -m $5 -f -c 32"
busco -i $2 -l $3 -o $4 -m $5 -f -c 32

###---How to run:
# ./../scripts_for_genome_announcement/buscorun.sh /scratch/rmallik1/PhD_EVILab/GenomeAnnouncement/PolypterusBichir/ polypterus_bichir_lapradei_nodups_noadapters_noemptylines.fsa actinopterygii_odb10 busco_bichir_actinopterygii
# ./../scripts_for_genome_announcement/buscorun.sh /scratch/rmallik1/PhD_EVILab/GenomeAnnouncement/PolypterusBichir/ RNA actinopterygii_odb10 busco_longnosegar_rna


##---how to run
##--Rerun script for gar RNA
## ./scripts_for_genome_announcement/buscorun.sh /scratch/rmallik1/PhD_EVILab/GenomeAnnouncement/longnosegar trinity_out_dir.Trinity.fasta actinopterygii_odb10 busco_rna_gar_rerun transcriptome

##-- how to run for bichir senegalus
