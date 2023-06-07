#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --mail-user=lucascarter2025@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -t 4:00:00
#SBATCH --job-name=get_contacts

# Set your working directory
cd /projects/b1042/BackmanLab/HiC2/opt/juicer/work/Lamin_HiC/juicer_analysis/Rep1/aligned

# Get contacts from Hi-C maps by chromosome:
chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
for i in "${chromosomes[@]}"
do
    java -jar /projects/b1042/BackmanLab/HiC2/opt/juicer/scripts/common/juicer_tools.jar dump observed VC inter_30.hic $i $i BP 25000 /projects/b1042/BackmanLab/HiC2/opt/juicer/work/ActD_HiC/contact_data/CTRL/mega/contacts/chr${i}-observed_25Kb.txt
done
