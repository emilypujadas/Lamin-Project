#!/bin/sh
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH -N 1
#SBATCH -n 16 ## Number of node
#SBATCH --time=12:00:00
#SBATCH --mem=50G
#SBATCH --mail-user=lucascarter2025@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --job-name=TopDom

cd /projects/b1042/BackmanLab/HiC2/opt/juicer/work/Lamin_HiC/contact_data

module load R/4.1.0

Rscript --vanilla --verbose mergeContacts.R
