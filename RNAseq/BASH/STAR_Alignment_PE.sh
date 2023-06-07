
#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=lucascarter2025@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --time 48:00:00
#SBATCH --job-name=STAR_align
#SBATCH --output=outlog
#SBATCH --error=errlog

# Load modules
module load python/anaconda3.6
module load STAR/2.6.0
module load rsem/1.3.3
module load samtools/1.6
module load deeptools/3.1.1
module load fastqc

# Set working directory
cd /projects/b1042/BackmanLab/Lucas/Adli_Lamin_RNA-seq/Adli_Lamin_RNA-seq_RawData

echo "$(date): Processing alignment of FASTQ sequence files"

echo "Generating directories."
# make sure the STAR and RSEM references are index references
genomefold=/projects/b1042/BackmanLab/ref_hg38ensembl/hg38.index/Gencode/STAR
gtffold=/projects/b1042/BackmanLab/ref_hg38ensembl/hg38.index/Gencode/gencode.v41.annotation.gtf
RSEMfold=/projects/b1042/BackmanLab/ref_hg38ensembl/hg38.index/Gencode/RSEM/RSEM

echo "STAR genome assembly is located in $genomefold"
echo "Annotations file is located in $gtffold"
echo "RSEM reference files are located in $RSEMfold"

outfold="$(dirname "$(pwd)")/${PWD##*/}_output"
[ ! -d $outfold ]&&mkdir $outfold

# Make QC directory
qcfold="$(dirname "$(pwd)")/${PWD##*/}_fastqc"
[ ! -d $qcfold ]&&mkdir $qcfold

# Make coverage directory
sortfold="$(dirname "$(pwd)")/${PWD##*/}_coverage"
[ ! -d $sortfold ]&&mkdir $sortfold

# Make .BAM directory
bamfold="$(dirname "$(pwd)")/${PWD##*/}_BAMs"
[ ! -d $bamfold ]&&mkdir $bamfold

# Make .OUT directory
outfold2="$(dirname "$(pwd)")/${PWD##*/}_OUT"
[ ! -d $outfold2 ]&&mkdir $outfold2

# Make .HTSEQ_COUNTS directory
cntfold="$(dirname "$(pwd)")/${PWD##*/}_HTSEQ_counts"
[ ! -d $cntfold ]&&mkdir $cntfold

# Make .RSEM_ COUNTS directory
cntfold2="$(dirname "$(pwd)")/${PWD##*/}_RSEM_counts"
[ ! -d $cntfold2 ]&&mkdir $cntfold2

# Define function align_data
function align_data {

R1=$1;
R2=$(echo $R1 | sed 's/R1/R2/g')

STAR --runThreadN 24 --quantMode TranscriptomeSAM --genomeDir $genomefold --readFilesIn ${R1} ${R2} --readFilesCommand zcat --outFileNamePrefix "$outfold/${R1%.*.*}_";
}

# Define function sam_tools
function sam_tools {
R1=$1;
BAM=$(echo $R1 | sed 's/_R1_001_Aligned.out.sam/_Aligned.out.bam/g')
BAMsort=$(echo $BAM | sed 's/_Aligned.out.bam/_Aligned.out.sort.bam/g')
BW=$(echo $BAMsort | sed 's/_Aligned.out.sort.bam/./g')

samtools view -b $R1 > $BAM
samtools sort $BAM > $BAMsort
samtools index $BAMsort

# Bam coverage to generate BigWig
bamCoverage --bam $BAMsort --normalizeUsing CPM --outFileName "${BW}-.bigWig" --filterRNAstrand reverse --binSize 1 --numberOfProcessors 60
bamCoverage --bam $BAMsort --normalizeUsing CPM --outFileName "${BW}+.bigWig" --filterRNAstrand forward --binSize 1 --numberOfProcessors 60

}

echo "Aligning files"
for file in *R1_001.fastq.gz;
do echo "${file}";
align_data $file;
done

echo "Running Fastqc on all FASTA files"
fastqc -t 12 *.fastq.gz
mv *.html $qcfold
mv *fastqc.zip $qcfold

# Change directory to $Outfold
cd $outfold

echo "Sorting alignments and generating BigWigs"
for file in *Aligned.out.sam;
do echo "${file}";
sam_tools $file;
done

echo "Gathering statistics on BAM files"
for file in *sort.bam;
do echo "$file";

samtools flagstat $file > \
${file}.FLAGSTAT.txt

done

echo "Generating HTSEQ count files from sample alignments"
for file in *.sam;
do echo "${file}";
output_file=$(echo $file | sed 's/_R1_001_Aligned.out.sam/.counts/g')
htseq-count -f sam -r pos -s no -t exon -i gene_id -m union "${file}" $gtffold > "${output_file}"
done

echo "Generating RSEM count files from sample alignments"
for file in *.toTranscriptome.out.bam;
do echo "${file}";
output_file=$(echo $file | sed 's/_R1_001_Aligned.toTranscriptome.out.bam/.rsem/g')
rsem-calculate-expression -p 24 --alignments --paired-end "${file}" $RSEMfold "${output_file}"
done

echo "Moving processed files to respective directories."

mv *FLAGSTAT.txt $qcfold
mv *.bam $bamfold
mv *.out $outfold2
mv *.counts $cntfold
mv *.results $cntfold2
mv *.bigWig $sortfold

echo "FASTQ processing complete!"
