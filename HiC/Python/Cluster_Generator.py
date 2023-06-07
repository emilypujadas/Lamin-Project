

## Load Libraries
import os
import bioframe

##-----------------------------Quest Allocation------------------------------##

allocationID = 'b1042'
partitionName = 'genomics'
email = 'lucascarter2025@u.northwestern.edu' #'Your_Email@u.northwestern.edu'

##----------------------------Directory Variables----------------------------##

## HiC files stored in /projects/b1042/BackmanLab/Emily/HiC

# Local machine dir here the scripts and directory structer is initially generated
output_path = '/Users/lucascarter/Documents/IBiS/Backman_Lab/projects/Emily/Lamins_HiC/HiC_Analysis_Dir/' # '/Users/Your_Name/Documents/Your_Working_Directory/'

# Location of juicer directory on cluster
juicer_dir = '/projects/b1042/BackmanLab/HiC2/opt/juicer' # '~/HiC/opt/juicer'

# Location of non-Juicer post-Hi-C analysis python/R code on cluster
analysis_code_dir = '/projects/b1042/BackmanLab/HiC2/opt/code_files/' # '~/HiC/code_files/'

# Location of working juicer directory on cluster where FASTQs will be processed
juicer_work_dir = '/projects/b1042/BackmanLab/HiC2/opt/juicer/work/' # '~/HiC/opt/juicer/work/'

##----------------------------Experiment Variables----------------------------##

# Directory name of the experiment you're analyzing
experiment_name = 'Lamin_HiC' # 'Name_Of_Experiment'

# A list of condition directories for each condition to be analyzed
conds = ['untreated','24hrAuxin', 'withdraw'] # ['Cond_1','Cond_2', ect..]

# A list of the number of replicates in each condition
reps_per_cond = [6,3] # [1,1]

# Name of genome to align to
genome = "hg19" # 'genome_name'

# Name of genome fa (if needed)
fasta_name = "GRCh19.primary_assembly.genome.fa" # 'name_of_fasta'

# Chromosomes to analyze contacts (depends on genome)
chroms = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X'}

# Name of txt file with restriction enzyme cuts of genome
restriction_sites =  'hg19_DpnII.txt' # "My_Restriction_Site_File.txt"

# Name of restriction enzyme
restriction_enzyme = 'DpnII' #'restriction_enzyme'

# Either "Rep" (analyze by replicate) or "mega" (pool all replicates together, analyze by condition)
keyword = "Rep" # 'Rep_or_Mega'


## Here, we create the directory structure for the HiC analysis ------------------------------------------##

##### Input: output file path, experiment name, conditions analyzed, replicates per condition

##### Output: Create correct directory/subdirectory structure

##--------------------------------------------------------------------------------------------------------##

def create_dir_struct(output_path,experiment_name,conds,reps_per_cond):

    experiment_dir = output_path+experiment_name
    if not os.path.isdir(experiment_dir):
        os.mkdir(experiment_dir)

    ## Run through each subdirectory and make sure it is there. If it isn't, mkdir
    subdirnames = ['juicer_analysis','contact_data','contact_domains','compartment_analysis']

    for subdirname in subdirnames:
        subdir = experiment_dir+'/'+subdirname+'/'

        if not os.path.isdir(subdir):
            os.mkdir(subdir)

## Here we generate the Bash scripts to process the fastq files using Juicer -----------------------------##

##### Input: output file path, Juicer software and working directories, experiment name, conditions analyzed, replicates per condition, aligned genome, file containing restriction sites for specified genome

##### Output: Make cluster files to perform juicer analysis on fastq files for each replicate in each condition

##--------------------------------------------------------------------------------------------------------##

def juicer_analysis_files(output_path,juicer_dir,juicer_work_dir,experiment_name,conds,reps_per_cond,genome,restriction_sites):
    juicer_subdir = "/juicer_analysis/"

    ## As before, check that directories exist. If not, make directories
    for i,cond in enumerate(conds):
        cond_dir = output_path+experiment_name+juicer_subdir+cond
        if not os.path.isdir(cond_dir):
            os.mkdir(cond_dir)

        for rep in range(1,reps_per_cond[i]+1):
            rep_dir = cond_dir+'/Rep'+str(rep)
            if not os.path.isdir(rep_dir):
                os.mkdir(rep_dir)

            rep_dir_cluster = juicer_work_dir+experiment_name+juicer_subdir+cond+"/Rep"+str(rep)

            if not os.path.isdir(rep_dir+'/fastq'):
                os.mkdir(rep_dir+'/fastq')

            ## Open and write file for each rep called "juicer_analysis.sh"
            f = open(rep_dir+'/juicer_analysis.sh','w')

            f.write(

"#!/bin/bash\n\
#SBATCH -A "+allocationID+"\n\
#SBATCH -p "+partitionName+"\n\
#SBATCH -N 1\n\
#SBATCH --mem=100GB\n\
#SBATCH --ntasks-per-node=16\n\
#SBATCH --mail-user="+email+"\n\
#SBATCH --mail-type=BEGIN\n\
#SBATCH --mail-type=END\n\
#SBATCH -t 24:00:00\n\
#SBATCH --job-name="+experiment_name+"_"+cond+"_Rep"+str(rep)+"_juicer_analysis\n\n\
# Load necessary modules\n\
module load bwa/0.7.17\n\n\
# Set your working directory\n\
cd "+rep_dir_cluster+"\n\n\
# Perform Juicer Analysis on fastq files:\n"
+juicer_dir+"/scripts/juicer.sh -D "+juicer_dir+" -g "+genome+" -y "+juicer_dir+"/restriction_sites/"+restriction_sites+" -p "+juicer_dir+"/chrom_sizes/"+genome+".chrom.sizes\n\n\
# Perform cleanup on analysis files:\n"
+juicer_dir+"/scripts/common/cleanup.sh\n\n\
# Unzip merged_nodups file to perform mega analysis:\n\
gunzip aligned/merged_nodups.txt.gz")

            f.close()

## This function merges replicates processed in into a single matrix by condition -------------------------##

##### Input: output file path, Juicer software and working directories, experiment name, conditions analyzed, replicates per condition, aligned genome, restriction enzyme

##### Output: Make cluster files to merge Juicer analysis for replicates to create mega file for each condition

##--------------------------------------------------------------------------------------------------------##

def mega_files(output_path,juicer_dir,juicer_work_dir,experiment_name,conds,reps_per_cond,genome,restriction_enzyme):
    juicer_subdir = "/juicer_analysis/"

    ## As before, check that directories exist. If not, make directories
    for i,cond in enumerate(conds):
        if (reps_per_cond[i]>1):
            cond_dir = output_path+experiment_name+juicer_subdir+cond
            if not os.path.isdir(cond_dir):
                os.mkdir(cond_dir)

            cond_dir_cluster = juicer_work_dir+experiment_name+juicer_subdir+cond

            for rep in range(1,reps_per_cond[i]+1):
                rep_dir = cond_dir+'/Rep'+str(rep)
                if not os.path.isdir(rep_dir):
                    os.mkdir(rep_dir)

                ## Open and write file for each rep called "mega.sh"
                f = open(cond_dir+'/mega.sh','w')

                f.write(
"#!/bin/bash\n\
#SBATCH -A "+allocationID+"\n\
#SBATCH -p "+partitionName+"\n\
#SBATCH -N 1\n\
#SBATCH --ntasks-per-node=16\n\
#SBATCH --mail-user="+email+"\n\
#SBATCH --mail-type=BEGIN\n\
#SBATCH --mail-type=END\n\
#SBATCH -t 48:00:00\n\
#SBATCH --job-name="+experiment_name+"_"+cond+"_mega\n\n\
# Set your working directory\n\
cd "+cond_dir_cluster+"\n\n\
# Merge hi-c analysis by replicate into mega analysis by condition:\n"
+juicer_dir+"/scripts/common/mega.sh -D "+juicer_dir+" -g "+genome+" -b "+restriction_enzyme+"\n\
gzip "+cond_dir_cluster+'/mega/aligned/merged')

                f.close()

## function to generate bash scripts that calculate map resolution for replicates and merged files -------##

##### Input: output file path, Juicer software and Juicer working directory, experiment name, conditions analyzed, replicates per condition

##### Output: Make cluster files to perform 3D-DNA analysis on merged_nodups.txt files from Juicer analysis for each replicate in each condition

##--------------------------------------------------------------------------------------------------------##

def calculate_map_resolution(output_path,juicer_dir,juicer_work_dir,experiment_name,conds,reps_per_cond):
    juicer_subdir = "/juicer_analysis/"

    ## As before, check that directories exist. If not, make directories
    for i,cond in enumerate(conds):
        cond_dir = output_path+experiment_name+juicer_subdir+cond
        if not os.path.isdir(cond_dir):
            os.mkdir(cond_dir)

        if (keyword =='Rep'):

            # Write get contacts file for each replicate
            for rep in range(1,reps_per_cond[i]+1):
                rep_dir = cond_dir+'/Rep'+str(rep)
                if not os.path.isdir(rep_dir):
                    os.mkdir(rep_dir)

                juicer_rep_dir = juicer_work_dir+experiment_name+juicer_subdir+cond+"/Rep"+str(rep)

                ## Open and write file for each rep called "calculate_map_resolution.sh"
                f = open(rep_dir+'/calculate_map_resolution.sh','w')

                ## This loop writes generates Bash scripts to calculate map resolution for replicates
                f.write(
"#!/bin/bash\n\
#SBATCH -A "+allocationID+"\n\
#SBATCH -p "+partitionName+"\n\
#SBATCH -N 1\n\
#SBATCH --ntasks-per-node=16\n\
#SBATCH --mail-user="+email+"\n\
#SBATCH --mail-type=BEGIN\n\
#SBATCH --mail-type=END\n\
#SBATCH -t 4:00:00\n\
#SBATCH --job-name="+experiment_name+"_"+cond+"_Rep"+str(rep)+"_calculate_map_resolution\n\n\
# Set your working directory\n\
cd "+juicer_rep_dir+"/aligned\n\n\
# Calculate map resolution:\n"
+juicer_dir+"/misc/calculate_map_resolution.sh merged_nodups.txt 50bp_coverage.txt\n")

                f.close()
## End script generator for replicate map calculate_map_resolution

        elif (keyword == "mega"):
            # Write get contacts file for mega file
            mega_dir = cond_dir+"/mega"
            if not os.path.isdir(mega_dir):
                os.mkdir(mega_dir)

            juicer_mega_dir = juicer_work_dir+experiment_name+juicer_subdir+cond+"/mega"

            f = open(mega_dir+'/calculate_map_resolution_mega.sh','w')

## This loop writes generates Bash scripts to calculate map resolution for merged files

            f.write(
"#!/bin/bash\n\
#SBATCH -A "+allocationID+"\n\
#SBATCH -p "+partitionName+"\n\
#SBATCH -N 1\n\
#SBATCH --ntasks-per-node=16\n\
#SBATCH --mail-user="+email+"\n\
#SBATCH --mail-type=BEGIN\n\
#SBATCH --mail-type=END\n\
#SBATCH -t 4:00:00\n\
#SBATCH --job-name="+experiment_name+"_"+cond+"_mega_calculate_map_resolution\n\n\
# Set your working directory\n\
cd "+juicer_mega_dir+"/aligned\n\n\
# Calculate map resolution:\n"
+juicer_dir+"/misc/calculate_map_resolution.sh merged_nodups.txt 50bp_coverage.txt\n")

            f.close()
## End script generator for mega(merged) map calculate_map_resolution

        else:
            print("Enter correct keyword")

## function to generate bash scripts that compress replicates --------------------------------------------##

##### Input: output file path, Juicer working directory, experiment name, conditions analyzed, replicates per condition

##### Output: Make cluster files to merge Juicer analysis for replicates to create mega file for each condition

##--------------------------------------------------------------------------------------------------------##

def compress_files(output_path,juicer_work_dir,experiment_name,conds,reps_per_cond):
    juicer_subdir = "/juicer_analysis/"

    ## As before, check that directories exist. If not, make directories
    for i,cond in enumerate(conds):
        cond_dir = output_path+experiment_name+juicer_subdir+cond
        if not os.path.isdir(cond_dir):
                os.mkdir(cond_dir)

        cond_dir_cluster = juicer_work_dir+experiment_name+juicer_subdir+cond

        for rep in range(1,reps_per_cond[i]+1):
            rep_dir = cond_dir+'/Rep'+str(rep)
            if not os.path.isdir(rep_dir):
                os.mkdir(rep_dir)

            f = open(cond_dir+'/compress_files.sh','w')

            ## This loop writes generates Bash scripts to compress files
            f.write(
"#!/bin/bash\n\
#SBATCH -A "+allocationID+"\n\
#SBATCH -p "+partitionName+"\n\
#SBATCH -N 1\n\
#SBATCH --ntasks-per-node=16\n\
#SBATCH --mail-user="+email+"\n\
#SBATCH --mail-type=BEGIN\n\
#SBATCH --mail-type=END\n\
#SBATCH -t 48:00:00\n\
#SBATCH --job-name="+experiment_name+"_"+cond+"_compress_files\n\n\
# Set your working directory\n\
cd "+cond_dir_cluster+"\n\n\
# Compress all stats files that were unzipped for mega analysis:\n\
gzip **/aligned/*.txt\n")

            f.close()

##------------- Call main method functions ----------------------------------------------------------------------##

## Generate the directory structure where each of the initial analysis scripts will be stored
create_dir_struct(output_path,experiment_name,conds,reps_per_cond)

## Generate Juicer analysis files that process FASTQs into initial .HiC files
juicer_analysis_files(output_path,juicer_dir,juicer_work_dir,experiment_name,conds,reps_per_cond,genome,restriction_sites)

## Generate scripts to merge replicates in single .HiC files to increase resolution
mega_files(output_path,juicer_dir,juicer_work_dir,experiment_name,conds,reps_per_cond,genome,restriction_enzyme)

## Generate scripts to calculate resolution of HiC data
calculate_map_resolution(output_path,juicer_dir,juicer_work_dir,experiment_name,conds,reps_per_cond)

## Generate scripts to compress files for ease of storage and use
compress_files(output_path,juicer_work_dir,experiment_name,conds,reps_per_cond)

##---------------------------------------------------------------------------------------------------------------##


##------------- Define Secondary Functions ----------------------------------------------------------------------##

## function to generate bash scripts for calling domains with arrowhead -----------------------------------------##

##### Input: output file path, Juicer software and working directories, experiment name, conditions analyzed, replicates per condition, keyword (either "Rep" or "mega"), matrix resolution and normalization method

##### Output: Make cluster files to obtain contact information from Hi-C maps generated from Juicer Analysis

##---------------------------------------------------------------------------------------------------------------##

def get_contact_domains_arrowhead(output_path,juicer_dir,juicer_work_dir,experiment_name,conds,reps_per_cond,keyword,res,matrix_norm):
    juicer_subdir = "/juicer_analysis/"
    cont_domain_subdir = "/contact_domains/arrowhead_domains/"
    if not os.path.isdir(output_path+experiment_name+cont_domain_subdir):
            os.mkdir(output_path+experiment_name+cont_domain_subdir)

    for i,cond in enumerate(conds):
        cond_dir = output_path+experiment_name+cont_domain_subdir+cond

        if not os.path.isdir(cond_dir):
            os.mkdir(cond_dir)

        if (keyword =='Rep'):
            # Write get contacts file for each replicate
            for rep in range(1,reps_per_cond[i]+1):
                rep_dir = cond_dir+'/Rep'+str(rep)
                if not os.path.isdir(rep_dir):
                    os.mkdir(rep_dir)

                juicer_rep_dir = juicer_work_dir+experiment_name+juicer_subdir+cond+"/Rep"+str(rep)
                cont_domain_rep_dir = juicer_work_dir+experiment_name+cont_domain_subdir+cond+"/Rep"+str(rep)

                f = open(rep_dir+'/contact_domains_'+matrix_norm+'_norm.sh','w')
                f.write("#!/bin/bash\n\
#SBATCH -A "+allocationID+"\n\
#SBATCH -p "+partitionName+"\n\
#SBATCH -N 1\n\
#SBATCH --ntasks-per-node=16\n\
#SBATCH --mail-user="+email+"\n\
#SBATCH --mail-type=BEGIN\n\
#SBATCH --mail-type=END\n\
#SBATCH -t 48:00:00\n\
#SBATCH --job-name="+experiment_name+"_"+cond+"_Rep"+str(rep)+"get_contact_domains\n\n\
# Set your working directory\n\
cd "+juicer_rep_dir+"/aligned\n\n\
# Arrowhead annotation of contact domains:\n"
+juicer_dir+"/scripts/common/juicer_tools arrowhead -r "+str(res*1000)+" -k "+matrix_norm+" inter_30.hic "+cont_domain_rep_dir+"/inter_30_contact_domains\n")
                f.close()

        elif (keyword == "mega"):
            # Get contact domains for mega file
            mega_dir = cond_dir+"/mega"
            if not os.path.isdir(mega_dir):
                os.mkdir(mega_dir)

            juicer_mega_dir = juicer_work_dir+experiment_name+juicer_subdir+cond+"/mega"
            cont_domain_mega_dir = juicer_work_dir+experiment_name+cont_domain_subdir+cond+"/mega"

            f = open(mega_dir+'/contact_domains_'+matrix_norm+'_norm_mega.sh','w')
            f.write("#!/bin/bash\n\
#SBATCH -A "+allocationID+"\n\
#SBATCH -p "+partitionName+"\n\
#SBATCH -N 1\n\
#SBATCH --ntasks-per-node=16\n\
#SBATCH --mail-user="+email+"\n\
#SBATCH --mail-type=BEGIN\n\
#SBATCH --mail-type=END\n\
#SBATCH -t 48:00:00\n\
#SBATCH --job-name="+experiment_name+"_"+cond+"_mega_get_contact_domains\n\n\
# Set your working directory\n\
cd "+juicer_mega_dir+"/aligned\n\n\
# Arrowhead annotation of contact domains:\n"
+juicer_dir+"/scripts/common/juicer_tools arrowhead -r "+str(res*1000)+" -k "+matrix_norm+" inter_30.hic "+cont_domain_mega_dir+"/inter_30_contact_domains\n")
            f.close()

        else:
            print("Enter correct keyword")

############### Experiment input #################################################################################
res = 20 #Matrix resolution in kbp
matrix_norm = 'KR' # matrix normalization method (OPTIONS: NONE/VC/VC_SQRT/KR)


##------------- Call second method functions --------------------------------------------------------------------##

get_contact_domains_arrowhead(output_path,juicer_dir,juicer_work_dir,experiment_name,conds,reps_per_cond,keyword,res,matrix_norm)

##---------------------------------------------------------------------------------------------------------------##
