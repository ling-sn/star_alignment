#!/usr/bin/env bash
#SBATCH --job-name=run_star
#SBATCH --mail-user=lingsn@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=run_star_%u_%A_%a.out
#SBATCH --array=0-11%2
#SBATCH --account=cweidman0
#SBATCH --time=4:00:00
#SBATCH --mem=40000m
#SBATCH --cpus-per-task=8
#SBATCH --partition=standard
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
################################################################################
# Edit the strings under 'declare -a tasks=(' to match your experiments.
#
# The #SBATCH --array variable above creates an array [0,1,2,3]. Change it so that the length
# is how many jobs you need (same as number of strings under $tasks).
#
# This script is submitted that many times, but only one line from $tasks is
# evaluated each time.
#
# For more info on #SBATCH variables, see https://arc.umich.edu/greatlakes/slurm-user-guide/
# and https://slurm.schedmd.com/sbatch.html
#
# This requires a conda environment with samtools and STAR (RNA-STAR), run create_env.sbatch 
# before running this file 
# 
# To call this script:
# sbatch run_star.sbatch
################################################################################

module purge
eval "$(conda shell.bash hook)"
conda activate ~/miniconda3/envs/RNA-STAR

declare -a tasks=(
"python3 run_star.py --input 7KO-Cyto-BS_processed_fastqs --genomeDir ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38 --runThreadN=8"
"python3 run_star.py --input 7KO-Cyto-NBS_processed_fastqs --genomeDir ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38 --runThreadN=8"
"python3 run_star.py --input 7KO-Nuc-BS_processed_fastqs --genomeDir ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38 --runThreadN=8"
"python3 run_star.py --input 7KO-Nuc-NBS_processed_fastqs --genomeDir ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38 --runThreadN=8"
"python3 run_star.py --input 7LKO-Cyto-BS_processed_fastqs --genomeDir ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38 --runThreadN=8"
"python3 run_star.py --input 7LKO-Cyto-NBS_processed_fastqs --genomeDir ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38 --runThreadN=8"
"python3 run_star.py --input 7LKO-Nuc-BS_processed_fastqs --genomeDir ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38 --runThreadN=8"
"python3 run_star.py --input 7LKO-Nuc-NBS_processed_fastqs --genomeDir ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38 --runThreadN=8"
"python3 run_star.py --input WT-Cyto-BS_processed_fastqs --genomeDir ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38 --runThreadN=8"
"python3 run_star.py --input WT-Cyto-NBS_processed_fastqs --genomeDir ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38 --runThreadN=8"
"python3 run_star.py --input WT-Nuc-BS_processed_fastqs --genomeDir ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38 --runThreadN=8"
"python3 run_star.py --input WT-Nuc-NBS_processed_fastqs --genomeDir ~~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38 --runThreadN=8"
)
eval ${tasks[$SLURM_ARRAY_TASK_ID]}
