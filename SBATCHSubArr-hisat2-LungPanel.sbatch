#!/usr/bin/env bash
#SBATCH --job-name=hisat2
#SBATCH --mail-user=cweidman@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=hisat2_%u_%A_%a.out
#SBATCH --array=0-19
#SBATCH --account=cweidman0
#SBATCH --time=8:00:00
#SBATCH --mem=16000m
#SBATCH --partition=standard
#SBATCH --ntasks-per-node=2
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
# This requires a conda environment for genome alignment, edit to your named version in the activate command
# 
# To call this script:
# sbatch SBATCH-arrXXX.sbatch
################################################################################

module purge
eval "$(conda shell.bash hook)"
conda activate ~/miniconda3/envs/RNA-SEQ

declare -a tasks=(
"run_hisat2.py --input Sample_6_processed_fastqs/A549-CDDO24hr-Rep1 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_16_processed_fastqs/A549-CDDO24hr-Rep2 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_5_processed_fastqs/A549-DMSO-Rep1 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_15_processed_fastqs/A549-DMSO-Rep2 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_4_processed_fastqs/H520-CDDO24hr-Rep1 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_14_processed_fastqs/H520-CDDO24hr-Rep2 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_3_processed_fastqs/H520-DMSO-Rep1 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_13_processed_fastqs/H520-DMSO-Rep2 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_2_processed_fastqs/HNY-CDDO24hr-Rep1 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_12_processed_fastqs/HNY-CDDO24hr-Rep3 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_1_processed_fastqs/HNY-DMSO-Rep1 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_11_processed_fastqs/HNY-DMSO-Rep3 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_10_processed_fastqs/HNY-KAV13-DMSO-RPMI-Blast-Rep2 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_20_processed_fastqs/HNY-KAV13-DMSO-RPMI-Blast-Rep3 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_9_processed_fastqs/HNY-KAV13-DMSO-RPMI-Rep2 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_19_processed_fastqs/HNY-KAV13-DMSO-RPMI-Rep3 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_8_processed_fastqs/SKMES1-CDDO24hr-Rep1 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_18_processed_fastqs/SKMES1-CDDO24hr-Rep2 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_7_processed_fastqs/SKMES1-DMSO-Rep1 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
"run_hisat2.py --input Sample_17_processed_fastqs/SKMES1-DMSO-Rep2 --index /home/cweidman/Software/genome_indices/hisat2_hg38/hg38_tran/genome_tran --annotation /home/cweidman/Software/genome_indices/homo_sapiens/Homo_sapiens.GRCh38.96.gtf --output alignments -I -C 2 -S /home/cweidman/Software/genome_indices/hisat2_hg38/Homo_sapiens.GRCh38.96.splicesites.txt"
)
eval ${tasks[$SLURM_ARRAY_TASK_ID]}



