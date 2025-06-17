#!/usr/bin/env bash
#SBATCH --job-name=hisat2
#SBATCH --mail-user=lingsn@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=hisat2_%u_%A_%a.out
#SBATCH --array=0-2
#SBATCH --account=cweidman0
#SBATCH --time=2:00:00
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
"run_hisat2.py --input 7KO-Cyto-BS_processed_fastqs/KEH-Rep1-7KO-HEK293T-Cyto-BS --index /home/lingsn/Software/genome_indices/hisat2_hg38/hg38p14_tran/GCF_000001405.40_GRCh38.p14_tran --annotation /home/lingsn/Software/genome_indices/hisat2_hg38/hg38p14_tran/GCF_000001405.40_GRCh38.p14_genomic.gtf --output alignments -I -C 2 -S /home/lingsn/Software/genome_indices/hisat2_hg38/hg38p14_tran/GCF_000001405.40_GRCh38.p14_splice_sites.txt"

"run_hisat2.py --input 7KO-Cyto-BS_processed_fastqs/KEH-Rep2-7KO-HEK293T-Cyto-BS --index /home/lingsn/Software/genome_indices/hisat2_hg38/hg38p14_tran/GCF_000001405.40_GRCh38.p14_tran --annotation /home/lingsn/Software/genome_indices/hisat2_hg38/hg38p14_tran/GCF_000001405.40_GRCh38.p14_genomic.gtf --output alignments -I -C 2 -S /home/lingsn/Software/genome_indices/hisat2_hg38/hg38p14_tran/GCF_000001405.40_GRCh38.p14_splice_sites.txt"

"run_hisat2.py --input 7KO-Cyto-BS_processed_fastqs/KEH-Rep3-7KO-HEK293T-Cyto-BS --index /home/lingsn/Software/genome_indices/hisat2_hg38/hg38p14_tran/GCF_000001405.40_GRCh38.p14_tran --annotation /home/lingsn/Software/genome_indices/hisat2_hg38/hg38p14_tran/GCF_000001405.40_GRCh38.p14_genomic.gtf --output alignments -I -C 2 -S /home/lingsn/Software/genome_indices/hisat2_hg38/hg38p14_tran/GCF_000001405.40_GRCh38.p14_splice_sites.txt"
)
eval ${tasks[$SLURM_ARRAY_TASK_ID]}



