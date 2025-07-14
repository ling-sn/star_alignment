#!/usr/bin/env bash
#SBATCH --job-name=realignGap
#SBATCH --mail-user=lingsn@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=realignGap_%u_%A_%a.out
#SBATCH --array=0-11%2
#SBATCH --account=cweidman0
#SBATCH --time=3:00:00
#SBATCH --mem=40000m
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
# sbatch realignGap.sbatch
################################################################################

module purge
eval "$(conda shell.bash hook)"
conda activate ~/miniconda3/envs/RNA-STAR

declare -a tasks=(
"python3 realignGap.py --folder_name 7KO-Cyto-BS_processed_fastqs --fasta_dir ~/scratch/star/star_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa --discard"

"python3 realignGap.py --folder_name 7KO-Cyto-NBS_processed_fastqs --fasta_dir ~/scratch/star/star_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa --discard"

"python3 realignGap.py --folder_name 7KO-Nuc-BS_processed_fastqs --fasta_dir ~/scratch/star/star_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa --discard"

"python3 realignGap.py --folder_name 7KO-Nuc-NBS_processed_fastqs --fasta_dir ~/scratch/star/star_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa --discard"

"python3 realignGap.py --folder_name 7LKO-Cyto-BS_processed_fastqs --fasta_dir ~/scratch/star/star_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa --discard"

"python3 realignGap.py --folder_name 7LKO-Cyto-NBS_processed_fastqs --fasta_dir ~/scratch/star/star_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa --discard"

"python3 realignGap.py --folder_name 7LKO-Nuc-BS_processed_fastqs --fasta_dir ~/scratch/star/star_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa --discard"

"python3 realignGap.py --folder_name 7LKO-Nuc-NBS_processed_fastqs --fasta_dir ~/scratch/star/star_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa --discard"

"python3 realignGap.py --folder_name WT-Cyto-BS_processed_fastqs --fasta_dir ~/scratch/star/star_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa --discard"

"python3 realignGap.py --folder_name WT-Cyto-NBS_processed_fastqs --fasta_dir ~/scratch/star/star_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa --discard"

"python3 realignGap.py --folder_name WT-Nuc-BS_processed_fastqs --fasta_dir ~/scratch/star/star_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa --discard"

"python3 realignGap.py --folder_name WT-Nuc-NBS_processed_fastqs --fasta_dir ~/scratch/star/star_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa --discard"
)
eval ${tasks[$SLURM_ARRAY_TASK_ID]}