#!/usr/bin/env bash
#SBATCH --job-name=create_env
#SBATCH --mail-user=YOUR_UNIQNAME@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/dev/null
#SBATCH --account=cweidman0
#SBATCH --time=5:00
#SBATCH --mem=5000m
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
# To call this script:
# sbatch create_env.sbatch
################################################################################

module purge
eval "$(conda shell.bash hook)"
conda create -n RNA-STAR python=3.9 samtools=1.17 star=2.7.11b pysam=0.21.0 parasail-python=1.3.1 numpy=1.24 pandas=1.4 openpyxl=3.1 scipy=1.11 dask=2024.8.0
