# **Running STAR alignment on processed fastp files**
## Necessary files
* `create_env.sbatch`
* STAR index
  * **Option 1 (Manual):** First, create the `star_hg38` folder according to the instructions below (see "_Creating the star_hg38 folder_"). Then, copy over `star_index.py` and `star_index.sbatch` and run the SBATCH file.
  * **Option 2 (Pre-Built):** Skip STAR index creation by directly using `--genomeDir /home/lingsn/scratch/star/star_hg38` in your `run_star.sbatch` file.
* `run_star.py` and `run_star.sbatch`
## Instructions
1. Run `create_env.sbatch` to create the RNA-STAR conda environment
3. Activate conda environment via `conda activate RNA-STAR`
4. Run `star_index.sbatch` one time to build hg38 genome index (_optional_)
5. Edit `run_star.sbatch` to match your experiments
6. Run `run_star.sbatch` to align each input folder using STAR
## Tools used in STAR alignment script
* **STAR** is used to align single-end reads (merged/unpaired) and paired reads (unmerged) according to the hg38 genome index
* **samtools** is used to merge single-end and paired .bam files together
## When do I use this pipeline?
This is used after running the fastp script on your raw data (fastq files). 
## Creating the star_hg38 folder
Read the following if you want to manually build the STAR index. Otherwise, feel free to skip this section.
1. Create an empty folder called `star_hg38` in the working directory containing your data folders
2. On Great Lakes Cluster, navigate to `/home/<uniqname>/umms-RNAlabDATA/Software/genome_indices/hisat2_hg38/hg38p14_tran/` and copy the following 4 files into `star_hg38`:
   * `GCF_000001405.40_GRCh38.p14_exons.txt`
   * `GCF_000001405.40_GRCh38.p14_genomic.fa`
   * `GCF_000001405.40_GRCh38.p14_genomic.gtf`
   * `GCF_000001405.40_GRCh38.p14_splice_sites.txt`
## Understanding the run_star SBATCH
```
python3 run_star.py --input 7KO-Cyto-BS_processed_fastqs --genomeDir /home/lingsn/scratch/star/star_hg38 --runThreadN=12
```
* **--input:** Absolute or relative path to folder containing merged, paired, and unpaired fastqs.
* **--genomeDir:** Path to hg38 genome index. If you are using the pre-built index, you can directly use `/home/lingsn/scratch/star/star_hg38`
* **--runThreadN=n:** Number of threads. By default n=12, and it is recommended to stay within the range of 8-12 threads for optimal results. However, if you choose to change the number of threads, then `#SBATCH --cpus-per-task=n` must also be changed accordingly within `run_star.sbatch`. 
## Additional information
* `star_index`
  * By default, overhang is set to 100 in `star_index.sbatch`. However, you can specify a custom overhang value if necessary. For more information, read about the `--sjdbOverhang` option [here](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).
  * If you are manually building the STAR index, the process may take up to 9 hours to complete.
* `run_star`
  * This script will manually reverse complement unpaired R2 fastqs prior to STAR alignment.
  * The walltime and memory in `run_star.sbatch` will need to be adjusted depending on the total size of any two data folders you have (since the SBATCH script can only run 2 concurrent tasks at once). For example, two folders that are 3.0 GB each could be processed with 50 min and 64 GB memory.
