## <ins>**PART I: Running STAR alignment on processed fastqs**</ins>
### Necessary files
<img src="https://github.com/user-attachments/assets/607f335f-5073-4c4b-b3d6-2e158d59eaed" width="400"/>

* STAR index
  * **Option 1 (Manual):** First, create the `star_index_hg38` folder according to the instructions below (see "_Creating the star_index_hg38 folder_"). Then, copy over `star_index.py` and `star_index.sbatch` and run the SBATCH file.
  * **Option 2 (Pre-Built; Recommended):** Skip STAR index creation by directly using `--genomeDir ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38` in your `run_star.sbatch` file.
* `create_env.sbatch`
* `run_star.py` and `run_star.sbatch`
* `tagXSstrandedData.awk`
### Instructions
1. Run `create_env.sbatch` to create the RNA-STAR conda environment
3. Activate conda environment via `conda activate RNA-STAR`
4. Run `star_index.sbatch` one time to build hg38 genome index (_optional_)
5. Edit `run_star.sbatch` to match your experiments
6. Run `run_star.sbatch` to align each input folder using STAR
### Tools used in STAR alignment script
* **STAR** is used to align single-end reads (merged/unpaired) and paired reads (unmerged) according to the hg38 genome index
* **samtools** is used to merge single-end and paired .bam files together
### When do I use this pipeline?
This is used after running the fastp script on your raw data (fastq files). 
### Creating the star_index_hg38 folder
Read the following if you want to manually build the STAR index. Otherwise, feel free to skip this section.
1. Create an empty folder called `star_index_hg38` in the working directory containing your data folders
2. On Great Lakes Cluster, navigate to `~/umms-RNAlabDATA/Software/genome_indices/hisat2_hg38/hg38p14_tran/` and copy the following 4 files into `star_index_hg38`:
   * `GCF_000001405.40_GRCh38.p14_exons.txt`
   * `GCF_000001405.40_GRCh38.p14_genomic.fa`
   * `GCF_000001405.40_GRCh38.p14_genomic.gtf`
   * `GCF_000001405.40_GRCh38.p14_splice_sites.txt`
### Understanding the run_star SBATCH
```
python3 run_star.py --input 7KO-Cyto-BS_processed_fastqs --genomeDir ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38 --runThreadN=12
```
* **--input:** Name of folder containing merged, paired, and unpaired fastqs. DO NOT INPUT A PATH.
* **--genomeDir:** Path to hg38 genome index. If you are using the pre-built index, you can directly use `~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38`
* **--runThreadN=n:** Number of threads. By default n=8, and it is recommended to stay within the range of 8-12 threads for optimal results. However, if you choose to change the number of threads, then `#SBATCH --cpus-per-task=n` must also be changed accordingly within `run_star.sbatch`. 
### Additional information
* `star_index`
  * By default, overhang is set to 100 in `star_index.sbatch`. However, you can specify a custom overhang value if necessary. For more information, read about the `--sjdbOverhang` option [here](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).
  * If you are manually building the STAR index, the process may take up to 9 hours to complete.
* `run_star`
  * This script will manually reverse complement unpaired R2 fastqs prior to STAR alignment.
  * The walltime in `run_star.sbatch` will need to be adjusted depending on the size of your data. For reference, it takes approximately 10m to process 3 GB of data, and 3h15m to process 42 GB of data.
### Citations
* Reverse complement code adapted from `run_hisat2.py` by Chase Weidmann
* `tagXSstrandedData.awk` sourced from [STAR Aligner](https://github.com/alexdobin/STAR/blob/master/extras/scripts/tagXSstrandedData.awk) by Alex Dobin
---
## <ins>**PART II: Realigning BAM files after STAR alignment**</ins>
### Necessary files
* `realignGap.py` and `realignGap.sbatch`
### Tools used in STAR realignment script
* **pysam** is used to iterate through BAM files and write in reads
* **parasail** is used to align reads with the Smith-Waterman algorithm
> "Alignment software usually uses a seed alignment algorithm to increase alignment speed; however, this also affects pairwise alignment accuracy, especially for bases near deletion signatures. To solve this, we integrated the Smith-Waterman local alignment algorithm into the pipeline for realignment. Reads that contained any mismatch, deletion, insertion, soft-clip or splicing were further processed by the realignment tool in the BID-pipe package. By setting the penalty of gap open and gap extension as −3 and −2, respectively, deletion signatures can have a higher priority in the alignment" (_Zhang et al., 520_). 
### When do I use this pipeline?
This is used after running the STAR alignment script (`run_star.py`). Start from the working directory that contains the `alignments` folder.
### Understanding the realignGap SBATCH
```
python3 realignGap.py --folder_name 7KO-Cyto-BS_processed_fastqs --fasta_dir ~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa --discard
```
* **--folder_name:** Name of processed_fastqs folder that you wish to realign. DO NOT INPUT A PATH.
* **--fasta_dir:** Path to FASTA file used to create hg38 genome index. If you used the pre-built index in PART I, you can directly use `~/umms-RNAlabDATA/Software/genome_indices/star_index_hg38/GCF_000001405.40_GRCh38.p14_genomic.fa`
* **--discard:** Writes discarded reads into a file for debugging. This is disabled by default with `--no-discard`, but it can be enabled with `--discard`
### Citations
* Realignment code adapted from [realignGap](https://github.com/y9c/pseudoU-BIDseq/blob/main/bin/realignGap) by Ye Chang
* Zhang et al. BID-seq for transcriptome-wide quantitative sequencing of mRNA pseudouridine at base resolution. _Nature Protocols_ 19, 517–538 (2024). https://doi.org/10.1038/s41596-023-00917-5

