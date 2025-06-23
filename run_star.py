from pathlib import Path
import traceback
import argparse
import subprocess

class StarAligner:
    def __init__(self, output_dir, input_name):
        self.output_folder = output_dir/input_name
        self.r1_filename = None
        self.r2_filename = None
        
    def single_reads(self, runThreadN, single, star_index, file):
        """
        Align single-end reads (merged/unpaired)
        """
        single_str = ",".join(single)
        prefix = self.output_folder/"single"
        
        try:
            cmd = ["STAR", "--runThreadN", str(runThreadN),
                   "--runMode", "alignReads",
                   "--readFilesIn", str(single_str),
                   "--readFilesCommand", "gunzip", "-c",
                   "--genomeDir", str(star_index),
                   "--outFileNamePrefix", str(prefix),
                   "--outSAMtype", "BAM", "SortedByCoordinate"]
            result = subprocess.run(cmd, 
                                    check = True, 
                                    capture_output = True, 
                                    text = True)
        except subprocess.CalledProcessError as e: ## error handling
            print(f"Failed to align merged or unpaired fastq file {file.name}: {e}")
            print("STDERR:", e.stderr)
            print("STDOUT:", e.stdout)
            traceback.print_exc()
            raise
        return result

    def paired_reads(self, runThreadN, r1_str, r2_str, star_index, file):
        """
        Align paired-end reads (unmerged)
        """
        paired_str = " ".join([r1_str, r2_str])
        prefix = self.output_folder/"paired"

        try:
            cmd = ["STAR", "--runThreadN", str(runThreadN),
                   "--runMode", "alignReads",
                   "--readFilesIn", str(paired_str),
                   "--readFilesCommand", "gunzip", "-c",
                   "--genomeDir", str(star_index),
                   "--outFileNamePrefix", str(prefix),
                   "--outSAMtype", "BAM", "SortedByCoordinate"]
            result = subprocess.run(cmd, 
                                    check = True, 
                                    capture_output = True, 
                                    text = True)
        except subprocess.CalledProcessError as e: ## error handling
            print(f"Failed to align unmerged fastq file {file.name}: {e}")
            print("STDERR:", e.stderr)
            print("STDOUT:", e.stdout)
            traceback.print_exc()
            raise
        return result
    
    def merge_bam(self, input_name):
        """
        Merges all .bam files, then 
        sorts and indexes into .bai
        """
        merged_bam = self.output_folder/f"{input_name}.bam"
        bam_list = [*self.output_folder.glob("*.bam")] # detect .bam files
        rm_list = [*self.output_folder.glob("*out.bam")]

        try:
            subprocess.run(["samtools", "merge", ## merge all .bam files into one
                            str(merged_bam), *map(str, bam_list)],
                            check = True, 
                            capture_output = True,
                            text = True)
            subprocess.run(["samtools", "index", str(merged_bam)], ## create .bai from .bam
                            check = True,
                            capture_output = True,
                            text = True)
            subprocess.run(["rm", *map(str, rm_list)], ## remove original .bam files
                            check = True,
                            capture_output = True,
                            text = True)
        except subprocess.CalledProcessError as e: ## error handling
            print(f"Failed to create {merged_bam.name} and convert to .bai: {e}")
            print("STDERR:", e.stderr)
            print("STDOUT:", e.stdout)
            traceback.print_exc()
            raise

def star_pipeline(folder_path, genomeDir, runThreadN):
    input_dir = Path(folder_path)
    star_index = Path(genomeDir)
    output_dir = input_dir/"alignments"
    output_dir.mkdir(exist_ok=True)
    input_name = input_dir.name

    ## initialize class
    aligner = StarAligner(output_dir, input_name)

    for subfolder in input_dir.iterdir(): ## amount of subfolders = number of replicates
        if subfolder.is_dir():
            single = []
            paired_r1 = []
            paired_r2 = []

            for file in subfolder.glob("*.fastq.gz"): ## iterate through indiv. files in subfolder
                try:
                    ## run star alignment functions
                    if "_merged" in file.name or "_unpaired" in file.name:
                        single.append(file)
                        aligner.single_reads(runThreadN, single, star_index, file)
                    elif "_unmerged" in file.name:
                        for r1_file in subfolder.glob("*unmerged_R1*"):
                            paired_r1.append(r1_file)
                            paired_r2.append(r1_file.with_name(r1_file.name.replace("_R1_", "_R2_")))
                        r1_str = ",".join(paired_r1)
                        r2_str = ",".join(paired_r2)
                        aligner.paired_reads(runThreadN, r1_str, r2_str, star_index, file)                    
                except Exception as e:
                    print(f"Failed to align {file.name} with STAR and produce .bam files: {e}")
                    traceback.print_exc()
                    continue
            
            ## merge bam files, convert to bai, & remove old files
            aligner.merge_bam(input_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Runs STAR alignment.")
    parser.add_argument("--input", help = "Path to directory with merged, paired, and unpaired fastqs", required = True)
    parser.add_argument("--genomeDir", help = "Path to genome index", required = True)
    parser.add_argument("--runThreadN", type = int, default = 12, help = "Number of CPU cores (default: 12)")
    args = parser.parse_args()

    print("Starting STAR alignment pipeline...")
    star_pipeline(args.input, args.genomeDir, args.runThreadN)
    print("Pipeline finished.")