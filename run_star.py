from pathlib import Path
import traceback
import argparse
import subprocess

class StarAligner:
    def __init__(self, output_dir, input_name):
        self.output_folder = output_dir/input_name
        self.r1_filename = None
        self.r2_filename = None
        
    def merged_reads(self, runThreadN, merged, star_index):
        """
        Align single-end reads (merged)
        """
        merged_str = ",".join(merged)
        prefix = self.output_folder/"merged"
        
        try:
            cmd = ["STAR", "--runThreadN", str(runThreadN),
                   "--runMode", "alignReads",
                   "--readFilesIn", str(merged_str),
                   "--readFilesCommand", "gunzip", "-c",
                   "--genomeDir", str(star_index),
                   "--outFileNamePrefix", str(prefix),
                   "--outSAMtype", "BAM", "SortedByCoordinate"]
            result = subprocess.run(cmd, 
                                    check = True, 
                                    capture_output = True, 
                                    text = True)
        except subprocess.CalledProcessError as e: ## error handling
            print(f"Failed to align merged fastqs with STAR: {e}")
            print("STDERR:", e.stderr)
            print("STDOUT:", e.stdout)
            traceback.print_exc()
            raise
        return result

    def unpaired_reads(self, runThreadN, unpaired, star_index):
        """
        Align single-end reads (unpaired)
        """
        unpaired_str = ",".join(unpaired)
        prefix = self.output_folder/"unpaired"
        
        try:
            cmd = ["STAR", "--runThreadN", str(runThreadN),
                   "--runMode", "alignReads",
                   "--readFilesIn", str(unpaired_str),
                   "--readFilesCommand", "gunzip", "-c",
                   "--genomeDir", str(star_index),
                   "--outFileNamePrefix", str(prefix),
                   "--outSAMtype", "BAM", "SortedByCoordinate"]
            result = subprocess.run(cmd, 
                                    check = True, 
                                    capture_output = True, 
                                    text = True)
        except subprocess.CalledProcessError as e: ## error handling
            print(f"Failed to align unpaired fastqs with STAR: {e}")
            print("STDERR:", e.stderr)
            print("STDOUT:", e.stdout)
            traceback.print_exc()
            raise
        return result

    def paired_reads(self, runThreadN, paired_r1, paired_r2, star_index):
        """
        Align paired-end reads (unmerged)
        """
        r1_str = ",".join(paired_r1)
        r2_str = ",".join(paired_r2)
        prefix = self.output_folder/"paired"

        try:
            cmd = ["STAR", "--runThreadN", str(runThreadN),
                   "--runMode", "alignReads",
                   "--readFilesIn", str(r1_str), str(r2_str),
                   "--readFilesCommand", "gunzip", "-c",
                   "--genomeDir", str(star_index),
                   "--outFileNamePrefix", str(prefix),
                   "--outSAMtype", "BAM", "SortedByCoordinate"]
            result = subprocess.run(cmd, 
                                    check = True, 
                                    capture_output = True, 
                                    text = True)
        except subprocess.CalledProcessError as e: ## error handling
            print(f"Failed to align unmerged fastq files: {e}")
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
        bam_list = [*self.output_folder.glob("*out.bam")] # detect .bam files
        rm_list = [*self.output_folder.glob("*out.bam")] # can also rm: *self.output_folder.glob("*.out"), *self.output_folder.glob("*out.tab")

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

def collect_files(subfolder, match_pattern, list):
    for i in subfolder.glob(match_pattern): ## used for finding files and appending them to a list; avoids redundant for loop later
        str_name = str(i)
        list.append(str_name)

def star_pipeline(folder_name, genomeDir, runThreadN):
    current_path = Path.cwd()
    input_dir = current_path/folder_name
    star_index = Path(genomeDir)
    output_dir = current_path/"alignments"
    output_dir.mkdir(exist_ok=True)
    input_name = input_dir.name

    ## initialize class
    aligner = StarAligner(output_dir, input_name)

    for subfolder in input_dir.iterdir(): ## amount of subfolders = number of replicates
        if subfolder.is_dir():
            merged = []
            unpaired = []
            paired_r1 = []
            paired_r2 = []

            for file in subfolder.glob("*.fastq.gz"): ## iterate through files and add to corresponding lsits
                try:
                    ## run star alignment functions
                    if "_merged" in file.name:
                        collect_files(subfolder, "*_merged*", merged)
                    elif "_unpaired" in file.name:
                        collect_files(subfolder, "*_unpaired*", unpaired)
                    elif "_unmerged" in file.name:
                        for r1_file in subfolder.glob("*_unmerged_R1*"):
                            r1_str_name = str(r1_file)
                            r2_file = r1_file.with_name(r1_file.name.replace("_R1_", "_R2_"))
                            r2_str_name = str(r2_file)
                            paired_r1.append(r1_str_name)
                            paired_r2.append(r2_str_name)           
                except Exception as e:
                    print(f"Failed to align {file.name} with STAR and produce .bam files: {e}")
                    traceback.print_exc()
                    continue
            
            ## run star alignment
            aligner.merged_reads(runThreadN, merged, star_index)
            aligner.unpaired_reads(runThreadN, unpaired, star_index)
            aligner.paired_reads(runThreadN, paired_r1, paired_r2, star_index) 

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