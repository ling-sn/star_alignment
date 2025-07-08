from pathlib import Path
import traceback
import argparse
import subprocess

class StarAligner:        
    def merged_reads(self, runThreadN, merged, star_index, processed_folder):
        """
        Align single-end reads (merged)
        """
        merged_str = ",".join(merged)
        prefix = processed_folder/"merged"
        
        try:
            cmd = ["STAR", "--runThreadN", str(runThreadN),
                   "--runMode", "alignReads",
                   "--readFilesIn", str(merged_str),
                   "--readFilesCommand", "gunzip", "-c",
                   "--genomeDir", str(star_index),
                   "--outFileNamePrefix", str(prefix),
                   "--outSAMtype", "BAM", "SortedByCoordinate",
                   "--outFilterType", "BySJout",
                   "--outSAMattributes", "NH", "NM", "AS", "MD"]
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

    def unpaired_reads(self, runThreadN, unpaired_r1, unpaired_r2, star_index, processed_folder):
        """
        Align single-end reads (unpaired)
        """
        r1_str = ",".join(unpaired_r1)
        r2_str = ",".join(unpaired_r2)
        prefix_1 = processed_folder/"unpaired_r1"
        prefix_2 = processed_folder/"unpaired_r2"
        
        try:
            cmd_1 = ["STAR", "--runThreadN", str(runThreadN),
                     "--runMode", "alignReads",
                     "--readFilesIn", str(r1_str),
                     "--readFilesCommand", "gunzip", "-c",
                     "--genomeDir", str(star_index),
                     "--outFileNamePrefix", str(prefix_1),
                     "--outSAMtype", "BAM", "SortedByCoordinate",
                     "--outFilterType", "BySJout",
                     "--outSAMattributes", "NH", "NM", "AS", "MD"]
            cmd_2 = ["STAR", "--runThreadN", str(runThreadN),
                     "--runMode", "alignReads",
                     "--readFilesIn", str(r2_str),
                     "--readFilesCommand", "gunzip", "-c",
                     "--genomeDir", str(star_index),
                     "--outFileNamePrefix", str(prefix_2),
                     "--outSAMtype", "BAM", "SortedByCoordinate",
                     "--outFilterType", "BySJout",
                     "--outSAMattributes", "NH", "NM", "AS", "MD"]
            subprocess.run(cmd_1, 
                           check = True, 
                           capture_output = True, 
                           text = True)
            subprocess.run(cmd_2, 
                           check = True, 
                           capture_output = True, 
                           text = True)
        except subprocess.CalledProcessError as e: ## error handling
            print(f"Failed to align unpaired fastqs with STAR: {e}")
            print("STDERR:", e.stderr)
            print("STDOUT:", e.stdout)
            traceback.print_exc()
            raise

    def paired_reads(self, runThreadN, paired_r1, paired_r2, star_index, processed_folder):
        """
        Align paired-end reads (unmerged)
        """
        r1_str = ",".join(paired_r1)
        r2_str = ",".join(paired_r2)
        prefix = processed_folder/"paired"

        try:
            cmd = ["STAR", "--runThreadN", str(runThreadN),
                   "--runMode", "alignReads",
                   "--readFilesIn", str(r1_str), str(r2_str),
                   "--readFilesCommand", "gunzip", "-c",
                   "--genomeDir", str(star_index),
                   "--outFileNamePrefix", str(prefix),
                   "--outSAMtype", "BAM", "SortedByCoordinate",
                   "--outFilterType", "BySJout",
                   "--outSAMattributes", "NH", "NM", "AS", "MD"]
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
    
    def tagXSstrandedData(self, awk_dir):
        """
        Add XS tags to every read in a stranded .bam file
        """
        tmp = f"{self.merged_bam.stem}_tmp.bam"

        try:
            file = subprocess.run(["samtools", "view", ## open up merged .bam
                                    "-h", str(self.merged_bam)],
                                    stdout = subprocess.PIPE,
                                    check = True, 
                                    capture_output = True,
                                    text = True)
            subprocess.run(["awk", "-v", "strType=2", ## run awk script
                            "-f", str(awk_dir)],
                            input = file.stdout, 
                            stdout = str(tmp),
                            check = True, 
                            capture_output = True,
                            text = True)
            subprocess.run(["mv", str(tmp), str(self.merged_bam)], ## rename tmp back to original filename
                            check = True, 
                            capture_output = True,
                            text = True)
        except subprocess.CalledProcessError as e: ## error handling
            print(f"Failed to add XS tags to {self.merged_bam.name}: {e}")
            print("STDERR:", e.stderr)
            print("STDOUT:", e.stdout)
            traceback.print_exc()
            raise
    
    def merge_bam(self, processed_folder, subfolder, awk_dir):
        """
        Merge all .bam files, then 
        sort and index into .bai
        """
        self.merged_bam = processed_folder/f"{subfolder.name}.bam"
        bam_list = [*processed_folder.glob("*out.bam")] # detect .bam files
        rm_list = [*processed_folder.glob("*out.bam")]

        try:
            subprocess.run(["samtools", "merge", ## merge all .bam files into one
                            str(self.merged_bam), *map(str, bam_list)],
                            check = True, 
                            capture_output = True,
                            text = True)
            self.tagXSstrandedData(awk_dir)
            subprocess.run(["samtools", "index", str(self.merged_bam)], ## create .bai from .bam
                            check = True,
                            capture_output = True,
                            text = True)
            subprocess.run(["rm", *map(str, rm_list)], ## remove original .bam files
                            check = True,
                            capture_output = True,
                            text = True)
        except subprocess.CalledProcessError as e: ## error handling
            print(f"Failed to create {self.merged_bam.name} and convert to .bai: {e}")
            print("STDERR:", e.stderr)
            print("STDOUT:", e.stdout)
            traceback.print_exc()
            raise

## MODIFIED CODE FROM run_hisat2.py

def complement_base(base):
    """
    Return complement of a single DNA base
    """
    if base == "A":
        return "T"
    elif base == "T":
        return "A"
    elif base == "C":
        return "G"
    elif base == "G":
        return "C"
    else:
        return base
    
def reverse_complement_fastq(file, output):
    """
    Manually reverse complement unpaired R2 fastq.gz file
    """
    # Decompress input fastq.gz file using gunzip
    file = Path(file)
    output = Path(output)
    unzipped_file = file.with_suffix("")
    
    try:        
        if not unzipped_file.exists():
            subprocess.run(["gunzip", "-k", str(file)],
                            check = True,
                            capture_output = True,
                            text = True)

        # Open file and reverse complement bases
        with open(unzipped_file, "r") as input_file, open(output, "w") as output_file:
            for i, line in enumerate(input_file):
                if i % 4 == 1:
                    # Reverse complement sequence
                    sequence = line.strip()
                    reverse_complement = "".join([complement_base(base) for base in sequence[::-1]])
                    output_file.write(reverse_complement + "\n")
                elif i % 4 == 3: 
                    # Reverse quality scores
                    quality_scores = line.strip()
                    output_file.write(quality_scores[::-1] + "\n") 
                else:  
                    # Non-sequence and non-quality score lines
                    output_file.write(line)

        subprocess.run(["gzip", str(output)],
                        check = True,
                        capture_output = True,
                        text = True)
    except subprocess.CalledProcessError as e:
        print(f"Failed to create output fastq {output}: {e}")
        print("STDERR:", e.stderr)
        print("STDOUT:", e.stdout)
        traceback.print_exc()
        raise
    if unzipped_file.exists():
        unzipped_file.unlink()

## MODIFIED CODE END

def collect_files(subfolder, match_pattern, list):
    for i in subfolder.glob(match_pattern): ## used for finding files and appending them to a list; avoids redundant for loop later
        str_name = str(i)
        list.append(str_name)

def star_pipeline(folder_name, genomeDir, runThreadN):
    current_path = Path.cwd()
    input_dir = current_path/folder_name
    input_name = input_dir.name
    star_index = Path(genomeDir)
    awk_dir = current_path/"tagXSstrandedData.awk"

    ## initialize class
    aligner = StarAligner()

    for subfolder in input_dir.iterdir(): ## amount of subfolders = number of replicates
        if subfolder.is_dir():
            merged = []
            unpaired_r1 = []
            unpaired_r2 = []
            paired_r1 = []
            paired_r2 = []
            processed_folder = current_path/"alignments"/input_name/f"{subfolder.name}_star"
            processed_folder.mkdir(exist_ok=True, parents=True)

            for file in subfolder.glob("*.fastq.gz"): ## iterate through files and add to corresponding lsits
                try:
                    ## run star alignment functions
                    if "_merged" in file.name:
                        collect_files(subfolder, "*_merged*", merged)
                    elif "_unpaired" in file.name:
                        if "_unpaired_R2" in file.name:
                            output = file.with_name(file.name.replace("_unpaired_", "_unpairedrc_"))
                            output_fastq = output.with_suffix("")
                            reverse_complement_fastq(file, output_fastq)
                            for r1_file in subfolder.glob("*_unpaired_R1*"):
                                r1_str_name = str(r1_file)
                                r2_str_name = str(output_fastq) + ".gz"
                                unpaired_r1.append(r1_str_name)
                                unpaired_r2.append(r2_str_name)
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
            aligner.merged_reads(runThreadN, merged, star_index, processed_folder)
            aligner.unpaired_reads(runThreadN, unpaired_r1, unpaired_r2, star_index, processed_folder)
            aligner.paired_reads(runThreadN, paired_r1, paired_r2, star_index, processed_folder) 

            ## merge bam files, convert to bai, & remove old files
            aligner.merge_bam(processed_folder, subfolder, awk_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Runs STAR alignment.")
    parser.add_argument("--input", help = "Path to directory with merged, paired, and unpaired fastqs", required = True)
    parser.add_argument("--genomeDir", help = "Path to genome index", required = True)
    parser.add_argument("--runThreadN", type = int, default = 8, help = "Number of CPU cores (default: 12)")
    args = parser.parse_args()

    print("Starting STAR alignment pipeline...")
    star_pipeline(args.input, args.genomeDir, args.runThreadN)
    print("Pipeline finished.")