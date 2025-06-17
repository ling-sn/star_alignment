from pathlib import Path
import subprocess
import traceback
import re

def star_index():
    """
    Before proceeding, the `star_hg38` folder should
    be in the starting directory containing the data 
    folders you want to process.
    """
    current_path = Path.cwd()
    genome_dir = f"{current_path}/star_hg38"
    genome_fasta = f"{genome_dir}/GCF_000001405.40_GRCh38.p14_genomic.fa"
    gtf_file = f"{genome_dir}/GCF_000001405.40_GRCh38.p14_genomic.gtf"
    threads = 2

    if len(list(Path(genome_dir).glob("*"))) != 4: ## ensures this only runs if index hasn't been created yet
        try:
            print("""--sjdbOverhang asks for the number READLENGTH-1, where READLENGTH is the length of your 
                  sequencing reads. For example, the ideal value for Illumina 2x100b paired-end reads is 100-1=99. 
                  For reads of varying length, the ideal value is max(ReadLength)-1. However, in most cases, the 
                  default value 100 will work as well as the ideal value.""")
            while True:
                overhang = int(input("Please input the overhang value.\nIf nothing is entered, the default value 100 will be used."))
                if re.search(r"\D", overhang):
                    print("Your input should be a number.")
                if overhang == "":
                    overhang = 100
                    break
                else:
                    break

            cmd = ["STAR", "--runThreadN", str(threads),
                   "--runMode", "genomeGenerate",
                   "--genomeDir", str(genome_dir),
                   "--genomeFastaFiles", str(genome_fasta),
                   "--sjdbGTFfile", str(gtf_file),
                   "--sjdbOverhang", str(overhang)]
            subprocess.run(cmd, 
                           check = True,
                           capture_output = True, 
                           text = True)
            
        except subprocess.CalledProcessError as e: ## error handling
            print(f"Failed to build STAR index: {e}")
            print("STDERR:", e.stderr)
            print("STDOUT:", e.stdout)
            traceback.print_exc()
            raise
    else:
        print("Index already exists.")

if __name__ == "__main__":
    print("Creating STAR index...")
    star_index()
    print("Index created.")