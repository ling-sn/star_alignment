#!/usr/bin/env python
'''
Created April 2023
@author: chasew and chat GPT
Runs hisat2 with specified parameters on a folder of paired fastqs and generates gene counts for differential expression analysis

'''
import os, sys, argparse
import subprocess
import gzip
import shutil
import pandas as pd


def complement_base(base):
	"""Return the complement of a single DNA base"""
	if base == 'A':
		return 'T'
	elif base == 'T':
		return 'A'
	elif base == 'C':
		return 'G'
	elif base == 'G':
		return 'C'
	else:
		return base
		
def reverse_complement_fastq(input_file, output_file):
	# Open input and output files
	with open(input_file, 'r') as in_f, open(output_file, 'w') as out_f:
		# Loop through input file by 4 lines at a time
		for i, line in enumerate(in_f):
			if i % 4 == 1:  # Sequence line
				# Reverse complement the sequence
				sequence = line.strip()
				reverse_complement = ''.join([complement_base(base) for base in sequence[::-1]])
				out_f.write(reverse_complement + '\n')
			elif i % 4 == 3:  # Quality score line
				# Reverse quality scores
				quality_scores = line.strip()
				out_f.write(quality_scores[::-1] + '\n')
			else:  # Non-sequence and non-quality score lines
				out_f.write(line)

def reverse_complement_fastq_gz(input_file, output_fastq):
	# Decompress the input fastq.gz file using gunzip, leaving the original behind
	subprocess.run(f'gunzip -k {input_file}', check=True, shell=True)

	unzipped_file = os.path.splitext(input_file)[0]
	reverse_complement_fastq(unzipped_file, output_fastq)

	# Compress the output fastq file using gzip
	subprocess.run(f'gzip {output_fastq}', check=True, shell=True)

	# Check if the output file was created
	if not os.path.exists(output_fastq+".gz"):
		raise RuntimeError(f"Failed to create output file {output_file + '.gz'}")

	else:
		os.remove(unzipped_file)

def hisat2_align_reads(genome_index, output_sam, reads_directory, num_threads=2, directional=True, transcriptome_splicesites = None):
	"""
	Aligns paired-end RNA-Seq reads to a genome and a transcriptome simultaneously using Hisat2.
	
	Parameters:
	genome_index (str): Path to the Hisat2 index file for the genome
	transcriptome_splicesites (str): Path to the file containing known splice sites from the transcriptome
	output_sam (str): Path to the output SAM file
	reads1 (str): Path to the first input fastq file with forward reads
	reads2 (str): Path to the second input fastq file with reverse reads
	num_threads (int): Number of threads to use for alignment (default=1)
	
	Returns:
	None
	"""

	fastqs = os.listdir(reads_directory)
	fastqs.sort() # should give [html, json, merged, unmerged1, unmerged2, unpaired1, unpaired2]

	if directional:
		strandedness = "RF"

		for filename in fastqs:
			if "_unpaired_R2" in filename:
				unpaired_R2 = os.path.join(reads_directory,os.path.splitext(filename.replace("unpaired_","unpairedrc_"))[0])
				reverse_complement_fastq_gz(os.path.join(reads_directory,filename), unpaired_R2)

		fastqs[6] = os.path.basename(unpaired_R2) + ".gz"
	
	else:
		strandedness = "FR"

	reads1 = os.path.join(reads_directory,fastqs[3])
	reads2 = os.path.join(reads_directory,fastqs[4])
	single_end_reads = [os.path.join(reads_directory,fastqs[2]),os.path.join(reads_directory,fastqs[5]),os.path.join(reads_directory,fastqs[6])]

	hisat2_cmd = ['hisat2',
				  '--dta-cufflinks',
				  '--rna-strandness', strandedness]
	
	if transcriptome_splicesites is not None:
		hisat2_cmd += ['--known-splicesite-infile', transcriptome_splicesites]

	hisat2_cmd += ['-p', str(num_threads),
			'-x', genome_index,
			'-1', reads1,
			'-2', reads2,
			'-U', ",".join(single_end_reads),
			'-S', output_sam]

	hisat2_cmd = " ".join(hisat2_cmd)

	subprocess.run(hisat2_cmd, check=True, shell=True)

	os.remove(os.path.join(reads_directory,fastqs[6]))

def convert_sam_to_bam(sam_file, output_dir):
	# Create the output directory if it doesn't exist
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	# Create a subdirectory within the output directory using the base name of the SAM file
	sam_base_name = os.path.splitext(os.path.basename(sam_file))[0]
	output_subdir = os.path.join(output_dir, sam_base_name)
	os.makedirs(output_subdir)

	# Define the output BAM file name
	bam_file = os.path.join(output_subdir, sam_base_name + ".bam")

	# Run samtools to convert SAM to BAM and sort the BAM file
	sam_to_bam_cmd = " ".join(["samtools", "view", "-b", "-o", bam_file, sam_file])
	sort_bam_cmd = " ".join(["samtools", "sort", "-o", bam_file, bam_file])
	subprocess.run(sam_to_bam_cmd, check=True, shell=True)
	subprocess.run(sort_bam_cmd, check=True, shell=True)

	# Generate BAM index file
	bam_index_cmd = " ".join(["samtools", "index", bam_file])
	subprocess.run(bam_index_cmd, check=True, shell=True)

	# Delete the original SAM file
	os.remove(sam_file)

	get_unpaired_bam_cmd = " ".join(["samtools", "view", "-b", "-F", "1", bam_file, ">", os.path.splitext(bam_file)[0] + "_unpaired.bam"])
	get_paired_bam_cmd = " ".join(["samtools", "view", "-b", "-f", "1", bam_file, ">", os.path.splitext(bam_file)[0] + "_paired.bam"])

	subprocess.run(get_unpaired_bam_cmd, check=True, shell=True)
	subprocess.run(get_paired_bam_cmd, check=True, shell=True)

	return bam_file

def run_feature_counts(bam_file, annotation, output_file, feature_type="gene",num_threads=2, directional=True,paired=True,remove_original=False):

	# featureCounts -a <annotation_file.gtf> -o counts.txt -T 4 -s 2 <aligned_reads.bam>
	
	fc_cmd = ['featureCounts','-a', annotation,'-o', output_file,'-t', feature_type, '-T', str(num_threads)]


	# The -s flag denotes strandedness. a value of 2 should be used for directional NEBNext libraries, while 0 should be used for nondirectional libraries
	
	if directional:
		direction = '2'
	else:
		direction = '0'

	fc_cmd += ['-s', direction]

	if paired:
		fc_cmd += ['-p', '--countReadPairs']

	fc_cmd += [bam_file]

	fc_cmd = " ".join(fc_cmd)

	subprocess.run(fc_cmd, check=True, shell=True)

	if remove_original:
		os.remove(bam_file)

def merge_counts_files(input_file1, input_file2, output_file):
	# Read input count files into pandas DataFrames
	df1 = pd.read_csv(input_file1, sep='\t', comment="#")
	df2 = pd.read_csv(input_file2, sep='\t', comment="#")

	# Extract the count column name from the last column header
	count_column1 = df1.columns[-1]
	count_column2 = df2.columns[-1]
	
	# Merge the counts column from both DataFrames based on the 'Geneid' column
	merged_df = pd.merge(df1, df2[['Geneid', count_column2]], on='Geneid', how='outer')

	# Add the counts from both files and assign the sum to the 'counts' column
	merged_df['counts'] = merged_df[count_column1] + merged_df[count_column2]
	
	# Drop the unnecessary columns
	merged_df.drop([count_column1, count_column2], axis=1, inplace=True)
	
	# Write the merged DataFrame to the output count file
	merged_df.to_csv(output_file, sep='\t', index=False)

	# remove redundant originals
	os.remove(input_file1)
	os.remove(input_file2)

if __name__ == '__main__': # allows another python script to import the functions

	parser = argparse.ArgumentParser(description="Runs fastP", usage='%(prog)s [options] --input in_dir --index ind_path',add_help=False) # a description of the function

	required = parser.add_argument_group('Required Input', 'Specify input directory, genome indices, and output directory.')

	required.add_argument("--input",help="path to directory with merged, paired and unpaired fastqs.", required=True)
	required.add_argument("--index",help="path to genome index.", required=True)
	required.add_argument("--annotation",help="path to GTF annotation file.", required=True)
	required.add_argument("--output",help="name of alignment output directory.", required=True)

	data_opt = parser.add_argument_group('Basic Arguments', 'These options can be used to change how the script is generally run.')
	data_opt.add_argument("-C","--threads",default=2,help="Specify number of cpus. Default is 2.")
	data_opt.add_argument("-D","--directional",action="count",help="Call this option if your reads are nondirectional")
	data_opt.add_argument("-I","--include_introns",action="count",help="Call this option if you expect/want to include intronic reads in your gene-level counts. Default is exonic-only")
	data_opt.add_argument("-S","--ss",default=None,help="Specify a known splice site txt file if not using a transcript-inclusive genome")

	help_opt = parser.add_argument_group('Help')
	help_opt.add_argument('-h', '--help', action="help", help="show this help message and exit")

	args = parser.parse_args()

	samout = os.path.basename(args.input)+".sam"

	hisat2_align_reads(args.index, samout, args.input, num_threads=args.threads, directional=(not args.directional), transcriptome_splicesites = args.ss)

	bam_file = convert_sam_to_bam(samout, args.output)

	bam_base = os.path.splitext(bam_file)[0]

	if args.include_introns:
		feature_type = 'gene'
	else:
		feature_type = 'exon'
	
	run_feature_counts(bam_base + "_unpaired.bam", args.annotation, bam_base + "_unpaired_counts.txt", feature_type=feature_type, num_threads=args.threads, directional=(not args.directional),paired=False, remove_original=True)
	run_feature_counts(bam_base + "_paired.bam", args.annotation, bam_base + "_paired_counts.txt", feature_type=feature_type, num_threads=args.threads, directional=(not args.directional),paired=True, remove_original=True)

	merge_counts_files(bam_base + "_paired_counts.txt", bam_base + "_unpaired_counts.txt", bam_base + "_counts.txt")

