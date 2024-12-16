################################################################################
#repression_test.py
################################################################################
import os

# This is a test file in which I develop code by which to reconstruct the
# fitting of the biochemical model/TargetScan to the HeLa and HEK93 data, with
# the ultimate purpose of checking whether the repression score can be
# improved with the insights gained from the three-prime AGO-RBNS analysis.

# Use Star to align the reads to the genome, in order to make the count data
# table.

## CONSTANT VARIABLE ###########################################################

# star_index = "/nfs/genomes/human_gp_feb_09/STAR_2.6.0/GRCh37.75_overhang_100/"
# I first used the above address, but the bottom one gives near-perfect (but
# not perfect) agreement with Kathy's files, in the folder:
# /lab/solexa_bartel/klin/miRNA_models_data/transfections/hela/compiled/ .

base_exp_dir = "/lab/solexa_bartel/mcgeary/transfections/"
star_index = "/nfs/genomes/human_gp_feb_09_no_random/STAR_2.6.0/GRCh37.75_overhang_100/"
annot_file = "%smetadata/annotations.gtf" %base_exp_dir
temp_dir   = "%sSTAR_alignment_temp" %base_exp_dir
raw_files_dir = "/lab/solexa_bartel/mcgeary/GEO_McGearyLin2019/raw_files/"


files_dir = "Repression/HeLa_samples.txt"


def main():	
	cell_line = "HeLa"
	with open(files_dir) as file:
		samples = file.readlines()
	for sample in samples:
		print("______________________________________")
		print(sample)
		mirna, rep, batch = sample.strip().split("\t")
		print(mirna)
		print(rep)
		print(batch)

		raw_file_name = "%s_%s_transfection_rep%s_batch%s.txt.gz" %(
			cell_line, mirna, rep, batch
		)

		fastq_file = raw_files_dir + raw_file_name
		print(fastq_file)
		out_dir = "%s%s/STAR_htseq/" %(base_exp_dir, cell_line)
		# Define the string that gives the script to be called by os.system:
		script = "bsub -q 14 bash Repression/NEXTflex_pipeline.sh %s %s %s %s %s" %(
			fastq_file, star_index, annot_file, temp_dir, out_dir
		)
		os.system(script)
	return()




if __name__ == "__main__":
    main()


# fastq_file = "/lab/solexa_bartel/mcgeary/temp_GEO/HeLa_lsy-6_rep2_batch5.fq"
# out_dir = "/lab/solexa_bartel/mcgeary/temp_output_repression_no_random"


# print(script)

