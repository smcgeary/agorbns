################################################################################
#MakeTransfectionCountData.py
################################################################################
import imp # Used to import general.py
imp.load_source("general",
				("/lab/bartel1_ata/mcgeary/computation/"
				 "AgoRBNS/general/general.py"
				)
			   )
imp.load_source("RBNS_methods",
				("/lab/bartel1_ata/mcgeary/computation/"
				 "AgoRBNS/general/RBNS_methods.py"
				)
			   )
from general import *
from RBNS_methods import *
from sitetypes import get_seq_site_map

# I first used the above address, but the bottom one gives near-perfect (but
# not perfect) agreement with Kathy's files, in the folder:
# /lab/solexa_bartel/klin/miRNA_models_data/transfections/hela/compiled/ .

base_exp_dir = "/lab/solexa_bartel/mcgeary/transfections/"
star_index = ("/nfs/genomes/human_gp_feb_09_no_random/STAR_2.6.0/"
			  "GRCh37.75_overhang_100/")
temp_dir   = "%sSTAR_alignment_temp" %base_exp_dir
raw_files_dir = "/lab/solexa_bartel/mcgeary/GEO_McGearyLin2019/raw_files/"




def main():
	time_start = time.time()
	# Parse command line arguments:
	arguments = ["cell_line"]
	args = parse_arguments(arguments)
	print(args)
	(cell_line,) = args
	print(cell_line)
	# Assign the cell-line argument.
	annot_file = "%smetadata/%s_annotations.gtf" %(base_exp_dir, cell_line)
	files_dir = "Repression/%s_samples.txt" %(cell_line)
	print(files_dir)
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
		out_dir = "%s%s/STAR_htseq/" %(base_exp_dir, cell_line)
		print(fastq_file)
		print(fastq_file)
		print(star_index)
		print(annot_file)
		print(temp_dir)
		print(out_dir)
		# Define the string that gives the script to be called by os.system:
		script = "bsub -q 14 bash Repression/NEXTflex_pipeline.sh %s %s %s %s %s" %(
			fastq_file, star_index, annot_file, temp_dir, out_dir
		)
		os.system(script)
	# End of loop, end of function.
	return()





if __name__ == "__main__":
	main()
