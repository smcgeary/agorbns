################################################################################
#MakeSiteTypeReadFiles.py
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



##########UTR READING###########################################################
hsa_UTR_dict = dict()
with open("/lab/bartel4_ata/kathyl/RNA_Seq/analysis/data/final/utr3.txt") as file_in:
	line = file_in.readline()
	while line:
		line_strip = line.strip().split("\t")
		hsa_UTR_dict[line_strip[0]] = line_strip[1]
		line = file_in.readline()

hsa_ORF_dict = dict()
with open("/lab/bartel4_ata/kathyl/RNA_Seq/analysis/data/final/orf_lengths.txt") as file_in:
	line = file_in.readline()
	line = file_in.readline()
	while line:
		line_strip = line.strip().split("\t")
		hsa_ORF_dict[line_strip[0]] = float(line_strip[1])/1000
		line = file_in.readline()




#Begin script

def main():
	missing_genes = 0
	count_df = pd.read_csv("RepressionData/Lin-Shi_transfection_data/counts.txt",
	                     sep="\t", header=None, index_col=0, names=None)
	# Make the output column names, incorporating mirna, rep, and batch.
	mirnas_original = list(count_df.loc["mirna"])
	reps = list(count_df.loc["rep"])
	batches = list(count_df.loc["batch"])
	col_strings = ["%s_rep%s_b%s" %(i[0], i[1], i[2]) for i in
				   zip(mirnas_original, reps, batches)]
	# Convert the mirna strings from Kathy to my representation, and then to
	# their seed sequences:
	mirnas = [re.sub("mir", "miR-", mirna) for mirna in mirnas_original]
	mirnas = [re.sub("lsy6", "lsy-6", mirna) for mirna in mirnas]
	mirnas = [re.sub("let7", "let-7a", mirna) for mirna in mirnas]
	mirna_seeds = [Mirna(mirna)["6mer"] for mirna in mirnas]
	# Make lists for the rpkms and the orf lenght / utr sequences.

	# Convert each row of count data to rpkm by dividing by the utr in kb,
	# assuming that it has both an orf length and a UTR sequence in the
	# dictionaries generated above.
	rpkm_list = []
	orf_utr_list = []
	for row in count_df.itertuples():
		index = row.Index
		if index not in ["mirna", "rep", "batch"]:
			if index in hsa_ORF_dict and index in hsa_UTR_dict and index != "NM_001025204":
				utr = hsa_UTR_dict[index]
				orf_length = hsa_ORF_dict[index]
				out_row = [index, orf_length, utr]
				orf_utr_list.append(out_row)
				rpkm_row = [(int(i) + 1)/orf_length for i in row[1:]]
				rpkm_list.append(rpkm_row)
			else:
				print(index)
				print(index in hsa_ORF_dict)
				print(index in hsa_UTR_dict)
				missing_genes += 1
	print(missing_genes)

	# Make a new count data frame that only has the genes which contain UTRs.
	filtered_genes = [i[0] for i in orf_utr_list]
	count_filtered_df = count_df.loc[["mirna", "rep", "batch"] + filtered_genes]
	# Make the orf length / utr dataframe:
	orf_utr_df = pd.DataFrame(orf_utr_list, columns=["gene", "orf_length", "utr"])
	orf_utr_df.index = orf_utr_df.iloc[:,0]
	orf_utr_df.index.name = None
	orf_utr_df = orf_utr_df.iloc[:,1:]
	# Make the rpkm dataframe, and convert this to a tpm dataframe:
	rpkm_df = pd.DataFrame(rpkm_list, columns=col_strings)
	col_sum_rpkm_df = rpkm_df.sum(0)/1e6
	tpm_df = rpkm_df.div(col_sum_rpkm_df, axis="columns")
	tpm_df.index = filtered_genes

	# Make the list that checks each UTR site for each miRNA seed:
	mirna_sites_list = []
	for row in orf_utr_df.itertuples():
		utr = row.utr
		sites = [int(seed in utr) for seed in mirna_seeds]
		mirna_sites_list.append(sites)
	# Convert the list into a dataframe:
	mirna_sites_df = pd.DataFrame(mirna_sites_list, columns = col_strings)
	mirna_sites_df.index = filtered_genes
	# print(count_df.iloc[:10,:5])
	# print(rpkm_df.iloc[:10, :5])
	# print(col_sum_rpkm_df)
	# print(tpm_df.iloc[:10,:5])
	# print(tpm_df.sum(0))
	# print(mirna_sites_df.iloc[:10, :5])

	mirna_sites_path = "RepressionData/Lin-Shi_transfection_data/utr_sites_filtered.txt"
	tpm_path = "RepressionData/Lin-Shi_transfection_data/tpm_filtered.txt"
	mirna_sites_df.to_csv(mirna_sites_path, sep="\t", header=True)
	tpm_df.to_csv(tpm_path, sep="\t", header=True)




################################################################################

if __name__ == "__main__":
    main()



