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
mmu_UTR_dict = dict()
with open("RepressionData/UTRs_mm10.txt") as file_in:
	line = file_in.readline()
	while line:
		line_strip = line.strip().split("\t")
		mmu_UTR_dict[line_strip[0]] = line_strip[1]
		line = file_in.readline()
hsa_UTR_dict = dict()
with open("/lab/bartel4_ata/kathyl/RNA_Seq/analysis/data/final/utr3.txt") as file_in:
	line = file_in.readline()
	while line:
		line_strip = line.strip().split("\t")
		hsa_UTR_dict[line_strip[0]] = line_strip[1]
		line = file_in.readline()

# #########MOUSE##################################################################
# # miR-155_______________________________________________________________________
# # 3T3 contact inhibited
# out = []
# missing_mmu = 0
# with open("RepressionData/Eichhorn_data_tables/GSE60426_miR-155_contact_inhibited_3T3_changes.txt", "rU") as file_in:
# 	for i in range(4):
# 		line = file_in.readline()
# 		print("pre line")
# 		print(line)
# 	while line:
# 		line_split = line.strip().split("\t")
# 		gene, rna_4, rna_8, rna_12, rpf_4, rpf_8, rpf_12, te_4, te_8, te_12 = line_split[:10]
# 		line = file_in.readline()
# 		if gene in mmu_UTR_dict:
# 			out.append([gene, rna_12, mmu_UTR_dict[gene]])
# 		else:
# 			missing_mmu += 1
# path_out = "RepressionData/Eichhorn_data_tables/processed/miR-155_mmu_3T3_contactinhibited_12h.txt"
# with open(path_out, "w") as file_out:
# 	file_out.write("\n".join(["\t".join(i) for i in out]))

# # 3T3 actively dividing
# out = []
# missing_mmu = 0
# with open("RepressionData/Eichhorn_data_tables/GSE60426_miR-155_actively_dividing_3T3_0h.txt", "rU") as file_in_1:
# 	with open("RepressionData/Eichhorn_data_tables/GSE60426_miR-155_actively_dividing_3T3_12h.txt", "rU") as file_in_2:
# 		for i in range(4):
# 			line1 = file_in_1.readline()
# 			line2 = file_in_2.readline()
# 			print("pre line")
# 			print(line1)
# 			print(line2)
# 		while line1:
# 			line1_split = line1.strip().split("\t")
# 			line2_split = line2.strip().split("\t")
# 			gene1, name1, rna_0,  rpf_0  = line1_split[:4]
# 			gene2, name2, rna_12, rpf_12 = line2_split[:4]
# 			print([rna_0, rna_12])
# 			line1 = file_in_1.readline()
# 			line2 = file_in_2.readline()
# 			if gene in mmu_UTR_dict:
# 				out.append([gene, rna_12, mmu_UTR_dict[gene]])
# 			else:
# 				missing_mmu += 1
# path_out = "RepressionData/Eichhorn_data_tables/processed/miR-155_mmu_3T3_actively_dividing_12h.txt"
# with open(path_out, "w") as file_out:
# 	file_out.write("\n".join(["\t".join(i) for i in out]))


# #########HUMAN##################################################################
# out = []
# missing_hsa_HeLa_miR_155 = 0
# with open("RepressionData/Eichhorn_data_tables/GSE60426_miR-155_HeLa_changes.txt", "rU") as file_in:
# 	for i in range(4):
# 		line = file_in.readline()
# 		print("pre line")
# 		print(line)
# 	for i in range(10):
# 		line_split = line.strip().split("\t")
# 		print(line_split)
# 		gene, rna, rpf, te, group = line_split[:5]
# 		line = file_in.readline()
# 		if gene in hsa_UTR_dict:
# 			out.append([gene, rna, hsa_UTR_dict[gene]])
# 		else:
# 			missing_hsa_HeLa_miR_155 += 1
# path_out = "RepressionData/Eichhorn_data_tables/processed/miR-155_hsa_HeLa.txt"
# with open(path_out, "w") as file_out:
# 	file_out.write("\n".join(["\t".join(i) for i in out]))

# out = []
# missing_hsa_HeLa_miR_1 = 0
# with open("RepressionData/Eichhorn_data_tables/GSE60426_miR-1_HeLa_changes.txt", "rU") as file_in:
# 	for i in range(4):
# 		line = file_in.readline()
# 		print("pre line")
# 		print(line)
# 	for i in range(10):
# 		line_split = line.strip().split("\t")
# 		print(line_split)
# 		gene, rna, rpf, te, group = line_split[:5]
# 		line = file_in.readline()
# 		if gene in hsa_UTR_dict:
# 			out.append([gene, rna, hsa_UTR_dict[gene]])
# 		else:
# 			missing_hsa_HeLa_miR_1 += 1
# path_out = "RepressionData/Eichhorn_data_tables/processed/miR-1_hsa_HeLa.txt"
# with open(path_out, "w") as file_out:
# 	file_out.write("\n".join(["\t".join(i) for i in out]))


# print(missing_mmu)
# print(missing_hsa_HeLa_miR_155)
# print(missing_hsa_HeLa_miR_1)

mir_batch_dict = {"mir144_rep1"  : 1,
				  "mir199a_rep1"  : 1,
				  "mir204_rep1"  : 1,

				  "mir137_rep1"  : 2,
				  "mir143_rep1"  : 2,
				  "mir144_rep2"  : 2,
				  "mir155_rep1"  : 2,
				  "mir199a_rep2" : 2,
				  "mir204_rep2"  : 2,
				  "mir205_rep1"  : 2,
				  "mir223_rep1"  : 2,

				  "mir7_rep1"    : 3,
				  "mir137_rep2"  : 3,
				  "mir139_rep1"  : 3,
				  "mir153_rep1"  : 3,
				  "mir182_rep1"  : 3,
				  "mir205_rep2"  : 3,
				  "mir216b_rep1" : 3,

				  "let7_rep2"    : 6,
				  "lsy6_rep2"    : 6,
				  "mir1_rep2"    : 6,
				  "mir124_rep2"  : 6,
				  "mir143_rep2"  : 6,
				  "mir153_rep2"  : 6,
				  "mir216b_rep2" : 6,
				  "mir223_rep2"  : 6,

				  "let7_rep1"   : 7,
				  "lsy6_rep1"   : 7,
				  "mir1_rep1"   : 7,
				  "mir7_rep2"   : 7,
				  "mir124_rep1" : 7,
				  "mir139_rep2" : 7,
				  "mir155_rep2" : 7,
				  "mir182_rep2" : 7
}

print(mir_batch_dict)

rep_dir = "/lab/bartel4_ata/kathyl/RNA_Seq/transfections/compiled/STAR_vikram"
rep_paths = ls(rep_dir)
print(rep_paths)

def MakeDataColumn(path):
	# Computes a single column of count data from the HeLa transfection
	# datasets.
	mirna_rep = path.split("_compiled.txt")[0]
	batch = mir_batch_dict[mirna_rep]
	mirna, rep = mirna_rep.split("_rep")
	header = pd.DataFrame([mirna, rep, batch], index=["mirna", "rep", "batch"])
	counts = pd.read_csv("%s/%s" %(rep_dir, path), header=None, index_col=0,
	                     sep="\t").iloc[:4312]
	counts = counts.rename_axis(None)
	counts.columns = [0]
	out = pd.concat([header, counts], axis=0)
	return(out)


def main():
	# Make the First Column:
	rep_df = MakeDataColumn(rep_paths[0])
	print(rep_df.iloc[:10])

	for path in rep_paths[1:]:
		rep_df_col = MakeDataColumn(path)
		rep_df = pd.concat([rep_df, rep_df_col], axis=1)
		print(rep_df.iloc[:5,:])

	print(rep_df.iloc[:5, :])
	# rep_df.to_csv("RepressionData/Lin-Shi_transfection_data/counts.txt", sep="\t", header=False)

################################################################################

if __name__ == "__main__":
    main()



