from collections import defaultdict
import pandas as pd

path = "/lab/solexa_bartel/magetz/Solexa/SmallRNASequencing_12-17/DPB220/"
file = "DPB220_22-23nt_5p3p-trim_q30_readsonly.txt"

full_path = path + file

unique_seqs = defaultdict(int)

print(unique_seqs)

with open(full_path) as file_in:
	sequence = file_in.readline().strip()
	while sequence:
		unique_seqs[sequence] += 1
		sequence = file_in.readline().strip()



# print(unique_seqs)


unique_seqs_df = pd.DataFrame.from_dict(unique_seqs, orient="index")
unique_seqs_df.columns = ["counts"]

unique_seqs_df.sort_values("counts", ascending=False, axis=0, inplace=True)

print(unique_seqs_df.iloc[:10,])




