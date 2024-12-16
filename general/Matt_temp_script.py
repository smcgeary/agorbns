import gzip

print("hello world")

path = "/lab/solexa_bartel/magetz/Solexa/SmallRNASequencing_12-17/DPB220/"

file = "DPB220_22-23nt_5p3p-trim_q30.txt.gz"

full_path = path + file

print(full_path)

j = 0
out_list = []
with gzip.open(full_path, "rb") as file_in:
	line = file_in.readline()
	while line:
		if j == 1:
			out_list += [line.strip()]
		# i += 1
		j += 1
		if j == 4:
			j = 0
		line = file_in.readline()
		

print(out_list)
path_out = path + "DPB220_22-23nt_5p3p-trim_q30_readsonly.txt"
print(path_out)
with open(path_out, "wb") as file_out:
	file_out.write("\n".join(out_list) + "\n")



