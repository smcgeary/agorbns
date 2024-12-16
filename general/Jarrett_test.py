# Test script for Jarrett.
# Going to make ECDF

# Path to file.
file_path = "/lab/solexa_bartel/jsmith/Experiments/E7/E7.1/js_i3_alignments/js_bed6_alignments/E7.1_js_i3/ATCACG-s_7_1_sequence_unique_intersect_bed.txt"

# Opening the file, creates the file object.
file = open(file_path)

# Loop in which we iterate over each line.

trans_dict = dict()

trans_dict["all"] = []

line = file.readline()

while line:
	line_split = line.split("\t")
	chrom = line_split[0]
	if chrom == "chrM":
		read_l = int(line_split[6])
		orf_l = int(line_split[13])
		orf_r = int(line_split[14])
		trans = line_split[15]
		if trans not in trans_dict.keys():
			trans_dict[trans] = []
		read_orf_l = read_l - orf_l
		orf_len = orf_r - orf_l + 1
		read_orf_l_norm = read_orf_l/orf_len
		trans_dict[trans].append(read_orf_l_norm)
		trans_dict["all"].append(read_orf_l_norm)
	line = file.readline()

# print(trans_dict)
# print(trans_dict.keys())
# print(len(trans_dict.keys()))

# for key, value in trans_dict.items():
# 	print("%s: %s" %(key, len(value)))

for key, value in trans_dict.items():
# 	print("____________________________________")
# 	print(value[:10])
	trans_dict[key] = sorted(value)
# 	print(trans_dict[key][:10])


for key in trans_dict.keys():
	total_reads = len(trans_dict[key])
	counts = [(i + 1)/total_reads for i in range(total_reads)]
	path = "temp_Jarrett_%s.txt" %key
	print(path)
	out_file = open(path, "w")
	for i in zip(trans_dict[key], counts):
		out_str = "%s\t%s\n" %i
		out_file.write(out_str)
	out_file.close()
	