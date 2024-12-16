

data_seqs = open("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/twist_reporter_assay_3p_tp/reads/duplex_series,1_counted_ranked.txt")

data_lines = data_seqs.readlines()

ref_seqs = open("ThreePrimePaperLet7ReporterAssay/libraries_210726/oligo_df_210726_library0.txt")

ref_lines = ref_seqs.readlines()

total_variants = len(ref_lines) - 1

# seq_dict_100 = dict()
seq_dict = dict()


for ref_line in ref_lines[1:]:
	entries = ref_line.split()
	seq_use = entries[-1]
	seq_use = seq_use[17:].split("TACCAATGCCCTGGCTC")[0]
	seq_dict[seq_use] = entries[:-1]


variants_encountered = []
ambiguous = []

for i, ref_line in enumerate(ref_lines[1:]):
	entries = ref_line.split()
	seq_use = entries[-1]
	expected_read = seq_use[17:][:100]
	total_keys = 0
	keys = []
	for key in seq_dict.keys():
		if expected_read in key:
			keys += [key]
			variant_num = seq_dict[key][0]
			if int(variant_num) not in variants_encountered:
				variants_encountered += [int(variant_num)]
	if len(keys) > 1:
		print("%s : %s" %(i, expected_read))
		print("total is %s" %len(keys))
		for key in keys:
			print(seq_dict[key])
		ambiguous += [[keys]]


print(len(variants_encountered))
print(total_variants)
print(len(ambiguous))

variants_encountered = []
ambiguous = []

# for i in range(1):
for i in range(len(data_lines)):
	if i%10000 == 0:
		print(i)
		print(len(variants_encountered))
		print(len(ambiguous))
	count, seq_i = data_lines[i].split()
	total_keys = 0
	keys = []
	for key in seq_dict.keys():
		if seq_i in key:
			keys += [key]
			variant_num = seq_dict[key][0]
			if int(variant_num) not in variants_encountered:
				variants_encountered += [int(variant_num)]
	if len(keys) > 1:
		print("%s : %s | %s" %(i, count, seq_i))
		print("total is %s" %len(keys))
		for key in keys:
			print(seq_dict[key])
			ambiguous += [[keys]]

# print(len(variants_encountered))
# print(total_variants)


	