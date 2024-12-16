import subprocess

geo_path = "ThreePrimeTargetingPaper/GEO_transfer/GEO_table_reporter_assay.txt"


def get_checksum(path):
    shell_call = "md5sum %s" %(path)
    p = subprocess.Popen([shell_call], shell=True, stdout = subprocess.PIPE)
    # Parse the output 
    output = p.communicate()[0].decode('utf-8').split("\n")[0]
    return(output)

# Get the checksums for the raw file:
with open(geo_path) as file:
	line = file.readline()
	line = file.readline()
	while line:
		filename = line.split("\t")[3].strip()
		full_path = "/lab/solexa_bartel/mcgeary/GEO_McGearyBisaria2021_3/raw_files/" + filename
		full_path = '"' + full_path + '"'
		checksum = get_checksum(full_path)
		print(checksum)
		line = file.readline()

# Get the checksums for the processed files.
full_path = "/lab/solexa_bartel/mcgeary/GEO_McGearyBisaria2021_3/processed_files/mpra_cpm.txt"
full_path = '"' + full_path + '"'
checksum = get_checksum(full_path)
print(checksum)


