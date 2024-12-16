import subprocess

raw_files_path = "ThreePrimeTargetingPaper/GEO_transfer/raw_files.txt"
processed_files_path = "ThreePrimeTargetingPaper/GEO_transfer/processed_files.txt"
processed_mismatch_files_path = "ThreePrimeTargetingPaper/GEO_transfer/processed_mismatch_files.txt"


def get_checksum(path):
    shell_call = "md5sum %s" %(path)
    p = subprocess.Popen([shell_call], shell=True, stdout = subprocess.PIPE)
    # Parse the output 
    output = p.communicate()[0].decode('utf-8').split("\n")[0]
    return(output)


with open(raw_files_path) as file:
	line = file.readline()
	while line:
		filename = line.strip()
		full_path = "/lab/solexa_bartel/mcgeary/GEO_McGearyBisaria2021/raw_files/" + filename
		full_path = '"' + full_path + '"'
		checksum = get_checksum(full_path)
		print(checksum)
		line = file.readline()

# with open(processed_files_path) as file:
# 	line = file.readline()
# 	while line:
# 		filename = line.strip()
# 		full_path = "/lab/solexa_bartel/mcgeary/GEO_McGearyBisaria2021/processed_files/" + filename
# 		full_path = '"' + full_path + '"'
# 		checksum = get_checksum(full_path)
# 		print(checksum)
# 		line = file.readline()

# with open(processed_mismatch_files_path) as file:
# 	line = file.readline()
# 	while line:
# 		filename = line.strip()
# 		full_path = "/lab/solexa_bartel/mcgeary/GEO_McGearyBisaria2021/processed_files/" + filename
# 		full_path = '"' + full_path + '"'
# 		checksum = get_checksum(full_path)
# 		print(checksum)
# 		line = file.readline()


