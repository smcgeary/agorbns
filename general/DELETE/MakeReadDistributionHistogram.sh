mirna=$1
exp=$2
cond=$3

DIR="/lab/solexa_bartel/mcgeary/AgoRBNS/"
indir=$DIR$mirna/$exp/reads
inpath=$indir/$cond.txt
outdir=$DIR$mirna/$exp/read_count_distributions
outpath=$outdir/$cond.txt
echo $outdir
if [ ! -d "$outdir" ]; then
	echo "out dir doesn't exist"
	mkdir $outdir
else
	echo "out dir exists"
fi
if [ ! -d "$indir" ]; then
	echo "in dir doesn't exist"
else
	echo "in dir exists"
fi

echo $outpath
echo $inpath
sort $inpath | uniq -c | sort -n | awk -F" " '{print $1}' | uniq -c | sort -n > $outpath
