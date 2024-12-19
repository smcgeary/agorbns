# for i in ls /lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/equilibrium/kmers/*buffer*
# do
# 	echo $i
# 	FILE="$(basename $i)"
# 	FILEPATH="$(echo $i | cut -d"/" -f 1-7)/kmers_cutoff_final/"
# 	NEWPATH="$FILEPATH$FILE"
# 	echo $NEWPATH
# 	cp $i $NEWPATH
# done
PATH="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-155/equilibrium/kmers/"
for i in I 40 4 0
do
	for j in 8 10
	do
		echo $i
		echo $j
		NEWPATH="$PATH$i""_""5""_""$j.txt"
		echo $NEWPATH
	done
done

# for i in ls /lab/solexa_bartel/mcgeary/AgoRBNS/miR-155/equilibrium/kmers/*buffer*
# do
# 	echo $i
# 	FILE="$(basename $i)"
# 	FILEPATH="$(echo $i | cut -d"/" -f 1-7)/kmers_cutoff_final/"
# 	NEWPATH="$FILEPATH$FILE"
# 	echo $NEWPATH
# 	cp $i $NEWPATH
# done

