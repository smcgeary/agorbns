j=0; \
for i in /lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/temp_pl_fold_2018_/*; do \
	if [[ $i = *"_temp"* ]]; then \
		j=$(( $j + 1 )); \
		rm -f $i; \
	fi
done
echo $j \
# echo $job; \
# python $(HOME)general/SubmitJob.py "$$job"
