# job=""; \
# echo $job; \
# for CON in I 40 12.6 4 1.26 0.4 0 0_combined; do \
# 	job+="python /lab/bartel1_ata/mcgeary/computation/AgoRBNS/"; \
# 	job+="AssignSiteTypes/AssignSites.py miR-1 equilibrium $CON 5 paper\n"; \
# done; \
# echo "$job"; \
# echo $job; \


job=""; \
for MIRNA in miR-7-23nt miR-7-24nt miR-7-25nt; do \
	for LEN_K in $(seq 8 11); do \
		if [ "$MIRNA" == "miR-7-23nt" ]; then \
			COND='I I_combined 40 12.6 1.26 0.4 0'; \
		else \
			COND='I I_combined 40 12.6 4 1.26 0.4 0'; \
		fi; \
		for CON in $COND; do \
			job="python  $MIRNA "; \
			job=$job" $CON 5 $LEN_K\n"; \
			echo $job; \
		done; \
	done; \
done; \
# echo $job; \
# python $(HOME)general/SubmitJob.py "$$job"
