# SHELL := /bin/bash

MIRNAe = let-7a miR-1 miR-155 miR-124 lsy-6 miR-7-23nt
MIRNAk = let-7a miR-1 miR-124 lsy-6
MIRNAeMiR7 = miR-7-22nt miR-7-23nt miR-7-24nt miR-7-25nt
MIRNAeMiR7_2 = miR-7-23nt miR-7-24nt miR-7-25nt
MIRNAeMiR7_3 = miR-7-24nt
MIRNAeThrP = let-7a-21nt miR-1 miR-155 let-7a_plus1 let-7a_minus1 let-7a_miR-155 miR-155_let-7a
MIRNAbc_HeLa_seventeen = miR-137 miR-205 miR-155 miR-223 miR-144 miR-143 miR-153 miR-216b miR-199a miR-204 miR-139 miR-182 miR-7 miR-1 miR-124 lsy-6 let-7a
MIRNAbc_HeLa_sixteen = miR-137 miR-205 miR-155 miR-223 miR-144 miR-143 miR-153 miR-216b miR-199a miR-204 miR-139 miR-182 miR-7 miR-1 miR-124 lsy-6
MIRNAbc_HeLa_six = let-7a miR-155 miR-7 miR-1 miR-124 lsy-6
MIRNAbc_HeLa_five = miR-155 miR-7 miR-1 miR-124 lsy-6
MIRNAbc_HeLa_four = miR-155 miR-7 miR-1 miR-124

DIR_pp = PreProcessReads/
DIR_plfold = AnalyzeStructure/
SCR_mp = MakeMiRNAReadFile.py
SCR_mp = MakeMiRNAReadFile.sh
SCR_pp = MakeReadFile.py
SCR_pp = MakeReadFile.sh
# SCR_pp_burge = MakeReadFileBurge.py

DIR_as = AssignSiteTypes/
SCR_as = AssignSites.py # Used for AGO-RBNS paper
SCR_as = AssignSites.sh # Used for AGO-RBNS paper

SCR_abs = AssignBipartiteSitesProgrammedLib.py # Used for ThrP paper
SCR_abms = AssignBipartiteMismatchSitesProgrammedLib.py # Used for ThrP paper
SCR_abmsr = AssignBipartiteMismatchSitesRandLib.py # Used for ThrP paper
# SCR_absn = AssignBipartiteSitesProgrammedLibNew.py # Used for ThrP paper
## USED TO GET THE KD VALUES WITH 1 to 1 correspondence with the sites in #
## the let-7a reporter library. ###########################################
SCR_absprl = AssignBipartiteSitesProgrammedLib_reporter_temp.py ###########
SCR_absr = AssignBipartiteSitesRandLib.py # Used for ThrP paper
SCR_akbo = AssessProgrammedKmers.py # Used for ThrP paper
SCR_ak = AssignKmers.py
# SCR_akr = AssignKmersNew.py
SCR_ack = AssignCompetitorKmers.py
SCR_apk = AssignPositionalKmers.py
SCR_akpl = AssignKmersProgrammedLib.py # Used for ThrP paper
SCR_af = AssignFlanks.py # Used for AGO-RBNS paper
SCR_af = AssignFlanks.sh # Used for AGO-RBNS paper


DIR_kd = SolveForKds/
SCR_sc = MakeSiteCountTable.py
SCR_msc = MakeMultiSiteCountTable.py
SCR_fc = MakeFlankCountTable.py
SCR_skd = FitSiteKds.R
SCR_ppskd = FitProgrammedPositionKds.R
SCR_pmmkd = FitProgrammedMismatchKds.R
SCR_12merkd = Fit12merKds.R
SCR_fkd = FitFlankKds.R
SCR_kdrepfile = MakeKdStringFile.py

DIR_plfold = AnalyzeStructure/
SCR_plfold = MakeSitePlfoldFiles_PAPER.py

COND_mirna = I 40 12.6 4 1.26 0.4 0



#### NEW DIRECTORIES RELATED TO FITTING BIOCHEM MODEL FOR THREEPRIME SCORE


HELA_TPM_KL_DIR = /lab/solexa_bartel/klin/miRNA_models_data/transfections/hela/$\
	processed/log_tpm_normed.txt
HELA_TPM_SM_DIR = /lab/solexa_bartel/mcgeary/transfections/HeLa/count_tables/$\
	logtpm_batchnormalized.txt

TPM_KL_DIR = /lab/solexa_bartel/klin/miRNA_models_data/transfections/
TPM_SM_DIR = /lab/solexa_bartel/mcgeary/transfections/

# DIR_rep = Repression/# Used for in vivo part of Threep paper.
# SCR_dts = DevelopThreePScore.py# Use for in vivo part of Threep paper.
# SCR_fts = kathy_scripts/target_scan/FitTargetScan.py# Use for in vivo part of Threep paper.


KL_FEAT_MEASURED_DIR = /lab/solexa_bartel/klin/miRNA_models_data_old/$\
	model_inputs/biochem/measured_kds/
KL_FEAT_PREDICTED_DIR = /lab/solexa_bartel/klin/miRNA_models_data_old/$\
	model_inputs/biochem/predicted_kds/feat1_allgenes_lr003_nodropout_batch50_$\
	rebalancekds2k_noaugment_repweight095_mask_w3_netpred/



# the $(filter $(exp), X Y) allows one to test if $(exp) is equal to string X or
# string Y. 
ifeq ($(exp), $(filter $(exp), equilibrium_tp equilibrium_tp_met))
	ifeq ($(mirna), miR-124)
		COND = I I_combined 40 12.6 4,1 4,2 1.26,1 1.26,2 0.4 0.126 0.04 0
		COND_P = I 40 12.6 0.4 0
		COND_P_DUP = 4 1.26
	else ifeq ($(mirna), miR-7-24nt)
		COND = I,1 I,2 40,1 40,2 40,3 40,4 12.6,1 12.6,2 4,1 4,2 4,3 4,4 1.26,1\
			1.26,2 0.4,1 0.4,2 0,1 0,2
		COND = I 40 12.6 4 1.26 0.4 0
		COND_P = I 12.6 1.26 0.4 0
		COND_P_DUP = 40 4
	else
		COND = I 40 12.6 4 1.26 0.4 0
		COND_P = $(COND)
	endif
else ifeq ($(exp), equilibrium_2_tp)
	COND = I,1 I,2 40,1 40,2 12.6,1 12.6,2 4,1 4,2 1.26,1 1.26,2 0.4,1 0.4,2\
		0,1 0,2
	COND = I I_combined 40 12.6 4 1.26 0.4 0
	COND_P =
else ifeq ($(exp), $(filter $(exp), equil_c2_nb equil_c2_alt_nb))
	COND = I 40 12.6 12.6_2 4 1.26 0.4 0
	COND_P = $(COND)
else ifeq ($(exp), equil_c_nb)
		COND = I 40 12.6 4 1.26 0.4 0
		COND_P = $(COND)
else ifeq ($(exp), equil_s_nb)
	COND = I 40 12.6 12.6_2 4 4_2 1.26 0.4 0
	COND_P = $(COND)
else ifeq ($(exp), $(filter $(exp), equil_sc_nb equil_sc_alt_nb))
	COND = I 40 12.6 4 1.26 0.4 0
	COND_P = $(COND)
else ifeq ($(exp), equilibrium3_nb)
	COND = I 40 12.6 4 1.26 0.4 0
	COND_P = $(COND)
else ifeq ($(exp), equilibrium2_nb)
		COND = I I_combined 40 12.6 1.26 0.4 0
		COND_P = I 40 12.6 1.26 0.4 0
else ifeq ($(exp), equil_pilot)
	COND = I L100A10
else
	COND = I I_combined 40 12.6 4 1.26 0.4 0
	COND_P = I 40 12.6 4 1.26 0.4 0
endif



ifdef buffer
	BUFFER3P = ' -buffer3p'
else
	BUFFER3P = ''
endif

ifdef unique
	UNIQUE = ' -uniq'
else
	UNIQUE = ''
endif

ifdef n_ex
	N_EX = ' -n_ex '$(n_ex)
	N_EX_ALL = ' n_ex='$(n_ex)
else
	N_EX = ''
	N_EX_ALL = ''
endif

# ifdef n_lag
# 	N_LAG = ' -n_lag '$(n_lag)
# 	N_LAG_ALL = ' n_lag='$(n_lag)
# else
# 	N_LAG = ''
# 	N_LAG_ALL = ''
# endif

ifdef n_s
	N_S = $(n_s)
else
	N_S = 0
endif

ifdef addcomp
	ADDCOMP = ' -addcomp '$(addcomp)
else
	ADDCOMP = ''
endif

ifdef wob
	WOB = ' -wob'
else
	WOB = ''
endif

# This conditional definition allows for an alternative competitor oligo to be
# queried (specifically used to figure out which let-7a is used in the second
# equilibrium_mmseed_nb experiment).
ifdef altmircomp
	ALTMIRCOMP = ' -alt_mirna_comp '$(altmircomp)
else
	ALTMIRCOMP = ''
endif


LEN_MAX = 11
LEN_MIN = 4





PreprocessAgoPurity :
	@(job=""; \
	for MIRNA in miR-1 miR-155; do \
		job="python $(DIR_pp)$(SCR_mp) $$MIRNA AGO_purity S1007_P "; \
		job=$$job" -test"; \
		echo $$job; \
		$$job;\
	done)




# takes functions $(mirna), $(exp), and $(test)

## Logic behind this is that is that some experiments have only reps 1 and 2.
## One experiment, that of miR-7-24nt/equilibrium_tp, has reps for all samples,
## and two more reps for samples 40 and 4, which are reps 3 and 4. This deals
## thta using the COND_P (the non-duplicated samples) and the COND_P_DUP (the 
## duplicated samples) using a different duplication strategy for that
## experiment.
PreprocessData :
	@(job_prefix="sbatch $(DIR_pp)$(SCR_pp) $(mirna) $(exp) "; \
	echo $(COND_P) ;\
	echo $(COND_P_DUP) ;\
	if [ "$(mirna)" = "miR-7-24nt" ] && [ "$(exp)" = "equilibrium_tp" ]; then \
		DUP_REPS="1 2 3 4"; \
		for CON in $(COND_P); do \
			for REP in 1 2; do \
				job=$$job_prefix"$$CON -rep $$REP$(test) -jobs 19"; \
				echo $$job; \
				# $$job; \
			done ;\
		done; \
	else \
		DUP_REPS="1 2"; \
		for CON in $(COND_P); do \
			job=$$job_prefix"$$CON$(test) -jobs 19"; \
			echo $$job; \
			# $$job; \
		done; \
	fi; \
	for CON in $(COND_P_DUP); do \
		for REP in $$DUP_REPS; do \
			job=$$job_prefix"$$CON -rep $$REP$(test) -jobs 19"; \
			echo $$job; \
			# $$job; \
		done; \
	done)



PreprocessAllEquilibrium :
	make mirna=miR-1 exp=equilibrium PreprocessData
	make mirna=let-7a exp=equilibrium PreprocessData
	make mirna=miR-155 exp=equilibrium PreprocessData
	make mirna=miR-124 exp=equilibrium PreprocessData
	make mirna=lsy-6 exp=equilibrium PreprocessData
	make mirna=miR-7-23nt exp=equilibrium2_nb PreprocessData
	make mirna=mir-1 exp=equil_pilot PreprocessData
	make mirna=miR-1 exp=equilibrium_tp PreprocessData
	make mirna=miR-124 exp=equilibrium_2_tp PreprocessData
	make mirna=miR-7-24nt exp=equilibrium_tp PreprocessData	
	job="sbatch $(DIR_pp)$(SCR_pp) miR-1 kin_pilot I_TGT -jobs 19"; \
	echo $$job; \
	$$job; \
	job="sbatch $(DIR_pp)$(SCR_pp) miR-7-24nt equilibrium3_nb I -jobs 19"; \
	echo $$job; \
	$$job


AssignSites :
	@(for CON in $(COND); do \
		job="sbatch $(DIR_as)$(SCR_as) $(mirna) $(exp) $$CON "; \
		job=$$job"$(n_constant) $(sitelist)"; \
		job=$$job$(BUFFER3P)" -jobs 19"; \
		echo $$job; \
		$$job; \
	done)

AssignFlanks :
	@(for CON in $(COND); do \
		job="sbatch $(DIR_as)$(SCR_af) $(mirna) $(exp) $$CON "; \
		job=$$job"$(n_constant) $(sitelist)"$(BUFFER3P)" -jobs 19"; \
		echo $$job; \
		$$job; \
	done)





AssignSitesJustCombinedInput :
# 	@(job="python $(HOME)$(DIR_as)$(SCR_as) $(mirna) $(exp) I_combined "; \
# 	  job=$$job"$(n_constant) $(sitelist)"$(BUFFER3P)$(UNIQUE)"\n"; \
# 	  echo $$job; \
# 	  bsub -q 18 -n 20 $$job)
	@(job="sbatch $(DIR_as)$(SCR_as) $(mirna) $(exp) I_combined $(n_constant) "; \
	  job=$$job"$(sitelist)"$(BUFFER3P)$(UNIQUE)" -jobs 19"; \
	  echo $$job; \
	  $$job)


AssignSitesAllEquilibrium :
	make mirna=miR-1 exp=equilibrium n_constant=$(n_constant) sitelist=$(sitelist) buffer=1 AssignSites
	make mirna=let-7a exp=equilibrium n_constant=$(n_constant) sitelist=$(sitelist) AssignSites
	make mirna=miR-155 exp=equilibrium n_constant=$(n_constant) sitelist=$(sitelist) AssignSites
	make mirna=miR-124 exp=equilibrium n_constant=$(n_constant) sitelist=$(sitelist) AssignSites
	make mirna=lsy-6 exp=equilibrium n_constant=$(n_constant) sitelist=$(sitelist) AssignSites
	make mirna=miR-7-23nt exp=equilibrium2_nb n_constant=$(n_constant) sitelist=$(sitelist) AssignSites

AssignSitesEquilibriumAllSiteLists :
	# make mirna=miR-1 exp=equilibrium n_constant=5 sitelist=canonical buffer=1 AssignSites
	# make n_constant=5 sitelist=resubmissionfinal AssignSitesAllEquilibrium
	make n_constant=5 sitelist=centered11 AssignSitesAllEquilibrium



# PreprocessEquilibrium :
# 	@(for CON in $(COND_mirna); \
# 		do { echo python $(DIR_pp)$(SCR_pp) $(mirna) $(exp) $$CON $(nb) $(test); \
# 		bsub -n 20 python $(DIR_pp)$(SCR_pp) $(mirna) $(exp) $$CON $(nb) $(test);} \
# 	done)

# PreprocessDataOld :
# 	@(job_prefix="python $(DIR_pp)$(SCR_pp) $(mirna) $(exp) "; \
# 	#job_prefix="bsub -q 18 -n 20 -R span[hosts=1] python $(DIR_pp)$(SCR_pp) $(mirna) $(exp) "; \
# 	if [ "$(mirna)" = "miR-7-24nt" ] && [ "$(exp)" = "equilibrium_tp" ]; then \
# 		DUP_REPS="1 2 3 4"; \
# 		for CON in $(COND_P); do \
# 			for REP in 1 2; do \
# 				job=$$job_prefix"$$CON -rep $$REP$(test) -jobs 19"; \
# 				#echo $$job; \
# 				$$job; \
# 			done ;\
# 		done; \
# 	else \
# 		DUP_REPS="1 2"; \
# 		for CON in $(COND_P); do \
# 			job=$$job_prefix"$$CON$(test) -jobs 19"; \
# 			#echo $$job; \
# 			$$job; \
# 		done; \
# 	fi; \
# 	for CON in $(COND_P_DUP); do \
# 		for REP in $$DUP_REPS; do \
# 			job=$$job_prefix"$$CON -rep $$REP$(test) -jobs 19"; \
# 			#echo $$job; \
# 			$$job; \
# 		done; \
# 	done)


PreprocessHighThroughputv2 :
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_v2 duplex -rep 1
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_v2 duplex -rep 1
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py miR-155 twist_reporter_assay_v2 duplex -rep 1
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py miR-124 twist_reporter_assay_v2 duplex -rep 1
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py lsy-6 twist_reporter_assay_v2 duplex -rep 1
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py miR-7 twist_reporter_assay_v2 duplex -rep 1
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_v2 duplex -rep 2
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_v2 duplex -rep 2
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py miR-155 twist_reporter_assay_v2 duplex -rep 2
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py miR-124 twist_reporter_assay_v2 duplex -rep 2
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py lsy-6 twist_reporter_assay_v2 duplex -rep 2
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py miR-7 twist_reporter_assay_v2 duplex -rep 2
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_v2 no_duplex -rep 1
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_v2 no_duplex -rep 1
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py miR-155 twist_reporter_assay_v2 no_duplex -rep 1
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py miR-124 twist_reporter_assay_v2 no_duplex -rep 1
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py lsy-6 twist_reporter_assay_v2 no_duplex -rep 1
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py miR-7 twist_reporter_assay_v2 no_duplex -rep 1
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_v2 no_duplex -rep 2
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_v2 no_duplex -rep 2
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py miR-155 twist_reporter_assay_v2 no_duplex -rep 2
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py miR-124 twist_reporter_assay_v2 no_duplex -rep 2
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py lsy-6 twist_reporter_assay_v2 no_duplex -rep 2
	bsub -q 18 -n 20 python PreProcessReads/MakeReadFile.py miR-7 twist_reporter_assay_v2 no_duplex -rep 2

CountReporterVariantsV2 :
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py let-7a twist_reporter_assay_v2 duplex -rep 1 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py let-7a twist_reporter_assay_v2 duplex -rep 2 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-1 twist_reporter_assay_v2 duplex -rep 1 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-1 twist_reporter_assay_v2 duplex -rep 2 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-155 twist_reporter_assay_v2 duplex -rep 1 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-155 twist_reporter_assay_v2 duplex -rep 2 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-124 twist_reporter_assay_v2 duplex -rep 1 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-124 twist_reporter_assay_v2 duplex -rep 2 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py lsy-6 twist_reporter_assay_v2 duplex -rep 1 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py lsy-6 twist_reporter_assay_v2 duplex -rep 2 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-7 twist_reporter_assay_v2 duplex -rep 1 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-7 twist_reporter_assay_v2 duplex -rep 2 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py let-7a twist_reporter_assay_v2 no_duplex -rep 1 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py let-7a twist_reporter_assay_v2 no_duplex -rep 2 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-1 twist_reporter_assay_v2 no_duplex -rep 1 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-1 twist_reporter_assay_v2 no_duplex -rep 2 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-155 twist_reporter_assay_v2 no_duplex -rep 1 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-155 twist_reporter_assay_v2 no_duplex -rep 2 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-124 twist_reporter_assay_v2 no_duplex -rep 1 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-124 twist_reporter_assay_v2 no_duplex -rep 2 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py lsy-6 twist_reporter_assay_v2 no_duplex -rep 1 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py lsy-6 twist_reporter_assay_v2 no_duplex -rep 2 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-7 twist_reporter_assay_v2 no_duplex -rep 1 -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-7 twist_reporter_assay_v2 no_duplex -rep 2 -jobs 19

CountReporterVariantsV2_mm :
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py let-7a twist_reporter_assay_v2 duplex -rep 1 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py let-7a twist_reporter_assay_v2 duplex -rep 2 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-1 twist_reporter_assay_v2 duplex -rep 1 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-1 twist_reporter_assay_v2 duplex -rep 2 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-155 twist_reporter_assay_v2 duplex -rep 1 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-155 twist_reporter_assay_v2 duplex -rep 2 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-124 twist_reporter_assay_v2 duplex -rep 1 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-124 twist_reporter_assay_v2 duplex -rep 2 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py lsy-6 twist_reporter_assay_v2 duplex -rep 1 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py lsy-6 twist_reporter_assay_v2 duplex -rep 2 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-7 twist_reporter_assay_v2 duplex -rep 1 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-7 twist_reporter_assay_v2 duplex -rep 2 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py let-7a twist_reporter_assay_v2 no_duplex -rep 1 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py let-7a twist_reporter_assay_v2 no_duplex -rep 2 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-1 twist_reporter_assay_v2 no_duplex -rep 1 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-1 twist_reporter_assay_v2 no_duplex -rep 2 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-155 twist_reporter_assay_v2 no_duplex -rep 1 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-155 twist_reporter_assay_v2 no_duplex -rep 2 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-124 twist_reporter_assay_v2 no_duplex -rep 1 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-124 twist_reporter_assay_v2 no_duplex -rep 2 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py lsy-6 twist_reporter_assay_v2 no_duplex -rep 1 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py lsy-6 twist_reporter_assay_v2 no_duplex -rep 2 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-7 twist_reporter_assay_v2 no_duplex -rep 1 -mm -jobs 19
	bsub -q 18 -n 20 python ReporterScreen/CountReporterVariants.py miR-7 twist_reporter_assay_v2 no_duplex -rep 2 -mm -jobs 19


# PreprocessEquilibriumPilot :
# 	@(job=""; \
# 	for CON in $(COND); do \
# 		job=$$job"python $(DIR_pp)$(SCR_pp) $(mirna) equil_pilot $$CON "; \
# 		job=$$job"$(nb) $(test)\n"; \
# 	done; \
# 	echo $$job; \
# 	python $(HOME)general/SubmitJob.py "$$job")


# PreprocessEquilibriumNB :
# 	@(for CON in $(COND); \
# 		do { echo python $(HOME)$(DIR_pp)$(SCR_pp) $(mirna) $(exp) $$CON -nb; \
# 		bsub -q bartel -n 20 python $(HOME)$(DIR_pp)$(SCR_pp) $(mirna) $(exp) $$CON -nb;} \
# 	done)





PreprocessAllEquilThrP :
	make mirna=let-7a-21nt exp=equil_c_nb PreprocessData
	#make mirna=let-7a-21nt exp=equil_s_nb PreprocessData
	#make mirna=let-7a-21nt exp=equil_c2_alt_nb PreprocessData
	#make mirna=miR-1 exp=equil_c_alt_nb PreprocessData
	#make mirna=miR-1 exp=equil_sc_alt_nb PreprocessData
	make mirna=let-7a-21nt exp=equil_c2_nb PreprocessData
	make mirna=miR-1 exp=equil_c_nb PreprocessData
	#make mirna=miR-1 exp=equi_sc_nb PreprocessData
	make mirna=let-7a_plus1 exp=equil_c_nb PreprocessData
	make mirna=let-7a_minus1 exp=equil_c_nb PreprocessData
	#make mirna=let-7a-21nt exp=equil_sc_nb PreprocessData
	make mirna=miR-155 exp=equil_sc_nb PreprocessData
	make mirna=let-7a_miR-155 exp=equil_c_nb PreprocessData
	make mirna=miR-155_let-7a exp=equil_c_nb PreprocessData







################### ALL OF THESE FUNCTIONS ARE FOR THE THREEPRIME PAPER

# FOR THREE PRIME PAPER
AssignBipartiteSites :
	@(for CON in $(COND); do \
		job="bsub -q 18 -n 20 -R span[hosts=1] python "; \
# 		job="python "; \
		job=$$job"$(HOME)$(DIR_as)$(SCR_abs) $(mirna) $(exp) $$CON "; \
		job=$$job"$(n_constant) -jobs 19"; \
		echo $$job; \
# 		$$job; \
	done)


AssignBipartiteSitesWithLet7Reporter :
	@(for CON in $(COND); do \
		job="bsub -q 18 -n 20 -R span[hosts=1] python "; \
		job=$$job"$(HOME)$(DIR_as)$(SCR_absprl) $(mirna) $(exp) $$CON "; \
		job=$$job"$(n_constant) -jobs 19"; \
		echo $$job; \
		$$job; \
	done)




AssignBipartiteMismatchSites :
	@(for CON in $(COND); do \
		job="bsub -q 18 -n 20 -R span[hosts=1] python "; \
		job=$$job"$(HOME)$(DIR_as)$(SCR_abms) $(mirna) $(exp) $$CON "; \
		job1=$$job"$(n_constant) $(start_mm) $(stop_mm) -jobs 19"; \
		job2=$$job"$(n_constant) $(start_mm) $(stop_mm) -new -jobs 19"; \
		echo $$job1; \
		echo $$job2; \
		$$job1; \
		$$job2; \
	done)


AssignRandomBipartiteMismatchSites :
	@(for CON in "I_combined" $(COND); do \
		job="bsub -q 18 -n 20 -R span[hosts=1] python "; \
		job=$$job"$(HOME)$(DIR_as)$(SCR_abmsr) $(mirna) $(exp) $$CON "; \
		job=$$job"$(n_constant) $(start_mm) $(stop_mm)"$(BUFFER3P)" -jobs 19"; \
		echo $$job; \
		$$job; \
	done)



AssignAllRandomMismatchSitesOneMirna :
	if [ "$(mirna)" = "miR-155" ]; then \
		LEN_MIR=$$(( 23 )); \
		EXP="equilibrium"; \
	elif [ "$(mirna)" = "miR-7-23nt" ]; then \
		LEN_MIR=$$(( 23 )); \
		EXP="equilibrium2_nb"; \
	else \
		LEN_MIR=22; \
		EXP="equilibrium"; \
	fi; \
	if [ "$(mirna)" = "miR-1" ]; then \
		BUFFER=" buffer=1"; \
	else \
		BUFFER=""; \
	fi; \
	for LEN in $$(seq 4 9); do \
		for START in $$(seq 9 $$(( $$LEN_MIR - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="make mirna=$(mirna) exp=$$EXP "; \
			job=$$job"n_constant=$(n_constant) start_mm=$$START "; \
			job=$$job"stop_mm=$$STOP$$BUFFER AssignRandomBipartiteMismatchSites"; \
			echo $$job; \
			$$job; \
		done; \
	done


FitAllRandomMismatchSitesOneMirna :
	@(if [ "$(mirna)" = "miR-155" ]; then \
		LEN_MIR=$$(( 23 )); \
		EXP="equilibrium"; \
	elif [ "$(mirna)" = "miR-7-23nt" ]; then \
		LEN_MIR=$$(( 23 )); \
		EXP="equilibrium2_nb"; \
		COMB=" -nocombI"; \
	else \
		LEN_MIR=22; \
		EXP="equilibrium"; \
		COMB=""; \
	fi; \
	if [ "$(mirna)" = "miR-1" ]; then \
		BUFFER=" -buffer"; \
		COMB=" -nocombI"; \
	else \
		BUFFER=""; \
	fi; \
	for LEN in $$(seq 4 9); do \
		for START in $$(seq 9 $$(( $$LEN_MIR - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="python SolveForKds/MakeSiteCountTable.py $(mirna) "; \
			job=$$job"$$EXP $(n_constant) randthrp_suppcomp "; \
			job=$$job"-start_mm $$START -stop_mm $$STOP$$BUFFER"; \
			echo $$job; \
			$$job; \
		done; \
	done; \
	for LEN in $$(seq 4 9); do \
		for START in $$(seq 9 $$(( $$LEN_MIR - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="bsub Rscript SolveForKds/FitSiteKds.R $(mirna) "; \
			job=$$job"$$EXP $(n_constant) randthrp_suppcomp "; \
			job=$$job"-start_mm $$START -stop_mm $$STOP$$COMB$$BUFFER"; \
			echo $$job; \
			$$job; \
		done; \
	done)




AssignAllMismatchSitesLet7a :
	@(for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 21 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="make mirna=let-7a-21nt exp=equil_c2_nb "; \
			job=$$job"n_constant=$(n_constant) start_mm=$$START "; \
			job=$$job"stop_mm=$$STOP AssignBipartiteMismatchSites"; \
			$$job; \
		done; \
	done)


AssignAllMismatchSitesMiR1 :
	@(for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 22 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="make mirna=miR-1 exp=equil_c_nb "; \
			job=$$job"n_constant=$(n_constant) start_mm=$$START "; \
			job=$$job"stop_mm=$$STOP AssignBipartiteMismatchSites"; \
			echo $$job; \
			$$job; \
		done; \
	done)

AssignAllMismatchSitesMiR155 :
	@(for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 23 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="make mirna=miR-155 exp=equil_sc_nb "; \
			job=$$job"n_constant=$(n_constant) start_mm=$$START "; \
			job=$$job"stop_mm=$$STOP AssignBipartiteMismatchSites"; \
			echo $$job; \
			$$job; \
		done; \
	done)


AssignAllMismatchSitesLet7aPlus1 :
	@(for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 21 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="make mirna=let-7a_plus1 exp=equil_c_nb "; \
			job=$$job"n_constant=$(n_constant) start_mm=$$START "; \
			job=$$job"stop_mm=$$STOP AssignBipartiteMismatchSites"; \
			$$job; \
		done; \
	done)

AssignAllMismatchSitesLet7aMinus1 :
	@(for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 21 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="make mirna=let-7a_minus1 exp=equil_c_nb "; \
			job=$$job"n_constant=$(n_constant) start_mm=$$START "; \
			job=$$job"stop_mm=$$STOP AssignBipartiteMismatchSites"; \
			$$job; \
		done; \
	done)

AssignAllMismatchSitesLet7aMiR155 :
	@(for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 23 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="make mirna=let-7a_miR-155 exp=equil_c_nb "; \
			job=$$job"n_constant=$(n_constant) start_mm=$$START "; \
			job=$$job"stop_mm=$$STOP AssignBipartiteMismatchSites"; \
			$$job; \
		done; \
	done)

AssignAllMismatchSitesMiR155Let7a :
	@(for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 21 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="make mirna=miR-155_let-7a exp=equil_c_nb "; \
			job=$$job"n_constant=$(n_constant) start_mm=$$START "; \
			job=$$job"stop_mm=$$STOP AssignBipartiteMismatchSites"; \
			$$job; \
		done; \
	done)






FitAllMismatchKdsLet7 :
	@(for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 21 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="python SolveForKds/MakeSiteCountTable.py let-7a-21nt "; \
			job=$$job"equil_c2_nb $(n_constant) progthrp_suppcomp "; \
			job=$$job"-start_mm $$START -stop_mm $$STOP -new "; \
			echo $$job; \
			$$job; \
		done; \
	done; \
	for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 21 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="bsub Rscript SolveForKds/FitSiteKds.R let-7a-21nt "; \
			job=$$job"equil_c2_nb $(n_constant) progthrp_suppcomp "; \
			job=$$job"-start_mm $$START -stop_mm $$STOP -new "; \
			echo $$job; \
			$$job; \
		done; \
	done)

FitAllMismatchKdsMiR1 :
	@(for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 22 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="python SolveForKds/MakeSiteCountTable.py miR-1 "; \
			job=$$job"equil_c_nb $(n_constant) progthrp_suppcomp "; \
			job=$$job"-start_mm $$START -stop_mm $$STOP -new "; \
			echo $$job; \
			$$job; \
		done; \
	done; \
	for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 22 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="bsub Rscript SolveForKds/FitSiteKds.R miR-1 "; \
			job=$$job"equil_c_nb $(n_constant) progthrp_suppcomp "; \
			job=$$job"-start_mm $$START -stop_mm $$STOP -new "; \
			echo $$job; \
			$$job; \
		done; \
	done)

FitAllMismatchKdsMiR155 :
	@(for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 23 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="python SolveForKds/MakeSiteCountTable.py miR-155 "; \
			job=$$job"equil_sc_nb $(n_constant) progthrp_suppcomp "; \
			job=$$job"-start_mm $$START -stop_mm $$STOP -new "; \
			echo $$job; \
			$$job; \
		done; \
	done; \
	for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 23 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="bsub Rscript SolveForKds/FitSiteKds.R miR-155 "; \
			job=$$job"equil_sc_nb $(n_constant) progthrp_suppcomp "; \
			job=$$job"-start_mm $$START -stop_mm $$STOP -new "; \
			echo $$job; \
			$$job; \
		done; \
	done)

FitAllMismatchKdsLet7aPlus1 :
	@(for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 21 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="python SolveForKds/MakeSiteCountTable.py let-7a_plus1 "; \
			job=$$job"equil_c_nb $(n_constant) progthrp_suppcomp "; \
			job=$$job"-start_mm $$START -stop_mm $$STOP -new "; \
			echo $$job; \
			$$job; \
		done; \
	done; \
	for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 21 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="bsub Rscript SolveForKds/FitSiteKds.R let-7a_plus1 "; \
			job=$$job"equil_c_nb $(n_constant) progthrp_suppcomp "; \
			job=$$job"-start_mm $$START -stop_mm $$STOP -new "; \
			echo $$job; \
			$$job; \
		done; \
	done)

FitAllMismatchKdsLet7aMinus1 :
	@(for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 21 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="python SolveForKds/MakeSiteCountTable.py let-7a_minus1 "; \
			job=$$job"equil_c_nb $(n_constant) progthrp_suppcomp "; \
			job=$$job"-start_mm $$START -stop_mm $$STOP -new "; \
			echo $$job; \
			$$job; \
		done; \
	done; \
	for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 21 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="bsub Rscript SolveForKds/FitSiteKds.R let-7a_minus1 "; \
			job=$$job"equil_c_nb $(n_constant) progthrp_suppcomp "; \
			job=$$job"-start_mm $$START -stop_mm $$STOP -new "; \
			echo $$job; \
			$$job; \
		done; \
	done)

FitAllMismatchKdsLet7aMiR155 :
	@(for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 21 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="python SolveForKds/MakeSiteCountTable.py let-7a_miR-155 "; \
			job=$$job"equil_c_nb $(n_constant) progthrp_suppcomp "; \
			job=$$job"-start_mm $$START -stop_mm $$STOP -new "; \
			echo $$job; \
			$$job; \
		done; \
	done; \
	for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 21 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="bsub Rscript SolveForKds/FitSiteKds.R let-7a_miR-155 "; \
			job=$$job"equil_c_nb $(n_constant) progthrp_suppcomp "; \
			job=$$job"-start_mm $$START -stop_mm $$STOP -new "; \
			echo $$job; \
			$$job; \
		done; \
	done)

FitAllMismatchKdsMiR155Let7a :
	@(for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 21 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="python SolveForKds/MakeSiteCountTable.py miR-155_let-7a "; \
			job=$$job"equil_c_nb $(n_constant) progthrp_suppcomp "; \
			job=$$job"-start_mm $$START -stop_mm $$STOP -new "; \
			echo $$job; \
			$$job; \
		done; \
	done; \
	for LEN in $$(seq 4 11); do \
		for START in $$(seq 9 $$(( 21 - $$LEN + 1 ))); do \
			STOP=$$(( $$START + $$LEN - 1 )); \
			job="bsub Rscript SolveForKds/FitSiteKds.R miR-155_let-7a "; \
			job=$$job"equil_c_nb $(n_constant) progthrp_suppcomp "; \
			job=$$job"-start_mm $$START -stop_mm $$STOP -new "; \
			echo $$job; \
			$$job; \
		done; \
	done)







AssignBipartiteSitesRandom :
	@(for CON in $(COND); do \
		job="bsub -q 18 -n 20 -R span[hosts=1] python "; \
		job=$$job"$(HOME)$(DIR_as)$(SCR_absr) $(mirna) $(exp) $$CON "; \
		job=$$job"$(n_constant)"$(BUFFER3P)" -jobs 19"; \
		echo $$job; \
		$$job; \
	done)


# FOR THREE PRIME PAPER
AssessProgrammedKmers :
	@(for CON in $(COND); do \
		job="bsub -q 18 -n 20 -R span[hosts=1] python "; \
		job=$$job"$(HOME)$(DIR_as)$(SCR_akbo) $(mirna) $(exp) $$CON "; \
		job=$$job"$(n_constant) -jobs 19"; \
		echo $$job; \
		$$job; \
	done)

# FOR THREE PRIME PAPER

PreprocessAllEquilProg :
	make mirna=let-7a-21nt exp=equil_c_nb PreProcessReads

AssignAllMirnasProgThrpSites :
	make mirna=let-7a-21nt exp=equil_c_nb n_constant=$(n_constant) AssignBipartiteSites
	make mirna=let-7a-21nt exp=equil_c2_nb n_constant=$(n_constant) AssignBipartiteSites
	#make mirna=let-7a-21nt exp=equil_s_nb n_constant=$(n_constant) AssignBipartiteSites
	#make mirna=let-7a-21nt exp=equil_sc_nb n_constant=$(n_constant) AssignBipartiteSites
	make mirna=miR-1 exp=equil_c_nb n_constant=$(n_constant) AssignBipartiteSites
	#make mirna=miR-1 exp=equil_sc_nb n_constant=$(n_constant) AssignBipartiteSites
	make mirna=let-7a_plus1 exp=equil_c_nb n_constant=$(n_constant) AssignBipartiteSites
	make mirna=let-7a_minus1 exp=equil_c_nb n_constant=$(n_constant) AssignBipartiteSites
	make mirna=miR-155 exp=equil_sc_nb n_constant=$(n_constant) AssignBipartiteSites
	make mirna=let-7a_miR-155 exp=equil_c_nb n_constant=$(n_constant) AssignBipartiteSites
	make mirna=miR-155_let-7a exp=equil_c_nb n_constant=$(n_constant) AssignBipartiteSites

AssignAllMirnasRandThrpSites :
	make mirna=miR-1 exp=equilibrium n_constant=$(n_constant) buffer="-buffer" AssignBipartiteSitesRandom
	make mirna=miR-1 exp=equilibrium n_constant=$(n_constant) buffer="-buffer" AssignBipartiteSitesRandJustCombinedInput
	make mirna=let-7a exp=equilibrium n_constant=$(n_constant) AssignBipartiteSitesRandom
	make mirna=let-7a exp=equilibrium n_constant=$(n_constant) AssignBipartiteSitesRandJustCombinedInput
	make mirna=let-7a-21nt exp=equilibrium_nb n_constant=$(n_constant) AssignBipartiteSitesRandom
	make mirna=let-7a-21nt exp=equilibrium_nb n_constant=$(n_constant) AssignBipartiteSitesRandJustCombinedInput
	make mirna=miR-155 exp=equilibrium n_constant=$(n_constant) AssignBipartiteSitesRandom
	make mirna=miR-155 exp=equilibrium n_constant=$(n_constant) AssignBipartiteSitesRandJustCombinedInput
	make mirna=miR-124 exp=equilibrium n_constant=$(n_constant) AssignBipartiteSitesRandom
	make mirna=miR-124 exp=equilibrium n_constant=$(n_constant) AssignBipartiteSitesRandJustCombinedInput
	make mirna=lsy-6 exp=equilibrium n_constant=$(n_constant) AssignBipartiteSitesRandom
	make mirna=lsy-6 exp=equilibrium n_constant=$(n_constant) AssignBipartiteSitesRandJustCombinedInput
	make mirna=miR-7-23nt exp=equilibrium2_nb n_constant=$(n_constant) AssignBipartiteSitesRandom
	make mirna=miR-7-23nt exp=equilibrium2_nb n_constant=$(n_constant) AssignBipartiteSitesRandJustCombinedInput



# FOR THREE PRIME PAPER
MakeAllMirnasProgThrpSiteCountTables :
	python SolveForKds/MakeSiteCountTable.py let-7a-21nt equil_c_nb $(n_constant) progthrp
	python SolveForKds/MakeSiteCountTable.py let-7a-21nt equil_c2_nb $(n_constant) progthrp
	python SolveForKds/MakeSiteCountTable.py let-7a-21nt equil_s_nb $(n_constant) progthrp
	python SolveForKds/MakeSiteCountTable.py let-7a-21nt equil_sc_nb $(n_constant) progthrp
	python SolveForKds/MakeSiteCountTable.py let-7a_plus1 equil_c_nb $(n_constant) progthrp
	python SolveForKds/MakeSiteCountTable.py let-7a_minus1 equil_c_nb $(n_constant) progthrp
	python SolveForKds/MakeSiteCountTable.py let-7a_miR-155 equil_c_nb $(n_constant) progthrp
	python SolveForKds/MakeSiteCountTable.py miR-1 equil_c_nb $(n_constant) progthrp
	python SolveForKds/MakeSiteCountTable.py miR-1 equil_sc_nb $(n_constant) progthrp
	python SolveForKds/MakeSiteCountTable.py miR-155 equil_sc_nb $(n_constant) progthrp
	python SolveForKds/MakeSiteCountTable.py miR-155_let-7a equil_c_nb $(n_constant) progthrp

	python SolveForKds/MakeSiteCountTable.py let-7a-21nt equil_c_nb $(n_constant) progthrp_suppcomp
	python SolveForKds/MakeSiteCountTable.py let-7a-21nt equil_c2_nb $(n_constant) progthrp_suppcomp
	python SolveForKds/MakeSiteCountTable.py let-7a-21nt equil_s_nb $(n_constant) progthrp_suppcomp
	python SolveForKds/MakeSiteCountTable.py let-7a-21nt equil_sc_nb $(n_constant) progthrp_suppcomp
	python SolveForKds/MakeSiteCountTable.py let-7a_plus1 equil_c_nb $(n_constant) progthrp_suppcomp
	python SolveForKds/MakeSiteCountTable.py let-7a_minus1 equil_c_nb $(n_constant) progthrp_suppcomp
	python SolveForKds/MakeSiteCountTable.py let-7a_miR-155 equil_c_nb $(n_constant) progthrp_suppcomp
	python SolveForKds/MakeSiteCountTable.py miR-1 equil_c_nb $(n_constant) progthrp_suppcomp
	python SolveForKds/MakeSiteCountTable.py miR-1 equil_sc_nb $(n_constant) progthrp_suppcomp
	python SolveForKds/MakeSiteCountTable.py miR-155 equil_sc_nb $(n_constant) progthrp_suppcomp
	python SolveForKds/MakeSiteCountTable.py miR-155_let-7a equil_c_nb $(n_constant) progthrp_suppcomp



MakeAllMirnasRandThrpSiteCountTables :
	python SolveForKds/MakeSiteCountTable.py miR-1 equilibrium $(n_constant) randthrp -buffer 
	python SolveForKds/MakeSiteCountTable.py let-7a equilibrium $(n_constant) randthrp
	python SolveForKds/MakeSiteCountTable.py miR-155 equilibrium $(n_constant) randthrp
	python SolveForKds/MakeSiteCountTable.py miR-124 equilibrium $(n_constant) randthrp
	python SolveForKds/MakeSiteCountTable.py lsy-6 equilibrium $(n_constant) randthrp
	python SolveForKds/MakeSiteCountTable.py miR-7-23nt equilibrium2_nb $(n_constant) randthrp

	python SolveForKds/MakeSiteCountTable.py miR-1 equilibrium $(n_constant) randthrp_suppcomp -buffer 
	python SolveForKds/MakeSiteCountTable.py let-7a equilibrium $(n_constant) randthrp_suppcomp
	python SolveForKds/MakeSiteCountTable.py miR-155 equilibrium $(n_constant) randthrp_suppcomp
	python SolveForKds/MakeSiteCountTable.py miR-124 equilibrium $(n_constant) randthrp_suppcomp
	python SolveForKds/MakeSiteCountTable.py lsy-6 equilibrium $(n_constant) randthrp_suppcomp
	python SolveForKds/MakeSiteCountTable.py miR-7-23nt equilibrium2_nb $(n_constant) randthrp_suppcomp

	python SolveForKds/MakeSiteCountTable.py miR-1 equilibrium $(n_constant) randthrp_comp -buffer 
	python SolveForKds/MakeSiteCountTable.py let-7a equilibrium $(n_constant) randthrp_comp
	python SolveForKds/MakeSiteCountTable.py miR-155 equilibrium $(n_constant) randthrp_comp
	python SolveForKds/MakeSiteCountTable.py miR-124 equilibrium $(n_constant) randthrp_comp
	python SolveForKds/MakeSiteCountTable.py lsy-6 equilibrium $(n_constant) randthrp_comp
	python SolveForKds/MakeSiteCountTable.py miR-7-23nt equilibrium2_nb $(n_constant) randthrp_comp


FitAllMirnasRandThrpKds :
	bsub Rscript SolveForKds/FitSiteKds.R miR-1 equilibrium $(n_constant) randthrp -nocombI -buffer
	bsub Rscript SolveForKds/FitSiteKds.R let-7a equilibrium $(n_constant) randthrp
	bsub Rscript SolveForKds/FitSiteKds.R miR-155 equilibrium $(n_constant) randthrp
	bsub Rscript SolveForKds/FitSiteKds.R miR-124 equilibrium $(n_constant) randthrp
	bsub Rscript SolveForKds/FitSiteKds.R lsy-6 equilibrium $(n_constant) randthrp
	bsub Rscript SolveForKds/FitSiteKds.R miR-7-23nt equilibrium2_nb $(n_constant) randthrp -nocombI

	bsub Rscript SolveForKds/FitSiteKds.R miR-1 equilibrium $(n_constant) randthrp_suppcomp -nocombI -buffer
	bsub Rscript SolveForKds/FitSiteKds.R let-7a equilibrium $(n_constant) randthrp_suppcomp
	bsub Rscript SolveForKds/FitSiteKds.R miR-155 equilibrium $(n_constant) randthrp_suppcomp
	bsub Rscript SolveForKds/FitSiteKds.R miR-124 equilibrium $(n_constant) randthrp_suppcomp
	bsub Rscript SolveForKds/FitSiteKds.R lsy-6 equilibrium $(n_constant) randthrp_suppcomp
	bsub Rscript SolveForKds/FitSiteKds.R miR-7-23nt equilibrium2_nb $(n_constant) randthrp_suppcomp -nocombI

	bsub Rscript SolveForKds/FitSiteKds.R miR-1 equilibrium $(n_constant) randthrp_suppcomp -nocombI -buffer -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R let-7a equilibrium $(n_constant) randthrp_suppcomp -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R miR-155 equilibrium $(n_constant) randthrp_suppcomp -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R miR-124 equilibrium $(n_constant) randthrp_suppcomp -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R lsy-6 equilibrium $(n_constant) randthrp_suppcomp -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R miR-7-23nt equilibrium2_nb $(n_constant) randthrp_suppcomp -nocombI -sumseed

	bsub Rscript SolveForKds/FitSiteKds.R miR-1 equilibrium $(n_constant) randthrp_comp -nocombI -buffer
	bsub Rscript SolveForKds/FitSiteKds.R let-7a equilibrium $(n_constant) randthrp_comp
	bsub Rscript SolveForKds/FitSiteKds.R miR-155 equilibrium $(n_constant) randthrp_comp
	bsub Rscript SolveForKds/FitSiteKds.R miR-124 equilibrium $(n_constant) randthrp_comp
	bsub Rscript SolveForKds/FitSiteKds.R lsy-6 equilibrium $(n_constant) randthrp_comp
	bsub Rscript SolveForKds/FitSiteKds.R miR-7-23nt equilibrium2_nb $(n_constant) randthrp_comp -nocombI

	bsub Rscript SolveForKds/FitSiteKds.R miR-1 equilibrium $(n_constant) randthrp_comp -nocombI -buffer -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R let-7a equilibrium $(n_constant) randthrp_comp -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R miR-155 equilibrium $(n_constant) randthrp_comp -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R miR-124 equilibrium $(n_constant) randthrp_comp -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R lsy-6 equilibrium $(n_constant) randthrp_comp -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R miR-7-23nt equilibrium2_nb $(n_constant) randthrp_comp -nocombI -sumseed


FitAllMirnasProgThrpKds :
	bsub Rscript SolveForKds/FitSiteKds.R let-7a-21nt equil_c_nb $(n_constant) progthrp
	bsub Rscript SolveForKds/FitSiteKds.R let-7a-21nt equil_c2_nb $(n_constant) progthrp
	bsub Rscript SolveForKds/FitSiteKds.R let-7a-21nt equil_s_nb $(n_constant) progthrp
	bsub Rscript SolveForKds/FitSiteKds.R let-7a-21nt equil_sc_nb $(n_constant) progthrp
	bsub Rscript SolveForKds/FitSiteKds.R let-7a_plus1 equil_c_nb $(n_constant) progthrp
	bsub Rscript SolveForKds/FitSiteKds.R let-7a_minus1 equil_c_nb $(n_constant) progthrp
	bsub Rscript SolveForKds/FitSiteKds.R let-7a_miR-155 equil_c_nb $(n_constant) progthrp
	bsub Rscript SolveForKds/FitSiteKds.R miR-1 equil_c_nb $(n_constant) progthrp
	bsub Rscript SolveForKds/FitSiteKds.R miR-1 equil_sc_nb $(n_constant) progthrp
	bsub Rscript SolveForKds/FitSiteKds.R miR-155 equil_sc_nb $(n_constant) progthrp
	bsub Rscript SolveForKds/FitSiteKds.R miR-155_let-7a equil_c_nb $(n_constant) progthrp

	bsub Rscript SolveForKds/FitSiteKds.R let-7a-21nt equil_c_nb $(n_constant) progthrp_suppcomp
	bsub Rscript SolveForKds/FitSiteKds.R let-7a-21nt equil_c2_nb $(n_constant) progthrp_suppcomp
	bsub Rscript SolveForKds/FitSiteKds.R let-7a-21nt equil_s_nb $(n_constant) progthrp_suppcomp
	bsub Rscript SolveForKds/FitSiteKds.R let-7a-21nt equil_sc_nb $(n_constant) progthrp_suppcomp
	bsub Rscript SolveForKds/FitSiteKds.R let-7a_plus1 equil_c_nb $(n_constant) progthrp_suppcomp
	bsub Rscript SolveForKds/FitSiteKds.R let-7a_minus1 equil_c_nb $(n_constant) progthrp_suppcomp
	bsub Rscript SolveForKds/FitSiteKds.R let-7a_miR-155 equil_c_nb $(n_constant) progthrp_suppcomp
	bsub Rscript SolveForKds/FitSiteKds.R miR-1 equil_c_nb $(n_constant) progthrp_suppcomp
	bsub Rscript SolveForKds/FitSiteKds.R miR-1 equil_sc_nb $(n_constant) progthrp_suppcomp
	bsub Rscript SolveForKds/FitSiteKds.R miR-155 equil_sc_nb $(n_constant) progthrp_suppcomp
	bsub Rscript SolveForKds/FitSiteKds.R miR-155_let-7a equil_c_nb $(n_constant) progthrp_suppcomp

	bsub Rscript SolveForKds/FitSiteKds.R let-7a-21nt equil_c_nb $(n_constant) progthrp_suppcomp -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R let-7a-21nt equil_c2_nb $(n_constant) progthrp_suppcomp -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R let-7a-21nt equil_s_nb $(n_constant) progthrp_suppcomp -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R let-7a-21nt equil_sc_nb $(n_constant) progthrp_suppcomp -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R let-7a_plus1 equil_c_nb $(n_constant) progthrp_suppcomp -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R let-7a_minus1 equil_c_nb $(n_constant) progthrp_suppcomp -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R let-7a_miR-155 equil_c_nb $(n_constant) progthrp_suppcomp -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R miR-1 equil_c_nb $(n_constant) progthrp_suppcomp -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R miR-1 equil_sc_nb $(n_constant) progthrp_suppcomp -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R miR-155 equil_sc_nb $(n_constant) progthrp_suppcomp -sumseed
	bsub Rscript SolveForKds/FitSiteKds.R miR-155_let-7a equil_c_nb $(n_constant) progthrp_suppcomp -sumseed





# FOR THREE PRIME PAPER
FitPositionalProgrammedSiteKds :
	if [ "$(mirna)" = "miR-155" ]; then \
		echo "hi"; \
		LEN_MIR=$$(( 23 )); \
		echo $$LEN_MIR; \
		MIRNA_SITE_FILE=$(mirna); \
	elif [ "$(mirna)" = "miR-155_let-7a" ]; then \
		LEN_MIR=21; \
		MIRNA_SITE_FILE="miR-155"; \
	elif [ "$(mirna)" = "miR-1" ]; then \
		LEN_MIR=22; \
		MIRNA_SITE_FILE=$(mirna); \
	elif [ "$(mirna)" = "let-7a_miR-155" ]; then \
		LEN_MIR=23; \
		MIRNA_SITE_FILE="let-7a-21nt"; \
	else \
		LEN_MIR=21; \
		MIRNA_SITE_FILE="let-7a-21nt"; \
		echo "$(mirna)"; \
	fi; \
	seed_site_path="/lab/bartel1_ata/mcgeary/computation/AgoRBNS"; \
	seed_site_path=$$seed_site_path"/AssignSiteTypes/"; \
	seed_site_path=$$seed_site_path"site_categories/programmedbase/sites."; \
	seed_site_path=$$seed_site_path"$$MIRNA_SITE_FILE"; \
	seed_site_path=$$seed_site_path"_programmedbase.txt"; \
	# The second condition of the OR gate allows the last line of the file to \
	# be read by the unix script, which otherwise wouldn't be read because \
	# there is no newline character. \
	while read -r SITE || [ -n "$$SITE" ]; do \
		job="bsub Rscript $(HOME)$(DIR_kd)$(SCR_ppskd) $(mirna) $(exp)"; \
		job=$$job" $(n_constant) $$SITE"; \
		echo $$job; \
		$$job; \
	done <$$seed_site_path; \
	echo $(LEN_MAX); \
	for LEN in $$(seq 4 11); do \
		NUM_K=$$(( $$LEN_MIR - 8 - $$LEN + 1 )); \
		for START in $$(seq 9 $$(( 9 + $$NUM_K - 1 ))); do \
			END=$$(( $$START + $$LEN - 1 )); \
			SITE=$$LEN"mer-m"$$START"."$$END; \
			job="bsub Rscript $(HOME)$(DIR_kd)$(SCR_ppskd) $(mirna) $(exp)"; \
			job=$$job" $(n_constant) $$SITE"; \
			echo $$job; \
			$$job; \
		done; \
	done


FitPositionalMismatchSiteKds :
	if [ "$(mirna)" = "miR-155" ]; then \
		echo "hi"; \
		LEN_MIR=$$(( 23 )); \
		echo $$LEN_MIR; \
		MIRNA_SITE_FILE=$(mirna); \
	elif [ "$(mirna)" = "miR-155_let-7a" ]; then \
		LEN_MIR=21; \
		MIRNA_SITE_FILE="miR-155"; \
	elif [ "$(mirna)" = "miR-1" ]; then \
		LEN_MIR=22; \
		MIRNA_SITE_FILE=$(mirna); \
	elif [ "$(mirna)" = "let-7a_miR-155" ]; then \
		LEN_MIR=23; \
		MIRNA_SITE_FILE="let-7a-21nt"; \
	else \
		LEN_MIR=21; \
		MIRNA_SITE_FILE="let-7a-21nt"; \
		echo "$(mirna)"; \
	fi; \
	seed_site_path="/lab/bartel1_ata/mcgeary/computation/AgoRBNS"; \
	seed_site_path=$$seed_site_path"/AssignSiteTypes/"; \
	seed_site_path=$$seed_site_path"site_categories/programmedbase/sites."; \
	seed_site_path=$$seed_site_path"$$MIRNA_SITE_FILE"; \
	seed_site_path=$$seed_site_path"_programmedbase.txt"; \
	# The second condition of the OR gate allows the last line of the file to \
	# be read by the unix script, which otherwise wouldn't be read because \
	# there is no newline character. \
	while read -r SITE || [ -n "$$SITE" ]; do \
		job="bsub Rscript $(HOME)$(DIR_kd)$(SCR_pmmkd) $(mirna) $(exp)"; \
		job=$$job" $(n_constant) $$SITE $(lambda)"; \
		echo $$job; \
		$$job; \
	done <$$seed_site_path; \
	echo $(LEN_MAX); \
	for LEN in $$(seq 4 11); do \
		NUM_K=$$(( $$LEN_MIR - 8 - $$LEN + 1 )); \
		for START in $$(seq 9 $$(( 9 + $$NUM_K - 1 ))); do \
			END=$$(( $$START + $$LEN - 1 )); \
			SITE=$$LEN"mer-m"$$START"."$$END; \
			job="bsub Rscript $(HOME)$(DIR_kd)$(SCR_pmmkd) $(mirna) $(exp)"; \
			job=$$job" $(n_constant) $$SITE $(lambda)"; \
			echo $$job; \
			$$job; \
		done; \
	done
	


AssignCompetitorOligo :
	@(for NUM in $$(seq 4 16); do \
		job="bsub -q 18 -n 20 -R span[hosts=1] python $(HOME)$(DIR_as)$(SCR_ack) $(mirna)"; \
		job=$$job" $(exp) $(cond) $(n_constant) $$NUM $(off)"$(ADDCOMP)$(WOB)$(ALTMIRCOMP)" -jobs 19"; \
		echo $$job; \
		$$job;\
	done)

AssignCompetitorOligoConditionRange :
	@(for CON in $(COND_mirna); do \
		make mirna=$(mirna) exp=$(exp) cond=$$CON n_constant=$(n_constant) \
		addcomp=$(addcomp) off=$(off) altmircomp=$(altmircomp) AssignCompetitorOligo; \
	done)



MakeSiteCountTables :
	@(job=""; \
	for MIRNA in $(MIRNAe); do \
		if [ "$$MIRNA" = "miR-7-23nt" ]; then \
			EXP='equilibrium2_nb'; \
		else \
			EXP='equilibrium'; \
		fi; \
		job=""; \
		echo $$MIRNA; \
		job=$$job"python $(HOME)$(DIR_kd)$(SCR_sc) $$MIRNA $$EXP "; \
		job=$$job"$(n_constant) $(sitelist)"; \
		echo $$job; \
		$$job; \
	done)

MakeAllPaperSiteCountTables :
	for SITELIST in canonical paperfinal paperextendedfinal baek centered11 bulge del; do \
		for MIRNA in $(MIRNAe); do \
			if [ "$$MIRNA" = "miR-7-23nt" ]; then \
				EXP='equilibrium2_nb'; \
				comb=" -nocombI"; \
				buffer=""; \
			elif [ "$$MIRNA" = "miR-1" ]; then \
				EXP='equilibrium'; \
				comb=""; \
				buffer=" -buffer"; \
			else \
				EXP='equilibrium'; \
				comb=""; \
				buffer=""; \
			fi; \
			job="python $(HOME)$(DIR_kd)$(SCR_sc) $$MIRNA $$EXP 5 $$SITELIST$$buffer"; \
			echo $$job; \
			$$job; \
			job="python $(HOME)$(DIR_kd)$(SCR_msc) $$MIRNA $$EXP 5 $$SITELIST$$buffer"; \
			echo $$job; \
		done; \
	done



FitAllPaperKds :
	for SITELIST in canonical paperfinal paperextendedfinal baek centered11 bulge del; do \
		for MIRNA in $(MIRNAe); do \
			if [ "$$MIRNA" = "miR-7-23nt" ]; then \
				EXP='equilibrium2_nb'; \
				comb=" -nocombI"; \
				buffer=""; \
			elif [ "$$MIRNA" = "miR-1" ]; then \
				EXP='equilibrium'; \
				comb=" -nocombI"; \
				buffer=" -buffer"; \
			else \
				EXP='equilibrium'; \
				comb=""; \
				buffer=""; \
			fi; \
			job="Rscript $(HOME)$(DIR_kd)$(SCR_skd) $$MIRNA $$EXP 5 $$SITELIST$$comb$$buffer"; \
			echo $$job; \
			bsub $$job; \
			job="Rscript $(HOME)$(DIR_kd)$(SCR_skd) $$MIRNA $$EXP 5 $$SITELIST$$comb$$buffer -single"; \
			echo $$job; \
			bsub $$job; \
		done; \
	done



MakeMultiSiteCountTables :
	@(job=""; \
	for MIRNA in $(MIRNAe); do \
		if [ "$$MIRNA" = "miR-7-23nt" ]; then \
			EXP='equilibrium2_nb'; \
		else \
			EXP='equilibrium'; \
		fi; \
		job=""; \
		echo $$MIRNA; \
		job=$$job"python $(HOME)$(DIR_kd)$(SCR_msc) $$MIRNA $$EXP "; \
		job=$$job"$(n_constant) $(sitelist)"; \
		echo $$job; \
		$$job; \
	done)




AssignBipartiteSitesRandJustCombinedInput :
	@(job="bsub -q 18 -n 20 -R span[hosts=1] python $(HOME)$(DIR_as)"; \
	  job=$$job"$(SCR_absr) $(mirna) $(exp) I_combined $(n_constant) "; \
	  job=$$job$(BUFFER3P)" -jobs 19"; \
	  echo $$job; \
	  $$job)


AssignRandomBipartiteMismatchSitesJustCombinedInput :
	@(job="bsub -q 18 -n 20 -R span[hosts=1] python $(HOME)$(DIR_as)"; \
	  job=$$job"$(SCR_abmsr) $(mirna) $(exp) I_combined $(n_constant) "; \
	  job=$$job"$(start_mm) $(stop_mm) "$(BUFFER3P)" -jobs 19"; \
	  echo $$job; \
	  $$job)

# 	python $(HOME)general/SubmitJob.py "$$job")

AssignSitesJustCombinedInputKinetics :
	@(job="python $(HOME)$(DIR_as)$(SCR_as) $(mirna) $(exp) I_combined "; \
	job=$$job"$(n_constant) $(sitelist)"$(BUFFER3P)$(N_EX)$(UNIQUE)"\n"; \
	echo $$job; \
	python $(HOME)general/SubmitJob.py "$$job")


# Assign12mers :
# 	for MIRSTART in $$(seq 1 5); do \
# 		for CON in $(COND); do \
# 			echo python $(HOME)$(DIR_as)$(SCR_as) $(mirna) $(exp) $$CON\
# 				$(n_constant) 12mers -mir_start $$MIRSTART;\
# 				bsub -n 21 python $(HOME)$(DIR_as)$(SCR_as) $(mirna) $(exp)\
# 				$$CON $(n_constant) 12mers -mir_start $$MIRSTART;\
# 		done; \
# 		echo python $(HOME)$(DIR_as)$(SCR_as) $(mirna)\
# 			$(exp) I_combined $(n_constant) 12mers -mir_start $$MIRSTART;\
# 		bsub -n 61 -m bigboy python $(HOME)$(DIR_as)$(SCR_as) $(mirna)\
# 			$(exp) I_combined $(n_constant) 12mers -mir_start $$MIRSTART;\
# 		done)



Assign12mers :
	@(job=""; \
	for MIRSTART in $$(seq 1 5); do \
		for CON in $(COND); do \
			job=$$job"python $(HOME)$(DIR_as)$(SCR_as) $(mirna) $(exp) $$CON "; \
			job=$$job"$(n_constant) 12mers -mir_start $$MIRSTART\n"; \
		done; \
		job=$$job"python $(HOME)$(DIR_as)$(SCR_as) $(mirna) $(exp) I_combined "; \
		job=$$job"$(n_constant) 12mers -mir_start $$MIRSTART\n"; \
	done; \
	echo $$job; \
	python $(HOME)general/SubmitJob.py "$$job")


Assign12mersSingle :
	@(job=""; \
	for CON in $(COND); do \
		job=$$job"python $(HOME)$(DIR_as)$(SCR_as) $(mirna) $(exp) $$CON "; \
		job=$$job"$(n_constant) 12mers -mir_start $(mir_start)\n"; \
	done; \
	job=$$job"python $(HOME)$(DIR_as)$(SCR_as) $(mirna) $(exp) I_combined "; \
	job=$$job"$(n_constant) 12mers -mir_start $(mir_start)\n"; \
	echo $$job; \
	python $(HOME)general/SubmitJob.py "$$job")



Collect12mers :
	for MIRSTART in $$(seq 1 5); do \
		python $(HOME)$(DIR_kd)$(SCR_sc) $(mirna) $(exp) $(n_constant) \
			12mers -mir_start $$MIRSTART; \
	done

CollectAll12mers :
	for MIRNA in $(MIRNAe); do \
		make mirna=$$MIRNA exp=equilibrium n_constant=$(n_constant) \
			Collect12mers; \
	done
	for MIRNA in $(MIRNAeMiR7_2); do \
		make mirna=$$MIRNA exp=equilibrium2_nb n_constant=$(n_constant) \
			Collect12mers; \
	done

Assign16mers :
	for START in $$(seq 1 5); do \
		for SPLIT in left right;\
			do { for CON in $(COND);\
				do { echo python $(HOME)$(DIR_as)$(SCR_as) $(mirna) $(exp)
					$$CON $(n_constant) 16mers -mir_start $$START -split16\
					$$SPLIT;\
				bsub -n 21 python $(HOME)$(DIR_as)$(SCR_as) $(mirna) $(exp)\
					$$CON $(n_constant) 16mers -mir_start $$START -split16\
					$SPLIT;}\
				done;\
				echo python $(HOME)$(DIR_as)$(SCR_as) $(mirna) $(exp)\
					I_combined $(n_constant) 16mers -mir_start $$START -split16\
					right;\
				bsub -n 61 -m bigboy python $(HOME)$(DIR_as)$(SCR_as) $(mirna)\
				$(exp) I_combined $(n_constant) 16mers -mir_start $$START\
				-split16 right;}\
			done

Collect16mers :
	for MIRSTART in $$(seq 1 5); do \
		for SPLIT in left right; do \
			python $(HOME)$(DIR_kd)$(SCR_sc) $(mirna) $(exp) $(n_constant) \
				16mers -mir_start $$MIRSTART -split16 $$SPLIT; \
		done \
	done

CollectAll16mers :
	for MIRNA in $(MIRNAe); do \
		make mirna=$$MIRNA exp=equilibrium n_constant=$(n_constant) \
			Collect16mers; \
	done
	for MIRNA in $(MIRNAeMiR7_2); do \
		make mirna=$$MIRNA exp=equilibrium2_nb n_constant=$(n_constant) \
			Collect16mers; \
	done




# AssignFlanksJustCombinedInput :
# 	@(job="python $(HOME)$(DIR_as)$(SCR_af) $(mirna) $(exp) I_combined "; \
# 	job=$$job"$(n_constant) $(sitelist)"$(BUFFER3P)"\n"; \
# 	echo $$job; \
# 	python $(HOME)general/SubmitJob.py "$$job")


AssignFlanksJustCombinedInput :
	@(job="bsub -q 18 -n 20 -R span[hosts=1] python $(HOME)$(DIR_as)"; \
	  job=$$job"$(SCR_af) $(mirna) $(exp) I_combined $(n_constant) "; \
	  job=$$job"$(sitelist)"$(BUFFER3P)" -jobs 19"; \
	  echo $$job; \
	  $$job)




MakeInputCountDistributions :
	for CON in $(COND); do \
		bsub $(HOME)general/MakeReadDistributionHistogram.sh $(mirna) $(exp) $$CON; \
	done;


AssignKmersSingle :
	@(for CON in $(COND); do \
		job="bsub -q 18 -n 20 -R span[hosts=1] python $(HOME)$(DIR_as)$(SCR_ak) $(mirna)"; \
		job=$$job" $(exp) $$CON $(n_constant) $(len_k) -n_ex $(n_ex) -jobs 19"; \
		echo $$job; \
		$$job;\
	done)



AssignKmersResub :
	COND="I 40"; \
	for CON in $$COND; do \
		job="bsub -q 18 -n 20 -R span[hosts=1] python "; \
		job=$$job"$(HOME)$(DIR_as)$(SCR_akr) $(mirna) $(exp) $$CON 5 10 "; \
		job=$$job"-n_ex $(n_ex)$(BUFFER3p) -jobs 19"; \
		echo $$job; \
		$$job; \
	done; \
	)

AssignKmersAllWithStart :
	if [ "$(mirna)" = "miR-7-23nt" ] || [ "$(mirna)" = "miR-7-24nt" ]; then \
		MIRNA_NAME="miR-7"; \
	else \
		MIRNA_NAME="$(mirna)"; \
	fi; \
	if [ "$(exp)" = "equilibrium_2_tp" ] || [ "$(exp)" = "equilibrium_tp" ]; then \
		CUTOFF="2"; \
	else \
		CUTOFF=""; \
	fi; \
	if [ "$(mirna)" = "miR-1" ] && [ "$(exp)" = "equilibrium" ]; then \
		BUFFER=" -buffer"; \
	else \
		BUFFER=""; \
	fi; \
	COND="I 40"; \
	SITES_PATH="$(HOME)$(DIR_as)site_categories/papercutoff"$$CUTOFF"/sites."; \
	SITES_PATH=$$SITES_PATH$$MIRNA_NAME"_papercutoff"$$CUTOFF".txt"; \
	echo $$SITES_PATH; \
	LEN_SITE=`wc -l $$SITES_PATH | cut -d " " -f 1`; \
	echo $$LEN_SITE; \
	LEN_SITE=$$(( $$LEN_SITE + 1 )); \
	echo $(n_start); \
	for N_X in $$(seq $(n_start) $$LEN_SITE); do \
	echo $$N_X; \
		for CON in $$COND; do \
			job="bsub -q 18 -n 20 python $(HOME)$(DIR_as)$(SCR_akr) $(mirna) $(exp)"; \
			job="bsub -q 18 -n 20 -R span[hosts=1] python "; \
			job=$$job"$(HOME)$(DIR_as)$(SCR_ak) $(mirna) $(exp) $$CON 5 10 "; \
			job=$$job"-n_ex $$N_X$$BUFFER -jobs 19"; \
			echo $$job; \
			$$job; \
		done; \
	done; \
	)


AssignKmersAllWithStartStop :
	if [ "$(mirna)" = "miR-7-23nt" ] || [ "$(mirna)" = "miR-7-24nt" ]; then \
		MIRNA_NAME="miR-7"; \
	else \
		MIRNA_NAME="$(mirna)"; \
	fi; \
	if [ "$(exp)" = "equilibrium_2_tp" ] || [ "$(exp)" = "equilibrium_tp" ]; then \
		CUTOFF="2"; \
	else \
		CUTOFF=""; \
	fi; \
	if [ "$(mirna)" = "miR-1" ] && [ "$(exp)" = "equilibrium" ]; then \
		BUFFER=" -buffer"; \
	else \
		BUFFER=""; \
	fi; \
	COND="I 40"; \
	SITES_PATH="$(HOME)$(DIR_as)site_categories/papercutoff"$$CUTOFF"/sites."; \
	SITES_PATH=$$SITES_PATH$$MIRNA_NAME"_papercutoff"$$CUTOFF".txt"; \
	echo $$SITES_PATH; \
	LEN_SITE=`wc -l $$SITES_PATH | cut -d " " -f 1`; \
	echo $$LEN_SITE; \
	LEN_SITE=$$(( $$LEN_SITE + 1 )); \
	echo $(n_start); \
	for N_X in $$(seq $(n_start) $(n_stop)); do \
	echo $$N_X; \
		for CON in $$COND; do \
			job="bsub -q 18 -n 20 -R span[hosts=1] python "; \
			job=$$job"$(HOME)$(DIR_as)$(SCR_ak) $(mirna) $(exp) $$CON "; \
			job=$$job"$(n_constant) $(len_k) -n_ex $$N_X$$BUFFER -jobs 19"; \
			echo $$job; \
			$$job; \
		done; \
	done; \
	)





AssignKmersAllResub :
	if [ "$(mirna)" = "miR-7-23nt" ] || [ "$(mirna)" = "miR-7-24nt" ]; then \
		MIRNA_NAME="miR-7"; \
	else \
		MIRNA_NAME="$(mirna)"; \
	fi; \
	if [ "$(exp)" = "equilibrium_2_tp" ] || [ "$(exp)" = "equilibrium_tp" ]; then \
		CUTOFF="2"; \
	else \
		CUTOFF=""; \
	fi; \
	if [ "$(mirna)" = "miR-1" ] && [ "$(exp)" = "equilibrium" ]; then \
		BUFFER=" -buffer"; \
	else \
		BUFFER=""; \
	fi; \
	COND="I 40"; \
	SITES_PATH="$(HOME)$(DIR_as)site_categories/papercutoff"$$CUTOFF"/sites."; \
	SITES_PATH=$$SITES_PATH$$MIRNA_NAME"_papercutoff"$$CUTOFF".txt"; \
	echo $$SITES_PATH; \
	LEN_SITE=`wc -l $$SITES_PATH | cut -d " " -f 1`; \
	echo $$LEN_SITE; \
	LEN_SITE=$$(( $$LEN_SITE + 1 )); \
	for N_X in $$(seq 0 $$LEN_SITE); do \
		for CON in $$COND; do \
			job="bsub -q 18 -n 20 python $(HOME)$(DIR_as)$(SCR_akr) $(mirna) $(exp)"; \
			job="bsub -q 18 -n 20 -R span[hosts=1] python "; \
			job=$$job"$(HOME)$(DIR_as)$(SCR_akr) $(mirna) $(exp) $$CON 5 10 "; \
			job=$$job"-n_ex $$N_X$$BUFFER -jobs 19"; \
			echo $$job; \
			$$job; \
		done; \
	done; \
	)
	






PreprocessThreePrimeHighThroughputLibrary :
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_3p_tp duplex_series -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_3p_tp duplex_series -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a-21nt twist_reporter_assay_3p_tp duplex_series -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_3p_tp no_duplex_series -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_3p_tp no_duplex_series -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a-21nt twist_reporter_assay_3p_tp no_duplex_series -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_3p_tp duplex_parallel -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_3p_tp duplex_parallel -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a-21nt twist_reporter_assay_3p_tp duplex_parallel -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_3p_tp no_duplex_parallel -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_3p_tp no_duplex_parallel -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a-21nt twist_reporter_assay_3p_tp no_duplex_parallel -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_3p_tp duplex_series -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_3p_tp duplex_series -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a-21nt twist_reporter_assay_3p_tp duplex_series -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_3p_tp no_duplex_series -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_3p_tp no_duplex_series -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a-21nt twist_reporter_assay_3p_tp no_duplex_series -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_3p_tp duplex_parallel -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_3p_tp duplex_parallel -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a-21nt twist_reporter_assay_3p_tp duplex_parallel -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_3p_tp no_duplex_parallel -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_3p_tp no_duplex_parallel -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a-21nt twist_reporter_assay_3p_tp no_duplex_parallel -rep 2

PreprocessThreePrimeHighThroughputLibrary2 :
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_3p_2_tp duplex_series -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_3p_2_tp duplex_series -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a-21nt twist_reporter_assay_3p_2_tp duplex_series -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_3p_2_tp no_duplex_series -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_3p_2_tp no_duplex_series -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a-21nt twist_reporter_assay_3p_2_tp no_duplex_series -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_3p_2_tp duplex_parallel -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_3p_2_tp duplex_parallel -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a-21nt twist_reporter_assay_3p_2_tp duplex_parallel -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_3p_2_tp no_duplex_parallel -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_3p_2_tp no_duplex_parallel -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a-21nt twist_reporter_assay_3p_2_tp no_duplex_parallel -rep 1
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_3p_2_tp duplex_series -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_3p_2_tp duplex_series -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a-21nt twist_reporter_assay_3p_2_tp duplex_series -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_3p_2_tp no_duplex_series -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_3p_2_tp no_duplex_series -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a-21nt twist_reporter_assay_3p_2_tp no_duplex_series -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_3p_2_tp duplex_parallel -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_3p_2_tp duplex_parallel -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a-21nt twist_reporter_assay_3p_2_tp duplex_parallel -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py miR-1 twist_reporter_assay_3p_2_tp no_duplex_parallel -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a twist_reporter_assay_3p_2_tp no_duplex_parallel -rep 2
	bsub -q 18 -n 20 -R span[hosts=1] python PreProcessReads/MakeReadFile.py let-7a-21nt twist_reporter_assay_3p_2_tp no_duplex_parallel -rep 2




CountThreePrimeHighThroughputLibraryVariants :
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py miR-1 twist_reporter_assay_3p_tp duplex_series 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a twist_reporter_assay_3p_tp duplex_series 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a-21nt twist_reporter_assay_3p_tp duplex_series 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py miR-1 twist_reporter_assay_3p_tp no_duplex_series 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a twist_reporter_assay_3p_tp no_duplex_series 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a-21nt twist_reporter_assay_3p_tp no_duplex_series 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py miR-1 twist_reporter_assay_3p_tp duplex_parallel 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a twist_reporter_assay_3p_tp duplex_parallel 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a-21nt twist_reporter_assay_3p_tp duplex_parallel 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py miR-1 twist_reporter_assay_3p_tp no_duplex_parallel 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a twist_reporter_assay_3p_tp no_duplex_parallel 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a-21nt twist_reporter_assay_3p_tp no_duplex_parallel 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py miR-1 twist_reporter_assay_3p_tp duplex_series 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a twist_reporter_assay_3p_tp duplex_series 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a-21nt twist_reporter_assay_3p_tp duplex_series 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py miR-1 twist_reporter_assay_3p_tp no_duplex_series 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a twist_reporter_assay_3p_tp no_duplex_series 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a-21nt twist_reporter_assay_3p_tp no_duplex_series 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py miR-1 twist_reporter_assay_3p_tp duplex_parallel 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a twist_reporter_assay_3p_tp duplex_parallel 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a-21nt twist_reporter_assay_3p_tp duplex_parallel 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py miR-1 twist_reporter_assay_3p_tp no_duplex_parallel 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a twist_reporter_assay_3p_tp no_duplex_parallel 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a-21nt twist_reporter_assay_3p_tp no_duplex_parallel 2


CountThreePrimeHighThroughputLibraryVariants2 :
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py miR-1 twist_reporter_assay_3p_2_tp duplex_series 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a twist_reporter_assay_3p_2_tp duplex_series 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a-21nt twist_reporter_assay_3p_2_tp duplex_series 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py miR-1 twist_reporter_assay_3p_2_tp no_duplex_series 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a twist_reporter_assay_3p_2_tp no_duplex_series 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a-21nt twist_reporter_assay_3p_2_tp no_duplex_series 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py miR-1 twist_reporter_assay_3p_2_tp duplex_parallel 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a twist_reporter_assay_3p_2_tp duplex_parallel 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a-21nt twist_reporter_assay_3p_2_tp duplex_parallel 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py miR-1 twist_reporter_assay_3p_2_tp no_duplex_parallel 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a twist_reporter_assay_3p_2_tp no_duplex_parallel 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a-21nt twist_reporter_assay_3p_2_tp no_duplex_parallel 1
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py miR-1 twist_reporter_assay_3p_2_tp duplex_series 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a twist_reporter_assay_3p_2_tp duplex_series 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a-21nt twist_reporter_assay_3p_2_tp duplex_series 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py miR-1 twist_reporter_assay_3p_2_tp no_duplex_series 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a twist_reporter_assay_3p_2_tp no_duplex_series 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a-21nt twist_reporter_assay_3p_2_tp no_duplex_series 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py miR-1 twist_reporter_assay_3p_2_tp duplex_parallel 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a twist_reporter_assay_3p_2_tp duplex_parallel 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a-21nt twist_reporter_assay_3p_2_tp duplex_parallel 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py miR-1 twist_reporter_assay_3p_2_tp no_duplex_parallel 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a twist_reporter_assay_3p_2_tp no_duplex_parallel 2
	bsub -q 18 -n 20 -R span[hosts=1] python ThreePrimePaperLet7ReporterAssay/CountReporterVariants3p.py let-7a-21nt twist_reporter_assay_3p_2_tp no_duplex_parallel 2




# AssignKmersEfficientSingle11 :
# 	@(if [ "$(mirna)" = "miR-7-23nt" ]; then \
# 		EXP="equilibrium2_nb"; \
# 		COND="I 12.6 40 0"; \
# 		MIRNA_NAME="miR-7"; \
# 	elif [ "$(mirna)" = "miR-124" ]; then \
# 		EXP="equilibrium_2_tp"; \
# 		COND="I I_combined 4 40 0"; \
# 		MIRNA_NAME="$(mirna)"; \
# 	else \
# 		EXP="equilibrium"; \
# 		COND="I I_combined 4 40 0"; \
# 		MIRNA_NAME="$(mirna)"; \
# 	fi; \
# 	echo $$MIRNA_NAME; \
# 	N_X=`wc -l $(HOME)$(DIR_as)site_categories/papercutoff/sites."$$MIRNA_NAME"_papercutoff.txt | cut -d " " -f 1`; \
# 	echo $$N_X; \
# # 	echo `expr "$$N_X"` \
# 	N_X=$$(( $$N_X + 1 )); \
# 	echo $$N_X; \
# 	for CON in $$COND; do \
# 		job="bsub -q 18 -n 20 python $(HOME)$(DIR_as)$(SCR_ak) $(mirna) $$EXP"; \
# 		job=$$job" $$CON $(n_constant) 11 -n_ex $$N_X"$(BUFFER3P)" -jobs 19"; \
# 		echo $$job; \
# 		$$job; \
# 	done; \
# 	)






# AssignKmers :
# 	@(job=""; \
# 	for CON in $(COND); do \
# 		job=$$job"python $(HOME)$(DIR_as)$(SCR_ak) $(mirna) $(exp) $$CON "; \
# 		job=$$job"$(n_constant) $(len_k)"$(N_EX)$(BUFFER3P)"\n"; \
# 	done; \
# 	echo $$job; \
# 	python $(HOME)general/SubmitJob.py "$$job")

# AssignKmersAll :
# 	if [ $(exp) = "kinetics" ]; then \
# 		for MIRNA in $(MIRNAk); do \
# 			line="make mirna=$$MIRNA exp=$(exp) n_constant=$(n_constant) "\
# 			line=$$line"len_k=$(len_k)"$(N_EX_ALL)" AssignKmers"; \
# 			echo $$line; \
# 			bsub $$line; \
# 		done; \
# 	fi; \

	# @(job=""; \
	# for CON in $(COND); do \
	# 	job=$$job"python $(HOME)$(DIR_as)$(SCR_ak) $(mirna) $(exp) $$CON "; \
	# 	job=$$job"$(n_constant) $(len_k)"$(N_EX)$(BUFFER3P)"\n"; \
	# done; \
	# echo $$job; \
	# python $(HOME)general/SubmitJob.py "$$job")


# AssignKmersFinal :
# 	@(job=""; \
# 	for MIRNA in $(MIRNAS); do \
# 		for LEN_K in $$(seq 8 10); do \
# 			for CON in $(COND); do \
# 				job=$$job"python $(HOME)$(DIR_as)$(SCR_ak) $$MIRNA $(exp)"; \
# 				job=$$job" $$CON $(n_constant) $$LEN_K \n"; \
# 			done; \
# 		done; \
# 	done; \
# 	python $(HOME)general/SubmitJob.py "$$job")

# AssignKmersBack :
# 	@(job=""; \
# 	for MIRNA in $(MIRNAS); do \
# 		if [ "$$MIRNA" = "miR-7-23nt" ]; then \
# 			COND='I 12.6'; \
# 			FILE="$(HOME)$(DIR_as)site_categories/paper/sites.miR-7_paper.txt"; \
# 		else \
# 			COND='I 4'; \
# 			FILE="$(HOME)$(DIR_as)site_categories/paper/sites.$(mirna)_paper.txt"; \
# 		fi; \
# 		NUMSITES=$$(($$(wc -l $$FILE | cut -d ' ' -f 1) + 1)); \
# 		for NUM in $$(seq 1 $$NUMSITES); do \
# 			for CON in $$COND; do \
# 				job=$$job"python $(HOME)$(DIR_as)$(SCR_ak) $$MIRNA $(exp)"; \
# 				job=$$job" $$CON $(n_constant) 8 -n_ex $$NUM\n"; \
# 			done; \
# 		done; \
# 	done; \
# 	echo $$job; \
# 	python $(HOME)general/SubmitJob.py "$$job")

# AssignKmersBack10 :
# 	@(job=""; \
# 	for MIRNA in $(MIRNAS); do \
# 		if [ "$$MIRNA" = "miR-7-23nt" ]; then \
# 			COND='I 12.6'; \
# 			FILE="$(HOME)$(DIR_as)site_categories/paper/sites.miR-7_paper.txt"; \
# 		else \
# 			COND='I 4'; \
# 			FILE="$(HOME)$(DIR_as)site_categories/paper/sites.$(mirna)_paper.txt"; \
# 		fi; \
# 		NUMSITES=$$(($$(wc -l $$FILE | cut -d ' ' -f 1) + 1)); \
# 		for NUM in $$(seq 1 $$NUMSITES); do \
# 			for CON in $$COND; do \
# 				job=$$job"python $(HOME)$(DIR_as)$(SCR_ak) $$MIRNA $(exp)"; \
# 				job=$$job" $$CON $(n_constant) 10 -n_ex $$NUM\n"; \
# 			done; \
# 		done; \
# 	done; \
# 	echo $$job; \
# 	python $(HOME)general/SubmitJob.py "$$job")



# AssignKmersQuick :
# 	job=""; \
# 	for MIRNA in $(MIRNAe); do \
# 		echo $$MIRNA; \
# 		if [ "$$MIRNA" = "miR-7-23nt" ]; then \
# 			COND='I 40'; \
# 			EXP='equilibrium2_nb'; \
# 		else \
# 			COND='I 40'; \
# 			EXP='equilibrium'; \
# 		fi; \
# 		echo $$COND; \
# 		for CON in $$COND; do \
# 			job=$$job"python $(HOME)$(DIR_as)$(SCR_ak) $$MIRNA $$EXP"; \
# 			job=$$job" $$CON -3 8 -n_ex $(n_ex)\n"; \
# 			job=$$job"python $(HOME)$(DIR_as)$(SCR_ak) $$MIRNA $$EXP"; \
# 			job=$$job" $$CON -3 10 -n_ex $(n_ex)\n"; \
# 		done; \
# 	done; \
# 	echo $$job; \
# 	python $(HOME)general/SubmitJob.py "$$job"

# AssignKmersQuickSingle :
# 	job=""; \
# 	if [ "$(mirna)" = "miR-7-23nt" ]; then \
# 		COND='I 40'; \
# 		EXP='equilibrium2_nb'; \
# 	else \
# 		COND='I 40'; \
# 		EXP='equilibrium'; \
# 	fi; \
# 	echo $$COND; \
# 	for CON in $$COND; do \
# 		job=$$job"python $(HOME)$(DIR_as)$(SCR_ak) $(mirna) $$EXP"; \
# 		job=$$job" $$CON -3 8 -n_ex $(n_ex)\n"; \
# 		job=$$job"python $(HOME)$(DIR_as)$(SCR_ak) $(mirna) $$EXP"; \
# 		job=$$job" $$CON -3 10 -n_ex $(n_ex)\n"; \
# 	done; \
# echo $$job; \
# python $(HOME)general/SubmitJob.py "$$job"

# AssignKmersmiR-1Comparison :
# 	job=""; \
# 	for CON in I 40; do \
# 		job=$$job"python $(HOME)$(DIR_as)$(SCR_ak) miR-1 equilibrium"; \
# 		job=$$job" $$CON $(n_constant) 9 -n_ex 0"$(BUFFER3P)"\n"; \
# 	done; \
# 	for CON in I L100A10; do \
# 		job=$$job"python $(HOME)$(DIR_as)$(SCR_ak) miR-1 equil_pilot"; \
# 		job=$$job" $$CON $(n_constant) 9 -n_ex 0"$(BUFFER3P)"\n"; \
# 	done; \
# echo $$job; \
# python $(HOME)general/SubmitJob.py "$$job"

AssignKmersKinetics :
	@(job=""; \
	for CON in $(COND); do \
		job=$$job"python $(HOME)$(DIR_as)$(SCR_ak) $(mirna) kinetics"; \
		job=$$job" $$CON $(n_constant)\n"; \
	done; \
	echo $$job; \
	python $(HOME)general/SubmitJob.py "$$job"; \
	)

AssignPositionalKmers :
	@(job=""; \
	for CON in $(COND); do \
		job=$$job"python $(HOME)$(DIR_as)$(SCR_apk) $(mirna) $(exp) $$CON "; \
		job=$$job"$(n_constant) $(len_k)"$(N_EX)$(BUFFER3P)"\n"; \
	done; \
	echo $$job; \
	python $(HOME)general/SubmitJob.py "$$job")


AssignPositionalKmersProgrammedLib :
	@(job_prefix="bsub -q 18 -n 20 -R span[hosts=1] python "; \
	job_prefix=$$job_prefix"$(HOME)$(DIR_as)$(SCR_akpl) $(mirna) $(exp) "; \
	for CON in $(COND); do \
		job=$$job_prefix" $$CON $(n_constant) $(len_k)"; \
		echo $$job; \
		$$job; \
	done)

AssignPositionalKmersProgrammedLibSeedEx :
	@(job_prefix="bsub -q 18 -n 20 -R span[hosts=1] python "; \
	job_prefix=$$job_prefix"$(HOME)$(DIR_as)$(SCR_akpl) $(mirna) $(exp) "; \
	for CON in $(COND); do \
		job=$$job_prefix" $$CON $(n_constant) $(len_k) -seedex"; \
		echo $$job; \
		$$job; \
	done)



FitFlankKds :
	@(job_py="python $(HOME)$(DIR_kd)$(SCR_fc) $(mirna) $(exp) $(n_constant) $(sitelist)"$(BUFFER3P); \
	echo $$job_py; \
	$$job_py; \
	if [ "$(mirna)" = "miR-7-23nt" ]; then \
		MIRSITE="miR-7"; \
	else \
		MIRSITE="$(mirna)"; \
	fi; \
	echo $$MIRSITE; \
	DIR="$(HOME)$(DIR_as)site_categories/$(sitelist)/sites."$$MIRSITE"_$(sitelist).txt"; \
	echo $$DIR; \
	for SITE in `cat $$DIR`; do \
		echo $$SITE; \
		job1="bsub Rscript $(HOME)$(DIR_kd)$(SCR_fkd) $(mirna) $(exp) $(n_constant) $(sitelist) $$SITE"$(BUFFER3P); \
		job2="bsub Rscript $(HOME)$(DIR_kd)$(SCR_fkd) $(mirna) $(exp) $(n_constant) $(sitelist) $$SITE -nocombI"$(BUFFER3P); \
		echo $$job1; \
		echo $$job2; \
		$$job1; \
		$$job2; \
	done)
		# job1="Rscript $(HOME)$(DIR_kd)$(SCR_fkd) $(mirna) $(exp)"; \
		# job1=$$job1" $(n_constant) $(sitelist) '$$SITE'"; \
		# echo $$job1; \
		# bsub $$job1; \
		# job2=$$job1" -nocombI"; \
		# echo $$job2; \
		# bsub $$job2; \

AssignAllSitesEquilibrium :
	echo $(MIRNAe)
	@(for MIRNA in $(MIRNAe); \
		do { make mirna=$$MIRNA exp=equilibrium n_constant=$(n_constant)\
			sitelist=$(sitelist) mir_start=$(mir_start) AssignSites;} \
	done)

AssignAll12mersEquilibrium :
	echo $(MIRNAe)
	@(for MIRNA in $(MIRNAe); \
		do { make mirna=$$MIRNA exp=equilibrium n_constant=$(n_constant)\
			sitelist=$(sitelist) mir_start=$(mir_start) Assign12mers;} \
	done)


AssignAllSitesEquilibriumMiR7 :
	echo $(MIRNAeMiR7)
	@(for MIRNA in $(MIRNAeMiR7); \
		do { make mirna=$$MIRNA exp=equilibrium_nb n_constant=$(n_constant)\
			sitelist=$(sitelist) AssignSites;} \
	done)

AssignAllSitesEquilibriumMiR7_2 :
	echo $(MIRNAeMiR7_2)
	@(for MIRNA in $(MIRNAeMiR7_2); \
		do { make mirna=$$MIRNA exp=equilibrium2_nb n_constant=$(n_constant)\
			sitelist=$(sitelist) AssignSites;} \
	done)

AssignAll12mersEquilibriumMiR7 :
	echo $(MIRNAeMiR7)
	@(for MIRNA in $(MIRNAeMiR7); \
		do { make mirna=$$MIRNA exp=equilibrium_nb n_constant=$(n_constant)\
			sitelist=12mers mir_start=$(mir_start) Assign12mers;} \
	done)


AssignAll12mersEquilibriumMiR7_2 :
	echo $(MIRNAeMiR7_2)
	@(for MIRNA in $(MIRNAeMiR7_2); \
		do { make mirna=$$MIRNA exp=equilibrium2_nb n_constant=$(n_constant)\
			sitelist=12mers mir_start=$(mir_start) Assign12mers;} \
	done)

AssignAllFlanksEquilibrium :
	echo $(MIRNAe)
	@(for MIRNA in $(MIRNAe); \
		do { make mirna=$$MIRNA exp=equilibrium n_constant=$(n_constant)\
			sitelist=$(sitelist) AssignFlanks;} \
	done)

CalculatePlfold_miRNA:
	@(SITES=`cat $(HOME)$(DIR_as)sites.$(mirna)_$(sitelist).txt`; \
	for SITE in $$SITES; \
		do { echo bsub -n 20 python $(HOME)$(DIR_plfold)$(SCR_plfold) $(mirna) $(exp) $(cond) \
			$(n_constant) $(sitelist) $$SITE; \
			bsub -n 21 python $(HOME)$(DIR_plfold)$(SCR_plfold) $(mirna) $(exp) $(cond) \
			$(n_constant) $(sitelist) $$SITE;} \
	done)

AssignSitesAllKinetics :
	echo $(MIRNAk)
	@(for MIRNA in $(MIRNAk); \
		do { make mirna=$$MIRNA exp=kinetics n_constant=$(n_constant)\
			sitelist=$(sitelist) AssignSites;} \
	done)

AssignFlanksKinetics :
	@(DIR="/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AssignSiteTypes/"; \
	SCR="MakeFlankingNucleotideReadFiles_PAPER.py"; \
	for CON in I I_TGT_combined I_ACA_combined 2,1 2,2 5,1 5,2 8,1 8,2 15,1 \
		15,2 30 30,1 30,2 90 300 900 2700 7200 7200,1 7200,2 0 0- Equil Equil-; \
		do { bsub -n 21 python $$DIR$$SCR $(mirna) kinetics $$CON\
			$(n_constant) $(sitelist);} \
	done)

AssignFlanksKineticsAll :
	for MIRNA in miR-1 let-7a miR-124 lsy-6; \
		do { make mirna=$$MIRNA n_constant=$(n_constant) sitelist=$(sitelist)\
			AssignFlanksKinetics;} \
	done

FitSiteKds :
	python $(HOME)$(DIR_kd)$(SCR_sc) $(mirna) $(exp) $(n_constant)\
		$(sitelist)
	bsub Rscript $(HOME)$(DIR_kd)$(SCR_skd) $(mirna) $(exp) $(n_constant)\
		$(sitelist) $(reps)

Fit12merKds :
	@(for MIR_START in $$(seq 1 5); \
		do { python $(HOME)$(DIR_kd)$(SCR_sc) $(mirna) $(exp) $(n_constant)\
				12mers -mir_start $$MIR_START; \
			bsub Rscript $(HOME)$(DIR_kd)$(SCR_12merkd) $(mirna) $(exp) $(n_constant)\
				12mers -mir_start $$MIR_START;} \
	done)

Fit16merKds :
	@(for MIR_START in $$(seq 1 5); \
		do { python $(HOME)$(DIR_kd)$(SCR_sc) $(mirna) $(exp) $(n_constant)\
				16mers -mir_start $$MIR_START -split16 left; \
			python $(HOME)$(DIR_kd)$(SCR_sc) $(mirna) $(exp) $(n_constant)\
				16mers -mir_start $$MIR_START -split16 right; \
			bsub Rscript $(HOME)$(DIR_kd)$(SCR_12merkd) $(mirna) $(exp) $(n_constant)\
				16mers -mir_start $$MIR_START -split16 left;} \
			bsub Rscript $(HOME)$(DIR_kd)$(SCR_12merkd) $(mirna) $(exp) $(n_constant)\
				16mers -mir_start $$MIR_START -split right;} \
	done)

Fit12merKdsNoComb :
	@(for MIR_START in $$(seq 2 5); \
		do { python $(HOME)$(DIR_kd)$(SCR_sc) $(mirna) $(exp) $(n_constant)\
				12mers -mir_start $$MIR_START $(reps); \
			bsub Rscript $(HOME)$(DIR_kd)$(SCR_12merkd) $(mirna) $(exp) $(n_constant)\
				12mers -nocombI -mir_start $$MIR_START $(reps);} \
	done)

FitAllSiteKds :
	@(for MIRNA in $(MIRNAe); \
		do { make mirna=$$MIRNA exp=$(exp) n_constant=$(n_constant)\
			sitelist=$(sitelist) FitSiteKds;} \
	done)

FitAll12merKds :
	for MIRNA in $(MIRNAe); \
		do { make mirna=$$MIRNA exp=$(exp) n_constant=$(n_constant)\
				Fit12merKds;} \
		done;} \
	done

FitAll12merKdsNoComb :
	for MIRNA in $(MIRNAe); \
		do { for MIR_START in $$(seq 2 5); \
			do { make mirna=$$MIRNA exp=$(exp) n_constant=$(n_constant)\
				mir_start=$$MIR_START Fit12merKdsNoComb;} \
		done;} \
	done

FitAllSiteKdsMiR7 :
	@(for MIRNA in $(MIRNAeMiR7_2); \
		do { make mirna=$$MIRNA exp=equilibrium2_nb n_constant=$(n_constant) sitelist=$(sitelist)\
			FitSiteKds;} \
	done)

FitAll12merKdsMiR7 :
	@(for MIRNA in $(MIRNAeMiR7_2); \
		do { make mirna=$$MIRNA exp=equilibrium2_nb n_constant=$(n_constant)\
			Fit12merKds;} \
	done)

FitAll12merKdsMiR7NoComb :
	@(for MIRNA in $(MIRNAeMiR7_2); \
		do { make mirna=$$MIRNA exp=equilibrium2_nb n_constant=$(n_constant)\
			Fit12merKdsNoComb;} \
	done)


FitAllFlankKds :
	@(for MIRNA in miR-1 let-7a miR-124 lsy-6; \
		do { make mirna=$$MIRNA exp=$(exp) n_constant=$(n_constant)\
		sitelist=$(sitelist) FitFlankKds;} \
	done)

MakeRepressionFiles :
	@(for MIRNA in miR-1 let-7a miR-124 lsy-6; \
		do { python $(DIR_kd)/$(SCR_kdrepfile) $$MIRNA $(exp)\
		$(n_constant) $(sitelist);} \
	done)

FitAllSiteKdsOld :
	make mirna=let-7a  exp=equilibrium     n_constant=$(n_constant) sitelist=$(sitelist) reps=$(reps) FitSiteKds
	make mirna=let-7a  exp=equilibrium_nb  n_constant=$(n_constant) sitelist=$(sitelist) reps=$(reps) FitSiteKds
	# make mirna=let-7a  exp=equil_mmseed_nb n_constant=$(n_constant) sitelist=$(sitelist) reps=$(reps) FitSiteKds
	# make mirna=let-7a  exp=equil_seed_nb   n_constant=$(n_constant) sitelist=$(sitelist) reps=$(reps) FitSiteKds
	make mirna=miR-1   exp=equilibrium     n_constant=$(n_constant) sitelist=$(sitelist) reps=$(reps) FitSiteKds
	make mirna=miR-155 exp=equilibrium     n_constant=$(n_constant) sitelist=$(sitelist) reps=$(reps) FitSiteKds
	make mirna=miR-124 exp=equilibrium     n_constant=$(n_constant) sitelist=$(sitelist) reps=$(reps) FitSiteKds
	make mirna=lsy-6   exp=equilibrium     n_constant=$(n_constant) sitelist=$(sitelist) reps=$(reps) FitSiteKds

FitSiteKdsOld :
	python /lab/bartel1_ata/mcgeary/computation/AgoRBNS/SolveForKds/GenerateSiteKds_PAPER.py $(mirna) $(exp) $(n_constant) $(sitelist)
	bsub Rscript /lab/bartel1_ata/mcgeary/computation/AgoRBNS/SolveForKds/GenerateSiteKds_PAPER.R $(mirna) $(exp) $(n_constant) $(sitelist) $(reps)


