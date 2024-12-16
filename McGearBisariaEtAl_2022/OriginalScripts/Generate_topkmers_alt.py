# This function takes concentrations, experiment_type, to generate
#the Kmer dict, enrichment Matrix, a heatmap of the kmers, and saves it to the kmer directory

# import imp
# imp.load_source('helpers', '/lab/bartel1_ata/nbisaria/AgoRBNS/kmer_analysis/find_all_kmer_positions/find_kmer_helpers.py')
# from helpers import generate_kmers
import matplotlib as mpl
mpl.use('Agg')
import imp
import numpy as np
import pandas as pd
import pickle
import itertools as it
import seaborn as sns
import matplotlib.pyplot as plt

import os
from scipy.spatial import distance
from scipy.cluster import hierarchy
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics import hamming_loss
import subprocess
from Bio.Seq import Seq
import csv
imp.load_source("plotfun",("/lab/bartel1_ata/nbisaria/AgoRBNS/general/plotfun.py"))
from plotfun import make_heatmap, clustermap

imp.load_source("general",
                ("/lab/bartel1_ata/nbisaria/"
                 "AgoRBNS/general/general.py"
                )
               )
from general import savelisttopickle,listtofasta,savelisttocsv,parse_arguments_general
# from skbio.sequence import DNASequence, RNASequence
#want to generate two kinds of data, which depends on enrichment values of kmers in general and kmers as a function of position. 

#will start by just stupidly adding all the values of the array and calculating an enrichment value for each kmer and outputing that 
#as the same kind of pandas matrix as the site-type


def getKmerdicts(kmer_path, conditions, extras='', extras2 = '', thresh=False, threshval = 5):
	os.chdir(kmer_path)
	alldicts={}
	for i, condition in enumerate(conditions):
		if len(extras2)>1:
			extras2_ = extras2[i]
		else:
			extras2_ = extras2
		file_in = condition +"_kmers_no6mer.txt"
		cmd = "df_%s = df" %("_".join(condition.split('.')))
		code = compile(cmd, "script", "exec")
		with open(file_in, 'rb') as f:
			kmer_dict = pickle.load(f, encoding="latin1")
			if thresh:
				df = pd.DataFrame(kmer_dict).clip(lower=threshval).replace(threshval,0)
			else:
				df = pd.DataFrame(kmer_dict)
			exec(cmd)
			dictname = 'df_%s' %("_".join(condition.split('.')))
			exec('alldicts[dictname]= df_%s' %("_".join(condition.split('.'))))
	return alldicts


def getEnrichMatrix(alldicts, conditions):
	TR= pd.DataFrame(0, columns = conditions, index = alldicts['df_I'].index)
	allenrichments = {}
	for condition in conditions:
		dictname = 'df_%s' %("_".join(condition.split('.')))
		exec("TR['%s'] = alldicts[dictname].sum(axis=1)" %(condition))
		print(TR)
		# enrichmatrix = pd.DataFrame(0, columns = alldicts['df_I'].columns, index = alldicts['df_I'].index)
		if condition != 'I':
			exec("allenrichments['%s'] = alldicts[dictname].divide(TR['%s'] , axis=0).divide(alldicts['df_I'].divide(TR['I'], axis=0))" %(condition, condition))
	return allenrichments

def allkmers_length(length):
	kmers = []
	kmers += ["".join(kmer) for kmer in list(it.product(["A","C","G","T"],repeat=length))]
	return kmers

def rm_nan(x):
	x = x[~np.isnan(x)]
	return x

# def make_heatmap(matrix, title, filename):
#     plt.figure(figsize = (16,16));
#     sns.heatmap(matrix)
#     plt.title(title)
#     plt.xlabel('Position in read')
#     plt.yticks(rotation=0)
#     plt.savefig(filename + '.pdf')
    # plt.show()

# def clustermap(matrix, title, filename):
# 	g = sns.clustermap(matrix,col_cluster=False,figsize=(16,16))
# 	#g = sns.clustermap(matrix,row_cluster=False)
# 	plt.title(title)
# 	plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
# 	plt.savefig('/lab/bartel1_ata/nbisaria/AgoRBNS/figures/'+filename + '.pdf')
# 	plt.show()


def revcomp(seq):
	comp = {'A': 'T', 'C':'G', 'G':'C', 'T':'A','N':'N', 'U':'A'}
	rev_com = "".join(comp.get(base,base) for base in reversed(seq))
	return rev_com

def makeClustmap(vals,title):
	matrix = vals.iloc[:,2:26].replace([np.nan, np.inf],0) 
	matrixnorm = matrix.div(matrix.sum(axis=1), axis=0)
	clustermap(matrix, title, title)
	clustermap(matrixnorm, title, title + '_norm')

def get_topkmers_pos(condition, thresh,allenrich, filtseeds=False, filtrandom = False):
	topvals = {}
	topkmers = []
	for i in allenrich[condition].index:
		topvals[str(i)] = allenrich[condition].iloc[i:i+1].T.sort_values([i], ascending=False).replace([np.inf, -np.inf], np.nan).dropna()
		test = topvals[str(i)] 
		topkmers.append(test[test>thresh].dropna().index.tolist())
	topkmers2 = [item for sublist in topkmers for item in sublist]
	setkmers = set(topkmers2)
	if filtseeds:
		substring1 = 'TACCT'
		# substring1 = 'ACCTC'
		substring2 = 'ACCTC'
		# substring2 = 'ACCTC'
		substring3 = 'CCTC'
		notseeds = [string for string in setkmers if substring2 not in string and substring1 not in string and substring3 not in string]
		idxs = pd.Index(notseeds)
		finalkmers = notseeds
	
	if filtrandom:
		# remove all kmers that correspond to 8mers enriched in the random library, probably need to reduce threshold and take top kmers at different concentraitons. 
		with open('top8mers.txt', 'r') as f:
			reader = csv.reader(f)
			top8mers = list(reader)[0]
		setkmers2 = list(setkmers)
		notvals = [val for val in setkmers2 if val not in top8mers]
		idxs = pd.Index(notvals)
		finalkmers = notvals
	else:
		idxs = pd.Index(setkmers)
		finalkmers = setkmers
	topvals_filt = pd.DataFrame(0,index = finalkmers, columns =allenrich[condition].index )
	for i in allenrich[condition].index:
		topvals_filt[i] = allenrich[condition].iloc[i,:][idxs]
	topvals_sorted = topvals_filt.replace([np.nan, np.inf], 0).sort_values([13], ascending=False)
	return topvals_sorted


def main():

	file_in1 = "/lab/solexa_bartel/nbisaria/analysis/let-7a/equilibrium/mix2/kmers/findupto8mers/I_kmers_aligned3.txt"
	file_in2 = "/lab/solexa_bartel/nbisaria/analysis/let-7a/equilibrium/mix2/kmers/findupto8mers/40_kmers_aligned3.txt"

	with open(file_in1, 'rb') as f1:
		kmer_dict1 = pickle.load(f1, encoding="latin1")
	with open(file_in2, "rb") as f2:
		kmer_dict2 = pickle.load(f2, encoding="latin1")

	print(len(kmer_dict1))
	print(len(kmer_dict2))

	# print(kmer_dict1["AACAA"])
	# print(kmer_dict1["TCACACTC"])

	# kmer_df_1 = pd.DataFrame.from_dict(kmer_dict1, orient="index")
	# kmer_df_2 = pd.DataFrame.from_dict(kmer_dict2, orient="index")

	# print(kmer_df_1)
	# print(kmer_df_2)

	all_dicts = {"df_I" : pd.DataFrame(kmer_dict1),
				 "df_40" : pd.DataFrame(kmer_dict2)}

	print(all_dicts)
	output = getEnrichMatrix(all_dicts, conditions=["I", "40"])["40"]

	top_kmers = {i : 0 for i in output.columns}
	print(top_kmers)
	kmers = list(top_kmers.keys())
	print(kmers[0])
	# kmer_df_1 = kmer_df_1/kmer_df_1.sum()
	# kmer_df_2 = kmer_df_2/kmer_df_2.sum()

	# print(kmer_df_1)
	# print(kmer_df_2)


	return




	experiment = 'equilibrium'
	arguments = ["miRNA", "library_type", "condition"]
	[mirna, library_type, condition] = parse_arguments_general(arguments)
	kmer_path = ("/lab/solexa_bartel/nbisaria/analysis/%s/equilibrium/%s/kmers/8mers_no6mers_pos_1_34/" % (mirna, library_type))
	print(kmer_path)
	# Eval = int(Eval)
	val = 5
	# condition = '4'
	condition = [condition]
	conditions  = ['I'] + condition
	alldicts = getKmerdicts(kmer_path,conditions,thresh=True, threshval=val)
	print(alldicts)
	return
	allenrich = getEnrichMatrix(alldicts, conditions)
	thresh = 10
	topvals = get_topkmers_pos(condition[0],thresh,allenrich)
	# topvals_filt_40_filt = get_topkmers_pos('40',thresh,allenrich_mix2_no6mer, filtseeds=True)
	# print('here2')
	# makeClustmap(topvals, mirna +'_' + condition[0] +'_' + str(thresh) + '_reads' + str(val) + 'no6mer')
	E_position = allenrich[condition[0]].T
	E_position = E_position.replace([np.inf, np.nan],0)

	# Eval =3 
	def topscore(row, numEvals):
	    val = row.sort_values(ascending=False)[0:numEvals].sum()
	    return val
	start = 3
	fin = 19
	# start = 3
	# fin = 17
	E_position['Score'] = E_position.apply(lambda row: topscore(row[start:fin],Eval), axis=1)
	E_position = E_position.sort_values(["Score"], ascending = False)
	make_heatmap(E_position.iloc[0:100,:-1], mirna, 'topp100kmers_conc' + str(condition[0]) + mirna + 'E' + str(Eval) + '_thresh' + str(val) + str(start) + "_" + str(fin), directory=kmer_path)
	clustermap(E_position.iloc[0:100,:-1],mirna,'cluster_topp100kmers_conc' + str(condition[0]) + mirna + 'E' + str(Eval) + '_thresh' + str(val) + str(start) + "_" + str(fin), directory=kmer_path )
	topkmers = E_position.iloc[0:100].index.tolist()
	listtofasta(topkmers, kmer_path + 'top100kmers_conc' + str(condition[0]) + str(start) +'_' +  str(fin) + '_' +'E' + str(Eval) + '_')
	savelisttocsv(topkmers, kmer_path + 'top100kmersrankedlist.txt')

if __name__ == "__main__":
    main()

