########use to get sitemap

import imp # Used to import general.py
import time # Used for time
imp.load_source("general",("/lab/bartel1_ata/nbisaria/AgoRBNS/general/general.py"))
imp.load_source("sitetypes",("/lab/bartel1_ata/nbisaria/AgoRBNS/AssignSiteTypes/sitetypes.py"))
from general import parse_arguments, print_time_elapsed, get_analysis_path, multiprocess_file, readline, readpicklefile, savelisttopickle, getseedlist, get_mirnaseq, get_rc
from sitetypes import get_seq_site_map_all, assign_site_type_to_read_streamline,get_seq_site_map_seedMMs,assign_site_type_to_read_regex
import pickle
from itertools import product
import pandas as pd

def makesitemaps(mirna, library_type, mm3p=False):


    filename = '/lab/bartel1_ata/nbisaria/AgoRBNS/AssignSiteTypes/Uniquesites_4to12_' + mirna + '.txt'
    sites = pd.read_csv(filename, sep='\t')
    print(sites)
    return
      
    # MMorbulgesites = (sites['length'] > 6) & (sites['length'] < 9) & ((sites['lesiontype'] == 'B') | (sites['lesiontype'] == 'MM')) & (sites['start']>0) & (sites['stop']<9)

    maxlen = len(get_mirnaseq(mirna))
    if (library_type != 'mix3') and (library_type != 'mix2_3') and (library_type != 'mix2_3_2') :
        MMorbulgesites = (sites['length'] > 7) & (sites['length'] < 9) & (sites['lesiontype'] == 'MM') & (sites['start']>0) & (sites['stop']<9)
        compsites = (sites['start']>0) & (sites['stop']<9) & (sites['length'] > 5)  & (sites['length'] < 10) & (sites['lesiontype'] == 'None')
        seedsites = sites[compsites | MMorbulgesites]
        seedsites = seedsites[~((seedsites['length'] == 8) &  ((seedsites['lesionloc'] == '0') | (seedsites['lesionloc'] == '7')) | (seedsites['length'] == 7) &  ((seedsites['lesionloc'] == '0') | (seedsites['lesionloc'] == '6')))]
        
    
    else:
        MMorbulgesites = (sites['length'] > 7) & (sites['length'] < 9) & (sites['lesiontype'] == 'MM') & (sites['start']>0) & (sites['stop']<9) | (sites['lesiontype'] == '2MM')
        # compsites = (sites['start']>0) & (sites['stop']<9) & ((sites['length'] == 6) | (sites['length'] == 8))  & (sites['length'] < 10) & (sites['lesiontype'] == 'None')
        compsites = (sites['start']>0) & (sites['stop']<9) & ((sites['length'] == 8))  & (sites['lesiontype'] == 'None')
        seedsites = sites[compsites | MMorbulgesites]


    if mm3p:
        # includes 7mers that have mismatches
        threepCsites =  ((sites['start']>8) & (sites['stop']<maxlen+1) & (sites['length'] == 7)  & (sites['lesiontype'] == 'None')) | ((sites['start']>8) & (sites['stop']<maxlen+1) & (sites['length'] ==7)  & (sites['lesiontype'] == 'MM') & ((sites['lesionloc'] == '0') | (sites['lesionloc'] == '6')))
        
    
    else:
        threepCsites =  (sites['start']>8) & (sites['stop']<maxlen+1) & (sites['length'] > 4)  & (sites['length'] < 11) & (sites['lesiontype'] == 'None')
    
    # if mirna[0:6] == 'let-7a':
    #     # extrasites = (sites['start']>8) & (sites['stop']<22) & (sites['length'] > 4)  & (sites['length'] < 10) & ((sites['targetseq'].str.contains('AGCC')) |(sites['targetseq'].str.contains('TCAAC'))|(sites['targetseq'].str.contains('TAACC'))|(sites['targetseq'].str.contains('TCAACA'))|(sites['targetseq'].str.contains('AACAACC')))
    #     extrasites = (sites['start']>8) & (sites['stop']<maxlen+1) & (sites['length'] > 4)  & (sites['length'] < 10) & ((sites['targetseq'].str.contains('AGCC')) |(sites['targetseq'].str.contains('TCAAC'))|(sites['targetseq'].str.contains('TAACC')))
        
    #     threepsites = sites[threepCsites | extrasites]
    # else:
    threepsites = sites[threepCsites]


    fileout = '/lab/bartel1_ata/nbisaria/AgoRBNS/AssignSiteTypes/sitelists/pd_seeds_exhaustive_registers.txt'
    with open(fileout, 'wb') as f:
        pickle.dump(seedsites, f)

    fileout = '/lab/bartel1_ata/nbisaria/AgoRBNS/AssignSiteTypes/sitelists/pd_threeps_exhaustive_registers.txt'
    with open(fileout, 'wb') as f:
        pickle.dump(threepsites, f)

   
    sitesseed = seedsites[['name', 'targetseq']]
    sites3p = threepsites[['name', 'targetseq']]

    seq_site_map_seed = sitesseed.set_index('name').T.to_dict('index')['targetseq']
    seq_site_map_3p = sites3p.set_index('name').T.to_dict('index')['targetseq']
    
    # sites.set_index("name", drop=True, inplace=True)
    return seq_site_map_seed, seq_site_map_3p




def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA","experient","condition", "library_type"]
    mirna, experiment, condition, library_type, replicate = parse_arguments(arguments)
    mm3p = False
    
    if mm3p:
        stradd = '_mms'
    else:
        stradd = ''
    seq_site_map_seed, seq_site_map_3p = makesitemaps(mirna, library_type, mm3p)
    print(seq_site_map_3p)
    fileout = '/lab/bartel1_ata/nbisaria/AgoRBNS/AssignSiteTypes/sitelists/map_seeds_%s_%s%s.txt' %(mirna, library_type,stradd)
    # savelisttopickle(seq_site_map_seed, fileout)
    finaldict = seq_site_map_seed.copy()
    finaldict.update(seq_site_map_3p)
    fileout = '/lab/bartel1_ata/nbisaria/AgoRBNS/AssignSiteTypes/sitelists/map_streamline_seedsand3psites_%s_%s%s.txt' %(mirna, library_type,stradd)
    # savelisttopickle(finaldict,fileout)

    site_seq_map_seed = {value: key for key, value in seq_site_map_seed.items()}
    site_seq_map_3p = {value: value for key, value in seq_site_map_3p.items()}
  
    seeds = getseedlist(library_type, mirna)
    mirnaseed = get_rc(get_mirnaseq(mirna)[0:8])

    counts_final = {}
    multi_counts_final = {}
    for seed in seeds:
        final3pseqs = []
        lengths = range(0,18)
        seqs = seq_site_map_3p.values()
        site_seq_map_3p = {}
        for length in lengths:
            for seq in seqs:
                finalseq = seq + '-N' + str(length) + '-' + seed
                final3pseqs.append(finalseq)

if __name__ == "__main__":
    main()
