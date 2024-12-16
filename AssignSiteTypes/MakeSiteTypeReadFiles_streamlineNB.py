################################################################################
#MakeSiteTypeReadFiles.py
################################################################################
import imp # Used to import general.py
import time # Used for time
imp.load_source("general",("/lab/bartel1_ata/nbisaria/AgoRBNS/general/general.py"))
imp.load_source("sitetypes",("/lab/bartel1_ata/nbisaria/AgoRBNS/AssignSiteTypes/sitetypes.py"))
from general import parse_arguments, print_time_elapsed, get_analysis_path, multiprocess_file, readline, readpicklefile, savelisttopickle, getseedlist, get_mirnaseq, get_rc
from sitetypes import get_seq_site_map_all, assign_site_type_to_read_streamline,get_seq_site_map_seedMMs,assign_site_type_to_read_regex
import pickle
from itertools import product
import pandas as pd
# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.
#here we have two site maps, one that is for normal sites that can be found in the random region
#the other is for sites that are found anywhere but require a regular expression search

##will need to edit to make more general for miR-1 and other miRNAs
def makesitemaps(mirna,library_type, mm3p=False):

    # if (mirna == 'let-7a_plus1') | (mirna == 'let-7a_minus1') | (mirna == 'let-7a_miR-155_chimera'):
    #     mirna = 'let-7a'
    # elif (mirna=='miR-155_let-7a_chimera'):
    #     mirna = 'miR-155' 
    filename = '/lab/bartel1_ata/nbisaria/AgoRBNS/AssignSiteTypes/Uniquesites_4to12_' + mirna + '.txt'
    sites = pd.read_csv(filename, sep='\t')
    
    
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


    # want sites that are atleast a 7 mer if you are including mismacthes
    #don't including the sites that are flanks of 7mer and 8mer
    # if mirna == 'miR-155':
    #     threepCsites =  (sites['start']>8) & (sites['stop']<25) & (sites['length'] > 4)  & (sites['length'] < 10) & (sites['lesiontype'] == 'None')
    # else:
    
    if mm3p:
        # includes 7mers that have mismatches
        #old version will try again
        # threepCsites =  ((sites['start']>8) & (sites['stop']<maxlen+1) & (sites['length'] > 4)  & (sites['length'] < 10) & (sites['lesiontype'] == 'None')) | ((sites['start']>8) & (sites['stop']<maxlen+1) & (sites['length'] ==7)  & (sites['lesiontype'] == 'MM') & ((sites['lesionloc'] == '0') | (sites['lesionloc'] == '6')))
        # threepCsites =  ((sites['start']>8) & (sites['stop']<maxlen+1) & (sites['length'] ==7)  & (sites['lesiontype'] == 'MM') & ((sites['lesionloc'] == '0') | (sites['lesionloc'] == '6')))
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

def getSiteRank():
    #get the siteranking of enrichment values 
    seq_site_map = makeseedmap()
    df = makeregexsites_all()[['name', 'seq']]
    allsites = seq_site_map.values() + df.name.values.tolist()
    
# def getseedlist(library_type, ):
#     seed = 'CTACCTCA'
#     express = []
#     mapMMs = {'A':['T','G','C'], 'T':['A','G','C'], 'G':['A','T','C'], 'C':['A', 'T', 'G']}
#     if library_type == 'mix2':
#         start = 1
#         stop = -1
#     else:
#         start = 0
#         stop = 0

#     for i in range(start, len(seed)+stop):
#         seed_list = list(seed)
#         mms = mapMMs[seed_list[i]]
#         # print(i)
#         for mm in mms:
#             seed_list[i] = mm
#             # print(mm)
#             seed_ = ''.join(seed_list)
#             express.append(seed_)
#     # express.append(seed)
#     l1 = express
        
#     fileout = '/lab/bartel1_ata/nbisaria/AgoRBNS/AssignSiteTypes/sitelists/seedlist.txt'
#     savelisttopickle(l1, fileout)
#     return l1


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
    savelisttopickle(seq_site_map_seed, fileout)
    finaldict = seq_site_map_seed.copy()
    finaldict.update(seq_site_map_3p)
    fileout = '/lab/bartel1_ata/nbisaria/AgoRBNS/AssignSiteTypes/sitelists/map_streamline_seedsand3psites_%s_%s%s.txt' %(mirna, library_type,stradd)
    savelisttopickle(finaldict,fileout)




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



        site_seq_map_3p = {value: value for key, value in seq_site_map_3p.items()}
        #first locate reads
        reads_path = get_analysis_path(mirna,experiment,condition,library_type,"reads", replicate)
        print(reads_path)
        reads_path = '/'.join(reads_path.split('/')[:-1]) + '/byseeds/' + seed + '_' + condition + '.txt'
        print(reads_path)
        site_counts_path = get_analysis_path(mirna,experiment,condition,library_type,"site_counts_streamline3p_longloop" + stradd, replicate)
        site_counts_path = '/'.join(site_counts_path.split('/')[:-1]) + '/' + seed + '_' + condition + '.txt'

        with open(reads_path,"rb") as file_in:
            results = multiprocess_file(
                file_in, readline, 10000000, assign_site_type_to_read_streamline, site_seq_map_seed,site_seq_map_3p,final3pseqs, seed, mirnaseed)

        dict_threads = [i[0] for i in results]
        dict_threads_multi = [i[1] for i in results]
       

        keys = set([j for i in dict_threads for j in i.keys()])
        values = [sum([thread[key] for thread in dict_threads if key in thread])
                  for key in keys]
        counts_sites_map = {key : value for (key,value) in zip(keys,values)}
        counts_final = { k: counts_final.get(k, 0) + counts_sites_map.get(k, 0) for k in set(counts_final) | set(counts_sites_map) }

        with open(site_counts_path,"wb") as file_out:
            for key in sorted(keys):
                value = counts_sites_map[key]
                file_out.write("%s:\t%s\n" % (key, value))

        with open(site_counts_path[:-4] + '_dict.txt', 'wb') as f:
            pickle.dump(counts_sites_map, f)

        keys = set([j for i in dict_threads_multi for j in i.keys()])
        values = [sum([thread[key] for thread in dict_threads_multi if key in thread])
                  for key in keys]
        multi_counts = {key : value for (key,value) in zip(keys,values)}
        multi_counts_final = { k: multi_counts_final.get(k, 0) + multi_counts.get(k, 0) for k in set(multi_counts_final) | set(multi_counts) }

        with open(site_counts_path[:-4] + '_multi.txt',"wb") as file_out:
            for key in sorted(keys):
                value = multi_counts[key]
                file_out.write("%s:\t%s\n" % (key, value))

        with open(site_counts_path[:-4] + '_multi_dict.txt', 'wb') as f:
            pickle.dump(multi_counts, f)

    site_counts_path = get_analysis_path(mirna,experiment,condition,library_type,"site_counts_streamline3p_longloop" + stradd, replicate)
    multi_counts_path = get_analysis_path(mirna,experiment,condition,library_type,"site_counts_streamline3p_longloop" + stradd + "_multi", replicate)

    keys = counts_final.keys()
    with open(site_counts_path[:-4] + '.txt',"wb") as file_out:
        for key in sorted(keys):
            value = counts_final[key]
            file_out.write("%s:\t%s\n" % (key, value))


    with open(site_counts_path[:-4] + '_dict.txt', 'wb') as f:
        pickle.dump(counts_final, f)

    keys = multi_counts_final.keys()
    with open(multi_counts_path,"wb") as file_out:
        for key in sorted(keys):
            value = multi_counts_final[key]
            file_out.write("%s:\t%s\n" % (key, value))


    with open(multi_counts_path[:-4] + '_dict.txt', 'wb') as f:
        pickle.dump(multi_counts_final, f)


    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

