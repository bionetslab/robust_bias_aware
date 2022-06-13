import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import mygene
import math



# bait_usages=[]
# node_degrees=[]
# G = nx.read_graphml("/Users/surya/Documents/RECORDRESULTS_STUDYBIAS/robust-eval-main-1/amim_test_suite/algorithms/robust/ms.graphml")

# # bait_usages = list(nx.get_node_attributes(G, 'bait_usage').values())

# for node, attr in G.nodes(data=True):
#     bait_usages.append(G.nodes[node]['bias_data'])
#     node_degrees.append(G.degree[node])
    
# data = np.column_stack([node_degrees, bait_usages])
# datafile_path = "NodeDegree_VS_BaitUsage_a_APID.txt"
# np.savetxt(datafile_path , data, fmt=['%s', '%s'])

#########################################################################################################

# BAIT=pd.read_csv('/Users/surya/Desktop/StudyBiasPaper/CODE/Robust/robust-main/bias-aware data/bait_usage.txt', sep=' ')
# bait_uniprot=BAIT['bait_uniprot']
# bait_symbol=BAIT['bait_symbol']
# bait_usage=BAIT['bait_usage']

# bait_uniprot_list=bait_uniprot.to_list()
# bait_symbol_list=bait_symbol.to_list()
# bait_usage_list=bait_usage.to_list()

# ZIP=zip(bait_uniprot_list, bait_usage_list)
# DICT=dict(ZIP)

# LENGTH_orig=len(bait_uniprot_list)
# list_status=np.ones((LENGTH_orig,1))
# UniprotsForMissingGeneSymbols=[]

# bait_uniprot_list1=[]
# bait_symbol_list1=[]
# bait_usage_list1=[]
# for i in range(LENGTH_orig):
#     if not(isinstance(bait_symbol_list[i], str)):
#         list_status[i]=0
#         UniprotsForMissingGeneSymbols.append(bait_uniprot_list[i])
#     else:
#         bait_uniprot_list1.append(bait_uniprot_list[i])
#         bait_symbol_list1.append(bait_symbol_list[i])
#         bait_usage_list1.append(bait_usage_list[i])

# noOfMissingGeneSymbols=len(UniprotsForMissingGeneSymbols)

# data = np.column_stack([bait_uniprot_list1, bait_usage_list1])
# datafile_path = "uniprotProteinID_baitUsage.txt"
# np.savetxt(datafile_path , data, fmt=['%s', '%s'])

# ----------------------------------------------------------------------------------------------------

# # with open("abs_uniprots.txt", "w") as output:
# #     output.write(str(abs_uniprots))

# # ------- SYNGO (no) -------

# # # mapping=pd.read_excel('idmap_query,symbol.xlsx')
# # # query=list(mapping['query'])
# # # symbol=list(mapping['symbol'])
# # # zipper=zip(query,symbol)
# # # dictionary=dict(zipper)

# # ------- mygene (yes) -------

# noOfGenes=len(UniprotsForMissingGeneSymbols)
# gene_entrezIDs_withFoundGeneNames=[]
# EntrezIDs_Status=np.zeros((noOfGenes,1))
# statuscounter=0
# cntr=0
# NotFoundGenes=[]
# try:
#     mg = mygene.MyGeneInfo()
#     out = mg.querymany(UniprotsForMissingGeneSymbols, scopes= 'uniprot', fields='symbol', species='human', verbose=False)
#     gene_names = []
#     for line in out:
#         try:
#             EntrezIDs_Status[statuscounter]=1
#             gene_names.append(line["symbol"])
#             gene_entrezIDs_withFoundGeneNames.append(UniprotsForMissingGeneSymbols[cntr])
#             statuscounter=statuscounter+1
#             cntr=cntr+1
#         except KeyError:
#             EntrezIDs_Status[statuscounter]=0
#             NotFoundGenes.append(UniprotsForMissingGeneSymbols[statuscounter])
#             statuscounter=statuscounter+1
#             pass
# except:
#     pass


# NumberOfFoundGenes=len(gene_names)
# # NumberOfNotFoundGenes=noOfGenes-NumberOfFoundGenes
# NumberOfNotFoundGenes=len(NotFoundGenes)

# for i in range(NumberOfFoundGenes):
#     bait_uniprot_list1.append(gene_entrezIDs_withFoundGeneNames)
#     bait_symbol_list1.append(gene_names)
#     bait_usage_list1.append(DICT[gene_entrezIDs_withFoundGeneNames[i]])
    

# data = np.column_stack([bait_uniprot_list1, bait_symbol_list1, bait_usage_list1])
# datafile_path = "BAIT_USAGE-FINAL-GENESYMBOLS.txt"
# np.savetxt(datafile_path , data, fmt=['%s', '%s', '%d'])


# ########################## Faulty (with NANs): ##########################


# bait_usage = pd.read_csv('/Users/surya/Desktop/UpdateFunction-StudyBiasProject-localCopy/bait_usage_intact.txt', sep='\t')

# all_genes = list(set(bait_usage.bait_symbol))
# gene_occurences = {gene: 0 for gene in all_genes}
# gene_bait_usage = {gene: 0 for gene in all_genes}
# for i in range(bait_usage.shape[0]):
#     gene = bait_usage.loc[i, 'bait_symbol']
#     gene_bait_usage[gene] += bait_usage.loc[i, 'bait_usage']
#     gene_occurences[gene] += 1
# gene_bait_usages = [gene_bait_usage[gene] for gene in all_genes]
# gene_bait_usage = pd.DataFrame({'gene': all_genes, 'bait_usage': gene_bait_usages})
# gene_bait_usage.set_index('gene',inplace=True)

# duplicates = [gene for gene in all_genes if gene_occurences[gene] > 1]

# ratio_of_duplicate_gene_symbols=len(duplicates) / len(all_genes)

# gene_bait_usage.to_csv('gene_bait_usage.csv')


################################### Correct (without NANs): ###################################


# bait_usage_data = pd.read_csv('/Users/surya/Desktop/UpdateFunction-StudyBiasProject-localCopy/bait_usage_intact.txt', sep='\t')

# bait_uniprot=bait_usage_data['bait_uniprot']
# bait_symbol=bait_usage_data['bait_symbol']
# bait_usage=bait_usage_data['bait_usage']

# bait_uniprot_list=bait_uniprot.to_list()
# bait_symbol_list=bait_symbol.to_list()
# bait_usage_list=bait_usage.to_list()


# # Bait usage data (with GENE_SYMBOL):
# # -----------------------------------

# LENGTH_orig=len(bait_uniprot_list)
# UniprotsForMissingGeneSymbols=[]

# bait_uniprot_list1=[]
# bait_symbol_list1=[]
# bait_usage_list1=[]
# for i in range(LENGTH_orig):
#     if not(isinstance(bait_symbol_list[i], str)):
#         UniprotsForMissingGeneSymbols.append(bait_uniprot_list[i])
#     else:
#         bait_uniprot_list1.append(bait_uniprot_list[i])
#         bait_symbol_list1.append(bait_symbol_list[i])
#         bait_usage_list1.append(bait_usage_list[i])


# bait_usage_data_NoMissingSymbols = pd.DataFrame([bait_uniprot_list1,bait_symbol_list1,bait_usage_list1])
# bait_usage_data_NoMissingSymbols = bait_usage_data_NoMissingSymbols.transpose()
# bait_usage_data_NoMissingSymbols.columns=['bait_uniprot','bait_symbol','bait_usage']


# all_genes = list(set(bait_usage_data_NoMissingSymbols.bait_symbol))
# gene_occurences = {gene: 0 for gene in all_genes}
# gene_bait_usage = {gene: 0 for gene in all_genes}
# for i in range(bait_usage_data_NoMissingSymbols.shape[0]):
#     gene = bait_usage_data_NoMissingSymbols.loc[i, 'bait_symbol']
#     gene_bait_usage[gene] += bait_usage_data_NoMissingSymbols.loc[i, 'bait_usage']
#     gene_occurences[gene] += 1
# gene_bait_usages = [gene_bait_usage[gene] for gene in all_genes]
# gene_bait_usage = pd.DataFrame({'gene': all_genes, 'bait_usage': gene_bait_usages})
# gene_bait_usage.set_index('gene',inplace=True)

# duplicates = [gene for gene in all_genes if gene_occurences[gene] > 1]
# ratio_of_duplicate_gene_symbols=len(duplicates) / len(all_genes)

# gene_bait_usage.to_csv('geneSymbol_baitUsage.csv')

#####################################################################################################

PAIR_FREQ_DATA=pd.read_csv('/Users/surya/Desktop/StudyBiasPaper/CODE/Robust/robust-main/bias-aware data/pair_study_frequency.txt', sep=' ')
IDs_interactor_A=PAIR_FREQ_DATA['IDs_interactor_A']
IDs_interactor_B=PAIR_FREQ_DATA['IDs_interactor_B']
freq=PAIR_FREQ_DATA['freq']
symbol_A=PAIR_FREQ_DATA['symbol_A']
symbol_B=PAIR_FREQ_DATA['symbol_B']

IDs_interactor_A=IDs_interactor_A.to_list()
IDs_interactor_B=IDs_interactor_B.to_list()
freq=freq.to_list()
symbol_A=symbol_A.to_list()
symbol_B=symbol_B.to_list()

LENGTH_orig=len(freq)

list_status_A=np.ones((LENGTH_orig,1))
list_status_B=np.ones((LENGTH_orig,1))

MissingIDs_A=[]
MissingIDs_B=[]

IDs_A=[]
IDs_B=[]
symbols_A=[]
symbols_B=[]
freqs=[]

for i in range(LENGTH_orig):
    if isinstance(symbol_A[i], str) and isinstance(symbol_B[i], str):
        IDs_A.append(IDs_interactor_A[i])
        IDs_B.append(IDs_interactor_B[i])
        symbols_A.append(symbol_A[i])
        symbols_B.append(symbol_B[i])
        freqs.append(freq[i])
    else:
        if not(isinstance(symbol_A[i], str)):
           list_status_A[i]=0
           MissingIDs_A.append(IDs_interactor_A[i])
        if not(isinstance(symbol_B[i], str)):
           list_status_B[i]=0
           MissingIDs_B.append(IDs_interactor_B[i]) 
        
noOfMissingGeneSymbols_A=len(MissingIDs_A)
noOfMissingGeneSymbols_B=len(MissingIDs_B)


pair_freqs_NoMissingSymbols = pd.DataFrame([IDs_A,IDs_B,freqs,symbols_A,symbols_B])
pair_freqs_NoMissingSymbols = pair_freqs_NoMissingSymbols.transpose()
pair_freqs_NoMissingSymbols.columns=['IDs_interactor_A','IDs_interactor_B','freq','symbol_A','symbol_B']


all_genes = list(set(pair_freqs_NoMissingSymbols.symbol_A).union(set(pair_freqs_NoMissingSymbols.symbol_B)))
study_attention = {gene: 0 for gene in all_genes}

for i in range(pair_freqs_NoMissingSymbols.shape[0]):
    gene = pair_freqs_NoMissingSymbols.loc[i,'symbol_A']
    if type(gene) == str:
        study_attention[gene] += 1
    gene = pair_freqs_NoMissingSymbols.loc[i,'symbol_B']
    if type(gene) == str:
        study_attention[gene] += 1

genes = [gene for gene, _ in study_attention.items()]
counts = [count for _, count in study_attention.items()]
study_att = pd.DataFrame({'gene': genes, 'study_attention': counts})
study_att.set_index('gene', inplace=True)
study_att.to_csv('study_attention.csv')

#########################################################################################################

PAIR_FREQ_DATA=pd.read_csv('/Users/surya/Desktop/StudyBiasPaper/CODE/Robust/robust-main/bias-aware data/pair_study_frequency.txt', sep=' ')
IDs_interactor_A=PAIR_FREQ_DATA['IDs_interactor_A']
IDs_interactor_B=PAIR_FREQ_DATA['IDs_interactor_B']
freq=PAIR_FREQ_DATA['freq']
symbol_A=PAIR_FREQ_DATA['symbol_A']
symbol_B=PAIR_FREQ_DATA['symbol_B']

IDs_interactor_A=IDs_interactor_A.to_list()
IDs_interactor_B=IDs_interactor_B.to_list()
freq=freq.to_list()
symbol_A=symbol_A.to_list()
symbol_B=symbol_B.to_list()

LENGTH_orig=len(freq)

list_status_A=np.ones((LENGTH_orig,1))
list_status_B=np.ones((LENGTH_orig,1))

MissingIDs_A=[]
MissingIDs_B=[]

IDs_A=[]
IDs_B=[]
symbols_A=[]
symbols_B=[]
freqs=[]

for i in range(LENGTH_orig):
    if isinstance(symbol_A[i], str) and isinstance(symbol_B[i], str):
        IDs_A.append(IDs_interactor_A[i])
        IDs_B.append(IDs_interactor_B[i])
        symbols_A.append(symbol_A[i])
        symbols_B.append(symbol_B[i])
        freqs.append(freq[i])
    else:
        if not(isinstance(symbol_A[i], str)):
           list_status_A[i]=0
           MissingIDs_A.append(IDs_interactor_A[i])
        if not(isinstance(symbol_B[i], str)):
           list_status_B[i]=0
           MissingIDs_B.append(IDs_interactor_B[i]) 
        
noOfMissingGeneSymbols_A=len(MissingIDs_A)
noOfMissingGeneSymbols_B=len(MissingIDs_B)


pair_freqs_NoMissingSymbols = pd.DataFrame([IDs_A,IDs_B,freqs,symbols_A,symbols_B])
pair_freqs_NoMissingSymbols = pair_freqs_NoMissingSymbols.transpose()
pair_freqs_NoMissingSymbols.columns=['IDs_interactor_A','IDs_interactor_B','freq','symbol_A','symbol_B']


all_genes = list(set(pair_freqs_NoMissingSymbols.symbol_A).union(set(pair_freqs_NoMissingSymbols.symbol_B)))
study_attention = {gene: 0 for gene in all_genes}

for i in range(pair_freqs_NoMissingSymbols.shape[0]):
    gene = pair_freqs_NoMissingSymbols.loc[i,'symbol_A']
    if type(gene) == str:
        study_attention[gene] += 1
    gene = pair_freqs_NoMissingSymbols.loc[i,'symbol_B']
    if type(gene) == str:
        study_attention[gene] += 1

genes = [gene for gene, _ in study_attention.items()]
counts = [count for _, count in study_attention.items()]
study_att = pd.DataFrame({'gene': genes, 'study_attention': counts})
study_att.set_index('gene', inplace=True)
study_att.to_csv('study_attention.csv')

#########################################################################################################


# bait_usage_data = pd.read_csv('bait_usage_intact.txt', sep='\t')

# bait_uniprot=bait_usage_data['bait_uniprot']
# bait_symbol=bait_usage_data['bait_symbol']
# bait_usage=bait_usage_data['bait_usage']

# bait_uniprot_list=bait_uniprot.to_list()
# bait_symbol_list=bait_symbol.to_list()
# bait_usage_list=bait_usage.to_list()


# # Bait usage data (with GENE_SYMBOL):
# # -----------------------------------

# LENGTH_orig=len(bait_uniprot_list)
# UniprotsForMissingGeneSymbols=[]

# bait_uniprot_list1=[]
# bait_symbol_list1=[]
# bait_usage_list1=[]
# for i in range(LENGTH_orig):
#     if not(isinstance(bait_symbol_list[i], str)):
#         UniprotsForMissingGeneSymbols.append(bait_uniprot_list[i])
#     else:
#         bait_uniprot_list1.append(bait_uniprot_list[i])
#         bait_symbol_list1.append(bait_symbol_list[i])
#         bait_usage_list1.append(bait_usage_list[i])


# bait_usage_data_NoMissingSymbols = pd.DataFrame([bait_uniprot_list1,bait_symbol_list1,bait_usage_list1])
# bait_usage_data_NoMissingSymbols = bait_usage_data_NoMissingSymbols.transpose()
# bait_usage_data_NoMissingSymbols.columns=['bait_uniprot','bait_symbol','bait_usage']


# all_genes = list(set(bait_usage_data_NoMissingSymbols.bait_symbol))
# gene_occurences = {gene: 0 for gene in all_genes}
# gene_bait_usage = {gene: 0 for gene in all_genes}
# for i in range(bait_usage_data_NoMissingSymbols.shape[0]):
#     gene = bait_usage_data_NoMissingSymbols.loc[i, 'bait_symbol']
#     gene_bait_usage[gene] += bait_usage_data_NoMissingSymbols.loc[i, 'bait_usage']
#     gene_occurences[gene] += 1
# gene_bait_usages = [gene_bait_usage[gene] for gene in all_genes]
# gene_bait_usage = pd.DataFrame({'gene': all_genes, 'bait_usage': gene_bait_usages})
# gene_bait_usage.set_index('gene',inplace=True)

# # duplicates = [gene for gene in all_genes if gene_occurences[gene] > 1]
# # ratio_of_duplicate_gene_symbols=len(duplicates) / len(all_genes)

# gene_bait_usage.to_csv('geneSymbol_baitUsage.csv')



