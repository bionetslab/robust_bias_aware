# -*- coding: utf-8 -*-
"""PlotResults.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1LJ-g1K9gDpiDHmoSuXwOZvCfyregEwas
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import mygene
import gseapy
# plt.rcParams.update(plt.rcParamsDefault)
# plt.rcParams['text.usetex'] = True

flatten = lambda l: [item for sublist in l for item in sublist]
# plt.rcParams['text.usetex'] = True

res_0point00=pd.read_csv('/Users/surya/Documents/RECORDRESULTS_STUDYBIAS/RESULTS/1. ADDITIVE/STUDY_ATTENTION/RESULTS-SUMMARIZED/all_results-0.00.csv')
res_0point00['lambda']=0.00

res_0point25=pd.read_csv('/Users/surya/Documents/RECORDRESULTS_STUDYBIAS/RESULTS/1. ADDITIVE/STUDY_ATTENTION/RESULTS-SUMMARIZED/all_results-0.25.csv')
res_0point25['lambda']=0.25

res_0point50=pd.read_csv('/Users/surya/Documents/RECORDRESULTS_STUDYBIAS/RESULTS/1. ADDITIVE/STUDY_ATTENTION/RESULTS-SUMMARIZED/all_results-0.50.csv')
res_0point50['lambda']=0.50

res_0point75=pd.read_csv('/Users/surya/Documents/RECORDRESULTS_STUDYBIAS/RESULTS/1. ADDITIVE/STUDY_ATTENTION/RESULTS-SUMMARIZED/all_results-0.75.csv')
res_0point75['lambda']=0.75

res_1point00=pd.read_csv('/Users/surya/Documents/RECORDRESULTS_STUDYBIAS/RESULTS/1. ADDITIVE/STUDY_ATTENTION/RESULTS-SUMMARIZED/all_results-1.00.csv')
res_1point00['lambda']=1.00

res = pd.concat([res_0point00, res_0point25, res_0point50, res_0point75, res_1point00])

index=pd.Series(range(len(res)))

# res = res.set_index(index)

res = res.set_index(index)

################################################################################################################################

SUM=[0]*13
COUNT=[0]*13

diseaseCondition=['ALS', 'LC', 'UC', 'HD', 'CD']
diseaseCodes=['hsa05014', 'hsa05223', 'hsa04060', 'hsa04630', 'hsa05321', 'hsa05016', 'hsa04621', 'hsa04140', 'C0002736', 'C1737250', 'C0009324', 'C0020179', 'C0021390']
diseaseIndex=[[0,8], [1,9], [2,3,4,10], [5,11], [6,2,3,4,7,12]]
diseaseCondition_diseaseIndex=dict(zip(diseaseCondition, diseaseIndex))

baitUsageData=pd.read_csv('../../../Get Disgenet and KEGG pathways/geneSymbol_baitUsage.txt', sep=' ')
geneSymbol_baitUsage=dict(zip(baitUsageData.gene, baitUsageData.bait_usage))

studyAttentionData=pd.read_csv('../../../Get Disgenet and KEGG pathways/geneSymbol_studyAttention.txt', sep=' ')
geneSymbol_studyAttention=dict(zip(studyAttentionData.gene, studyAttentionData.study_attention))

all_genes = list(res.result_genes)
genes = []
for i in range(len(all_genes)):
    try:
        genes.append(all_genes[i].split(","))
    except:
        pass

genes = list(set(flatten(genes)))
mg = mygene.MyGeneInfo()
out = mg.querymany(genes, scopes="entrezgene", fields='symbol', species='human', verbose=False)
mapping = dict()
for line in out:
    try:
        mapping[line["query"]] = line["symbol"]
    except KeyError:
        print("{0} was not mapped to any gene name".format(line["query"]))
        mapping[line["query"]] = line["query"]


AvgBaitUsage=[]
AvgStudyAttention=[]
# kegg_pval = []

for i in res.index :
    gs = str(res.result_genes[i])
    gs = gs.split(",")
    # gs = [mapping[x] for x in gs]
    GS=[]
    for k in gs:
        try:
            GS.append(mapping[k])
        except:
            pass
    for m in diseaseCondition_diseaseIndex[res.condition_name[i]]:
        SUM[m]+=res.mean_degree_result[i]*res.num_genes_result[i]
        COUNT[m]+=res.num_genes_result[i]
    
    sum_baitUsages=0
    len_baitUsages=0
    for j1 in GS:
        try:
            sum_baitUsages+=geneSymbol_baitUsage[j1]
            len_baitUsages+=1
        except:
            pass
    avgBaitUsage=sum_baitUsages/len_baitUsages
    AvgBaitUsage.append(avgBaitUsage)
    
    sum_studyAttentions=0
    len_studyAttentions=0
    for j2 in GS:
        try:
            sum_studyAttentions+=geneSymbol_studyAttention[j2]
            len_studyAttentions+=1
        except:
            pass
    avgStudyAttention=sum_baitUsages/len_baitUsages
    AvgStudyAttention.append(avgStudyAttention)
    
        
res['AverageBaitUsage']=AvgBaitUsage
res['AverageStudyAttention']=AvgStudyAttention


MeanNodeDegree = [i / j for i, j in zip(SUM, COUNT)]
    
diseaseCodes_MeanNodeDegree=dict(zip(diseaseCodes, MeanNodeDegree))

################################################################################################################################


kegg_disgenet_statistics = pd.read_csv('../../../Get Disgenet and KEGG pathways/kegg_disgenet_statistics.csv')

MeanNodeDegree_arranged=[]

for i in kegg_disgenet_statistics.index:
    MeanNodeDegree_arranged.append(MeanNodeDegree[diseaseCodes.index(kegg_disgenet_statistics.gene_set[i])])

kegg_disgenet_statistics['average_node_degree']=MeanNodeDegree_arranged

kegg_disgenet_statistics.to_csv('KEGGandDISGENETpathways_AvgNodeDegree.csv')

################################################################################################################################

original = res[res.network_generator_name == 'ORIGINAL']

# reference = pd.read_csv('KEGGandDISGENETpathways_AvgNodeDegree-BAITUSAGE(ROBUST2).csv')
# reference

reference=kegg_disgenet_statistics

# # -np.log10(0.05)

# ax.set(xlabel='x-axis label', ylabel='y-axis label')

fig, axes = plt.subplots(4, 2, figsize=(5, 9))
fig.suptitle('Results for additive edge costs based on study attention',fontsize=14)

sns.boxplot(ax=axes[0,0], data=original, x='lambda', y='disgenet_overlap')
ax1=axes[0,0]
ax1.set(xlabel='lambda', ylabel='DisGeNet overlap')

sns.boxplot(ax=axes[0,1], data=original, x='lambda', y='neg_log_gsea_p_value')
ax2=axes[0,1]
ax2.set(xlabel='lambda', ylabel='-log10 P(KEGG GSEA)')

sns.boxplot(ax=axes[1,0], data=original, x='lambda', y='mean_degree_result')
ax3=axes[1,0]
ax3.set(xlabel='lambda', ylabel='Mean node degree')

sns.boxplot(ax=axes[1,1], data=original, x='lambda', y='AverageBaitUsage')
ax4=axes[1,1]
ax4.set(xlabel='lambda', ylabel='Mean bait usage')

sns.boxplot(ax=axes[2,0], data=original, x='lambda', y='AverageStudyAttention', showfliers = False)
ax5=axes[2,0]
ax5.set(xlabel='lambda', ylabel='Mean study attention')

sns.boxplot(ax=axes[2,1], data=reference, x='data_source', y='average_node_degree')
ax6=axes[2,1]
ax6.set(xlabel='data source', ylabel='Mean node degree')

sns.boxplot(ax=axes[3,0], data=reference, x='data_source', y='mean_bait_usage')
ax7=axes[3,0]
ax7.set(xlabel='data source', ylabel='Mean bait usage')

sns.boxplot(ax=axes[3,1], data=reference, x='data_source', y='mean_study_attention')
ax8=axes[3,1]
ax8.set(xlabel='data source', ylabel='Mean study attention')

plt.tight_layout()
plt.show()

# for i in range(4):
#     for j in range(2):
#         if j == 0 or i == 0:
#             axes[i,j].set_xlabel(r'$\lambda$',fontsize=12)
# axes[0,0].set_ylabel('DisGeNET overlap',fontsize=12)
# axes[0,1].set_ylabel(r'$-\log_{10}P$ (KEGG GSEA)',fontsize=12)
# axes[1,0].set_ylabel('Mean node degree',fontsize=12)
# axes[1,1].set_ylabel('Mean node degree',fontsize=12)
# axes[1,1].set_xlabel('Data source',fontsize=12)
# axes[2,0].set_ylabel('Mean bait usage',fontsize=12)
# axes[2,1].set_ylabel('Mean bait usage',fontsize=12)
# axes[2,1].set_xlabel('Data source',fontsize=12)
# axes[3,0].set_ylabel('Mean study attention',fontsize=12)
# axes[3,1].set_ylabel('Mean study attention',fontsize=12)
# axes[3,1].set_xlabel('Data source',fontsize=12)
# fig.tight_layout()

# fig.savefig('results_overview.pdf')

# fig, axes = plt.subplots(1, 2, figsize=(5, 3))
# fig.suptitle('Results for additive edge costs based on bait usage',fontsize=14)
# sns.boxplot(ax=axes[0], data=original, x='lambda', y='mean_degree_result')
# sns.boxplot(ax=axes[1], data=reference, x='data_source', y='AverageNodeDegree')
            
            
# axes[0].set_xlabel(r'$\lambda$',fontsize=12)
# axes[0].set_ylabel('Mean node degree',fontsize=12)
# axes[1].set_ylabel('Mean node degree',fontsize=12)
# axes[1].set_xlabel('Data source',fontsize=12)
# fig.tight_layout()
# fig.savefig('node_degree.pdf')



# original.head()

# sns.boxplot(data=original, x='lambda', y='mean_degree_result')

# sns.relplot(data=original, x='lambda', y='neg_log_gsea_p_value', kind='line')

# sns.relplot(data=res, x='lambda', y='neg_log_gsea_p_value', col='network_generator_name', col_wrap=3, kind='line')

# sns.relplot(data=res, x='lambda', y='disgenet_overlap', col='network_generator_name', col_wrap=3, kind='line')