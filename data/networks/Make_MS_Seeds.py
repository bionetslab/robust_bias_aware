import pandas as pd
import networkx as nx
import mygene
import random

ms_seeds=pd.read_csv('../ms_seeds.txt', header=None)
ms_seeds_list=ms_seeds[0].tolist()

uniprots=list(set(ms_seeds_list)) # seeds in uniprot
gene_symbols=[]
entrez=[]

# ---------------------------------- Convert seeds to gene_symbol ----------------------------------

try:
    mg = mygene.MyGeneInfo()
    # out = mg.querymany(nodes_, scopes= 'symbol', fields='entrezgene', species='human', verbose=False)
    out = mg.querymany(uniprots, scopes= 'uniprot', fields='symbol', species='human', verbose=False)
except:
    pass

NODES_uniprot=[]
NODES_genesymbol=[]
for i in range(len(out)):
    try:
        NODES_genesymbol.append(out[i]['symbol'])
        NODES_uniprot.append(out[i]['query'])
    except:
        pass

uniprot_genesymbol_DICT=dict(zip(NODES_uniprot, NODES_genesymbol))

# ---------------------------------- Convert seeds to entrez ----------------------------------

try:
    mg = mygene.MyGeneInfo()
    # out = mg.querymany(nodes_, scopes= 'symbol', fields='entrezgene', species='human', verbose=False)
    out = mg.querymany(uniprots, scopes= 'uniprot', fields='entrezgene', species='human', verbose=False)
except:
    pass

NODES_uniprot=[]
NODES_entrez=[]
for i in range(len(out)):
    try:
        NODES_entrez.append(out[i]['entrezgene'])
        NODES_uniprot.append(out[i]['query'])
    except:
        pass

uniprot_entrez_DICT=dict(zip(NODES_uniprot, NODES_entrez))

# ------------------------------------- APID (MS seeds) -------------------------------------

UNIPROT_APID=nx.read_graphml('networks/UNIPROT/UNIPROT_APID.txt')
UNIPROT_APID_nodes=list(set(UNIPROT_APID.nodes))
UNIPROT_APID_ms_seeds=[]
for i in uniprots:
    if i in UNIPROT_APID_nodes:
        UNIPROT_APID_ms_seeds.append(i)
with open(r'UNIPROT_APID_ms_seeds.txt', 'w') as fp:
    for item in UNIPROT_APID_ms_seeds:
        fp.write("%s\n" % item)
try:
    UNIPROT_APID_ms_seeds_20=random.choices(UNIPROT_APID_ms_seeds, k=20)
    with open(r'UNIPROT_APID_ms_seeds_20.txt', 'w') as fp:
        for item in UNIPROT_APID_ms_seeds_20:
            fp.write("%s\n" % item)
    
    UNIPROT_APID_ms_seeds_10=random.choices(UNIPROT_APID_ms_seeds, k=10)
    with open(r'UNIPROT_APID_ms_seeds_10.txt', 'w') as fp:
        for item in UNIPROT_APID_ms_seeds_10:
            fp.write("%s\n" % item)
    
    UNIPROT_APID_ms_seeds_5=random.choices(UNIPROT_APID_ms_seeds, k=5)
    with open(r'UNIPROT_APID_ms_seeds_5.txt', 'w') as fp:
        for item in UNIPROT_APID_ms_seeds_5:
            fp.write("%s\n" % item)
except:
    pass


GENE_SYMBOL_APID=nx.read_graphml('networks/GENE_SYMBOL/GENE_SYMBOL_APID.txt')
GENE_SYMBOL_APID_nodes=list(set(GENE_SYMBOL_APID.nodes))
GENE_SYMBOL_APID_ms_seeds=[]
for i, j in uniprot_genesymbol_DICT:
    if j in GENE_SYMBOL_APID_nodes:
        GENE_SYMBOL_APID_ms_seeds.append(j)
with open(r'GENE_SYMBOL_APID_ms_seeds.txt', 'w') as fp:
    for item in GENE_SYMBOL_APID_ms_seeds:
        fp.write("%s\n" % item)

try:
    GENE_SYMBOL_APID_ms_seeds_20=random.choices(GENE_SYMBOL_APID_ms_seeds, k=20)
    with open(r'GENE_SYMBOL_APID_ms_seeds_20.txt', 'w') as fp:
        for item in GENE_SYMBOL_APID_ms_seeds_20:
            fp.write("%s\n" % item)
    
    GENE_SYMBOL_APID_ms_seeds_10=random.choices(GENE_SYMBOL_APID_ms_seeds, k=10)
    with open(r'GENE_SYMBOL_APID_ms_seeds_10.txt', 'w') as fp:
        for item in GENE_SYMBOL_APID_ms_seeds_10:
            fp.write("%s\n" % item)
    
    GENE_SYMBOL_APID_ms_seeds_5=random.choices(GENE_SYMBOL_APID_ms_seeds, k=5)
    with open(r'GENE_SYMBOL_APID_ms_seeds_5.txt', 'w') as fp:
        for item in GENE_SYMBOL_APID_ms_seeds_5:
            fp.write("%s\n" % item)
except:
    pass


ENTREZ_APID=nx.read_graphml('networks/ENTREZ/ENTREZ_APID.txt')
ENTREZ_APID_nodes=list(set(ENTREZ_APID.nodes))
ENTREZ_APID_ms_seeds=[]
for i, j in uniprot_entrez_DICT:
    if j in ENTREZ_APID_nodes:
        ENTREZ_APID_ms_seeds.append(j)
with open(r'ENTREZ_APID_ms_seeds.txt', 'w') as fp:
    for item in ENTREZ_APID_ms_seeds:
        fp.write("%s\n" % item)

try:
    ENTREZ_APID_ms_seeds_20=random.choices(ENTREZ_APID_ms_seeds, k=20)
    with open(r'ENTREZ_APID_ms_seeds_20.txt', 'w') as fp:
        for item in ENTREZ_APID_ms_seeds_20:
            fp.write("%s\n" % item)
    
    ENTREZ_APID_ms_seeds_10=random.choices(ENTREZ_APID_ms_seeds, k=10)
    with open(r'ENTREZ_APID_ms_seeds_10.txt', 'w') as fp:
        for item in ENTREZ_APID_ms_seeds_10:
            fp.write("%s\n" % item)
    
    ENTREZ_APID_ms_seeds_5=random.choices(ENTREZ_APID_ms_seeds, k=5)
    with open(r'ENTREZ_APID_ms_seeds_5.txt', 'w') as fp:
        for item in ENTREZ_APID_ms_seeds_5:
            fp.write("%s\n" % item)
    
except:
    pass
        
# ------------------------------------- BioGRID (MS seeds) -------------------------------------

UNIPROT_BioGRID=nx.read_graphml('networks/UNIPROT/UNIPROT_BioGRID.txt')
UNIPROT_BioGRID_nodes=list(set(UNIPROT_BioGRID.nodes))
UNIPROT_BioGRID_ms_seeds=[]
for i in uniprots:
    if i in UNIPROT_BioGRID_nodes:
        UNIPROT_BioGRID_ms_seeds.append(i)
with open(r'UNIPROT_BioGRID_ms_seeds.txt', 'w') as fp:
    for item in UNIPROT_BioGRID_ms_seeds:
        fp.write("%s\n" % item)
try:
    UNIPROT_BioGRID_ms_seeds_20=random.choices(UNIPROT_BioGRID_ms_seeds, k=20)
    with open(r'UNIPROT_BioGRID_ms_seeds_20.txt', 'w') as fp:
        for item in UNIPROT_BioGRID_ms_seeds_20:
            fp.write("%s\n" % item)
    
    
    UNIPROT_BioGRID_ms_seeds_10=random.choices(UNIPROT_BioGRID_ms_seeds, k=10)
    with open(r'UNIPROT_BioGRID_ms_seeds_10.txt', 'w') as fp:
        for item in UNIPROT_BioGRID_ms_seeds_10:
            fp.write("%s\n" % item)
            
    
    UNIPROT_BioGRID_ms_seeds_5=random.choices(UNIPROT_BioGRID_ms_seeds, k=5)
except:
    pass


GENE_SYMBOL_BioGRID=nx.read_graphml('networks/GENE_SYMBOL/GENE_SYMBOL_BioGRID.txt')
GENE_SYMBOL_BioGRID_nodes=list(set(GENE_SYMBOL_BioGRID.nodes))
GENE_SYMBOL_BioGRID_ms_seeds=[]
for i, j in uniprot_genesymbol_DICT:
    if j in GENE_SYMBOL_BioGRID_nodes:
        GENE_SYMBOL_BioGRID_ms_seeds.append(j)
try:
    GENE_SYMBOL_BioGRID_ms_seeds_20=random.choices(GENE_SYMBOL_BioGRID_ms_seeds, k=20)
    GENE_SYMBOL_BioGRID_ms_seeds_10=random.choices(GENE_SYMBOL_BioGRID_ms_seeds, k=10)
    GENE_SYMBOL_BioGRID_ms_seeds_5=random.choices(GENE_SYMBOL_BioGRID_ms_seeds, k=5)
except:
    pass


ENTREZ_BioGRID=nx.read_graphml('networks/ENTREZ/ENTREZ_BioGRID.txt')
ENTREZ_BioGRID_nodes=list(set(ENTREZ_BioGRID.nodes))
ENTREZ_BioGRID_ms_seeds=[]
for i, j in uniprot_entrez_DICT:
    if j in ENTREZ_BioGRID_nodes:
        ENTREZ_BioGRID_ms_seeds.append(j)
with open(r'ENTREZ_BioGRID_ms_seeds.txt', 'w') as fp:
    for item in ENTREZ_BioGRID_ms_seeds:
        fp.write("%s\n" % item)
try:
    ENTREZ_BioGRID_ms_seeds_20=random.choices(ENTREZ_BioGRID_ms_seeds, k=20)
    with open(r'ENTREZ_BioGRID_ms_seeds_20.txt', 'w') as fp:
        for item in ENTREZ_BioGRID_ms_seeds_20:
            fp.write("%s\n" % item)
    ENTREZ_BioGRID_ms_seeds_10=random.choices(ENTREZ_BioGRID_ms_seeds, k=10)
    with open(r'ENTREZ_BioGRID_ms_seeds_10.txt', 'w') as fp:
        for item in ENTREZ_BioGRID_ms_seeds_10:
            fp.write("%s\n" % item)
    ENTREZ_BioGRID_ms_seeds_5=random.choices(ENTREZ_BioGRID_ms_seeds, k=5)
    with open(r'ENTREZ_BioGRID_ms_seeds_5.txt', 'w') as fp:
        for item in ENTREZ_BioGRID_ms_seeds_5:
            fp.write("%s\n" % item)
except:
    pass

# ------------------------------------- HPRD (MS seeds) -------------------------------------

UNIPROT_HPRD=nx.read_graphml('networks/UNIPROT/UNIPROT_HPRD.txt')
UNIPROT_HPRD_nodes=list(set(UNIPROT_HPRD.nodes))
UNIPROT_HPRD_ms_seeds=[]
for i in uniprots:
    if i in UNIPROT_HPRD_nodes:
        UNIPROT_HPRD_ms_seeds.append(i)
try:
    UNIPROT_HPRD_ms_seeds_20=random.choices(UNIPROT_HPRD_ms_seeds, k=20)
    UNIPROT_HPRD_ms_seeds_10=random.choices(UNIPROT_HPRD_ms_seeds, k=10)
    UNIPROT_HPRD_ms_seeds_5=random.choices(UNIPROT_HPRD_ms_seeds, k=5)
except:
    pass


GENE_SYMBOL_HPRD=nx.read_graphml('networks/GENE_SYMBOL/GENE_SYMBOL_HPRD.txt')
GENE_SYMBOL_HPRD_nodes=list(set(GENE_SYMBOL_HPRD.nodes))
GENE_SYMBOL_HPRD_ms_seeds=[]
for i, j in uniprot_genesymbol_DICT:
    if j in GENE_SYMBOL_HPRD_nodes:
        GENE_SYMBOL_HPRD_ms_seeds.append(j)
try:
    GENE_SYMBOL_HPRD_ms_seeds_20=random.choices(GENE_SYMBOL_HPRD_ms_seeds, k=20)
    GENE_SYMBOL_HPRD_ms_seeds_10=random.choices(GENE_SYMBOL_HPRD_ms_seeds, k=10)
    GENE_SYMBOL_HPRD_ms_seeds_5=random.choices(GENE_SYMBOL_HPRD_ms_seeds, k=5)
except:
    pass


ENTREZ_HPRD=nx.read_graphml('networks/ENTREZ/ENTREZ_HPRD.txt')
ENTREZ_HPRD_nodes=list(set(ENTREZ_HPRD.nodes))
ENTREZ_HPRD_ms_seeds=[]
for i, j in uniprot_entrez_DICT:
    if j in ENTREZ_HPRD_nodes:
        ENTREZ_HPRD_ms_seeds.append(j)
try:
    ENTREZ_HPRD_ms_seeds_20=random.choices(ENTREZ_HPRD_ms_seeds, k=20)
    ENTREZ_HPRD_ms_seeds_10=random.choices(ENTREZ_HPRD_ms_seeds, k=10)
    ENTREZ_HPRD_ms_seeds_5=random.choices(ENTREZ_HPRD_ms_seeds, k=5)
except:
    pass

# ------------------------------------- STRING (MS seeds) -------------------------------------

UNIPROT_STRING=nx.read_graphml('networks/UNIPROT/UNIPROT_STRING.txt')
UNIPROT_STRING_nodes=list(set(UNIPROT_STRING.nodes))
UNIPROT_STRING_ms_seeds=[]
for i in uniprots:
    if i in UNIPROT_STRING_nodes:
        UNIPROT_STRING_ms_seeds.append(i)
try:
    UNIPROT_STRING_ms_seeds_20=random.choices(UNIPROT_STRING_ms_seeds, k=20)
    UNIPROT_STRING_ms_seeds_10=random.choices(UNIPROT_STRING_ms_seeds, k=10)
    UNIPROT_STRING_ms_seeds_5=random.choices(UNIPROT_STRING_ms_seeds, k=5)
except:
    pass


GENE_SYMBOL_STRING=nx.read_graphml('networks/GENE_SYMBOL/GENE_SYMBOL_STRING.txt')
GENE_SYMBOL_STRING_nodes=list(set(GENE_SYMBOL_STRING.nodes))
GENE_SYMBOL_STRING_ms_seeds=[]
for i, j in uniprot_genesymbol_DICT:
    if j in GENE_SYMBOL_STRING_nodes:
        GENE_SYMBOL_STRING_ms_seeds.append(j)
try:
    GENE_SYMBOL_STRING_ms_seeds_20=random.choices(GENE_SYMBOL_STRING_ms_seeds, k=20)
    GENE_SYMBOL_STRING_ms_seeds_10=random.choices(GENE_SYMBOL_STRING_ms_seeds, k=10)
    GENE_SYMBOL_STRING_ms_seeds_5=random.choices(GENE_SYMBOL_STRING_ms_seeds, k=5)
except:
    pass

ENTREZ_STRING=nx.read_graphml('networks/ENTREZ/ENTREZ_STRING.txt')
ENTREZ_STRING_nodes=list(set(ENTREZ_STRING.nodes))
ENTREZ_STRING_ms_seeds=[]
for i, j in uniprot_entrez_DICT:
    if j in ENTREZ_STRING_nodes:
        ENTREZ_STRING_ms_seeds.append(j)
try:
    ENTREZ_STRING_ms_seeds_20=random.choices(ENTREZ_STRING_ms_seeds, k=20)
    ENTREZ_STRING_ms_seeds_10=random.choices(ENTREZ_STRING_ms_seeds, k=10)
    ENTREZ_STRING_ms_seeds_5=random.choices(ENTREZ_STRING_ms_seeds, k=5)
except:
    pass

# ---------------------------------------------------------------------------------------------------------------
