import json
import ndex2
import networkx as nx
import pandas as pd
import mygene


client = ndex2.client.Ndex2()

# (1.1) APID Human Interactome:
# APID = client.get_network_as_cx_stream('35fcb572-c566-11e8-aaa6-0ac135e8bacf')

# (1.2) APID Human Interactome (only human proteins):
APID = client.get_network_as_cx_stream('9c38ce6e-c564-11e8-aaa6-0ac135e8bacf')


# (2) BioGRID:
BioGRID = client.get_network_as_cx_stream('becec556-86d4-11e7-a10d-0ac135e8bacf')

# (3) HPRD:
HPRD = client.get_network_as_cx_stream('1093e665-86da-11e7-a10d-0ac135e8bacf')

# (4) STRING:
STRING = client.get_network_as_cx_stream('cfcd4cdb-86da-11e7-a10d-0ac135e8bacf')

# (5) HIPPIE:
HIPPIE = client.get_network_as_cx_stream('89dd3925-3718-11e9-9f06-0ac135e8bacf')

# (6) HuRI:
HuRI=client.get_network_as_cx_stream('73bc2c06-5fb2-11e9-9f06-0ac135e8bacf')


APID = ndex2.create_nice_cx_from_raw_cx(json.loads(APID.content))
BioGRID = ndex2.create_nice_cx_from_raw_cx(json.loads(BioGRID.content))
HPRD = ndex2.create_nice_cx_from_raw_cx(json.loads(HPRD.content))
STRING = ndex2.create_nice_cx_from_raw_cx(json.loads(STRING.content))
HIPPIE = ndex2.create_nice_cx_from_raw_cx(json.loads(HIPPIE.content))
HuRI = ndex2.create_nice_cx_from_raw_cx(json.loads(HuRI.content))


APID = APID.to_networkx(mode='default')
BioGRID = BioGRID.to_networkx(mode='default')
HPRD = HPRD.to_networkx(mode='default')
STRING = STRING.to_networkx(mode='default')
HIPPIE = HIPPIE.to_networkx(mode='default')
HuRI = HuRI.to_networkx(mode='default')

####################################################################### APID: #######################################################################

src=[]
dest=[]
for u, v in nx.get_edge_attributes(APID,'name').items():
    nodes_=v.split(' (interacts with) ')
    src.append(nodes_[0])
    dest.append(nodes_[1])
APID_UNIPROT = {'node1':src,'node2':dest}
APID_UNIPROT = pd.DataFrame(APID_UNIPROT)
APID_UNIPROT.to_csv('UNIPROT/APID.txt', sep=' ', index=False)

# ------------------------------------------------------------------------

src_set=set(src)
dest_set=set(dest)
nodes_set=src_set.union(dest_set)
nodes_=list(nodes_set)

try:
    mg = mygene.MyGeneInfo()
    # out = mg.querymany(nodes_, scopes= 'symbol', fields='entrezgene', species='human', verbose=False)
    out = mg.querymany(nodes_, scopes= 'uniprot', fields='symbol', species='human', verbose=False)
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

SRC=[]
DEST=[]

for i in src:
    try:
        SRC.append(uniprot_genesymbol_DICT[i])
    except:
        SRC.append(i)

for i in dest:
    try:
        DEST.append(uniprot_genesymbol_DICT[i])
    except:
        DEST.append(i)

APID_GENE_SYMBOL={'node1':SRC,'node2':DEST}
APID_GENE_SYMBOL = pd.DataFrame(APID_GENE_SYMBOL)
APID_GENE_SYMBOL.to_csv('GENE_SYMBOL/APID.txt', sep=' ', index=False)

# ----------------------------------------------------------------------

src_set=set(src)
dest_set=set(dest)
nodes_set=src_set.union(dest_set)
nodes_=list(nodes_set)

try:
    mg = mygene.MyGeneInfo()
    out = mg.querymany(nodes_, scopes= 'uniprot', fields='entrezgene', species='human', verbose=False)
except:
    pass


NODES_uniprot=[]
NODES_entrezgene=[]
for i in range(len(out)):
    try:
        NODES_entrezgene.append(out[i]['entrezgene'])
        NODES_uniprot.append(out[i]['query'])
    except:
        pass

uniprot_entrezgene_DICT=dict(zip(NODES_uniprot, NODES_entrezgene))

SRC=[]
DEST=[]

for i in src:
    try:
        SRC.append(uniprot_entrezgene_DICT[i])
    except:
        SRC.append(i)

for i in dest:
    try:
        DEST.append(uniprot_entrezgene_DICT[i])
    except:
        DEST.append(i)

APID_ENTREZGENE={'node1':SRC,'node2':DEST}
APID_ENTREZGENE = pd.DataFrame(APID_ENTREZGENE)
APID_ENTREZGENE.to_csv('ENTREZ/APID.txt', sep=' ', index=False)

# ####################################################################### BIOGRID: #######################################################################

nodes_=[]
node_attributes=[]
for u, v in nx.get_node_attributes(BioGRID,'represents').items():
    symbols_=v.split('hgnc.symbol:')
    symbol=symbols_[1]
    nodes_.append(u)
    node_attributes.append(symbol)
nodes_attributes_DICT=dict(zip(nodes_, node_attributes))

edges_=BioGRID.edges
edges_list=list(edges_)
LIST=list(map(list, zip(*edges_list)))

src=[]
dest=[]
for i in LIST[0]:
    src.append(nodes_attributes_DICT[i])
for i in LIST[1]:
    dest.append(nodes_attributes_DICT[i])

BioGRID_GENE_SYMBOL={'node1':src,'node2':dest}
BioGRID_GENE_SYMBOL = pd.DataFrame(BioGRID_GENE_SYMBOL)
BioGRID_GENE_SYMBOL.to_csv('GENE_SYMBOL/BioGRID.txt', sep=' ', index=False)

# ----------------------------------------------------------------------

src_set=set(src)
dest_set=set(dest)
nodes_set=src_set.union(dest_set)
nodes_=list(nodes_set)

try:
    mg = mygene.MyGeneInfo()
    # out = mg.querymany(nodes_, scopes= 'symbol', fields='entrezgene', species='human', verbose=False)
    out = mg.querymany(nodes_, scopes= 'symbol', fields='uniprot', species='human', verbose=False)
except:
    pass

nodes_genesymbol=[]
nodes_uniprot=[]
for i in range(len(out)):
    try:
        nodes_uniprot.append(out[i]['uniprot'])
        nodes_genesymbol.append(out[i]['query'])
    except:
        pass

NODES_genesymbol=[]

NODES_uniprot=[]

for i in range(len(nodes_uniprot)):
    try:
        NODES_uniprot.append(nodes_uniprot[i]['Swiss-Prot'])
        NODES_genesymbol.append(nodes_genesymbol[i])
    except:
        pass

genesymbol_uniprot_DICT=dict(zip(NODES_genesymbol, NODES_uniprot))

SRC=[]
DEST=[]

for i in src:
    try:
        SRC.append(genesymbol_uniprot_DICT[i])
    except:
        SRC.append(i)

for i in dest:
    try:
        DEST.append(genesymbol_uniprot_DICT[i])
    except:
        DEST.append(i)

BioGRID_UNIPROT={'node1':SRC,'node2':DEST}
BioGRID_UNIPROT = pd.DataFrame(BioGRID_UNIPROT)
BioGRID_UNIPROT.to_csv('UNIPROT/BioGRID.txt', sep=' ', index=False)

# ----------------------------------------------------------------------

src_set=set(src)
dest_set=set(dest)
nodes_set=src_set.union(dest_set)
nodes_=list(nodes_set)

try:
    mg = mygene.MyGeneInfo()
    out = mg.querymany(nodes_, scopes= 'symbol', fields='entrezgene', species='human', verbose=False)
except:
    pass


NODES_genesymbol=[]
NODES_entrezgene=[]
for i in range(len(out)):
    try:
        NODES_entrezgene.append(out[i]['entrezgene'])
        NODES_genesymbol.append(out[i]['query'])
    except:
        pass

genesymbol_entrezgene_DICT=dict(zip(NODES_genesymbol, NODES_entrezgene))

SRC=[]
DEST=[]

for i in src:
    try:
        SRC.append(genesymbol_entrezgene_DICT[i])
    except:
        SRC.append(i)

for i in dest:
    try:
        DEST.append(genesymbol_entrezgene_DICT[i])
    except:
        DEST.append(i)

HPRD_ENTREZGENE={'node1':SRC,'node2':DEST}
HPRD_ENTREZGENE = pd.DataFrame(HPRD_ENTREZGENE)
HPRD_ENTREZGENE.to_csv('ENTREZ/BioGRID.txt', sep=' ', index=False)

# ####################################################################### HPRD: #######################################################################

nodes_=[]
node_attributes=[]
for u, v in nx.get_node_attributes(HPRD,'represents').items():
    symbols_=v.split('hgnc.symbol:')
    symbol=symbols_[1]
    nodes_.append(u)
    node_attributes.append(symbol)
nodes_attributes_DICT=dict(zip(nodes_, node_attributes))

edges_=HPRD.edges
edges_list=list(edges_)
LIST=list(map(list, zip(*edges_list)))

src=[]
dest=[]
for i in LIST[0]:
    src.append(nodes_attributes_DICT[i])
for i in LIST[1]:
    dest.append(nodes_attributes_DICT[i])

HPRD_GENE_SYMBOL={'node1':src,'node2':dest}
HPRD_GENE_SYMBOL = pd.DataFrame(HPRD_GENE_SYMBOL)
HPRD_GENE_SYMBOL.to_csv('GENE_SYMBOL/HPRD.txt', sep=' ', index=False)

# ----------------------------------------------------------------------

src_set=set(src)
dest_set=set(dest)
nodes_set=src_set.union(dest_set)
nodes_=list(nodes_set)

try:
    mg = mygene.MyGeneInfo()
    # out = mg.querymany(nodes_, scopes= 'symbol', fields='entrezgene', species='human', verbose=False)
    out = mg.querymany(nodes_, scopes= 'symbol', fields='uniprot', species='human', verbose=False)
except:
    pass

nodes_genesymbol=[]
nodes_uniprot=[]
for i in range(len(out)):
    try:
        nodes_uniprot.append(out[i]['uniprot'])
        nodes_genesymbol.append(out[i]['query'])
    except:
        pass

NODES_genesymbol=[]

NODES_uniprot=[]

for i in range(len(nodes_uniprot)):
    try:
        NODES_uniprot.append(nodes_uniprot[i]['Swiss-Prot'])
        NODES_genesymbol.append(nodes_genesymbol[i])
    except:
        pass

genesymbol_uniprot_DICT=dict(zip(NODES_genesymbol, NODES_uniprot))

SRC=[]
DEST=[]

for i in src:
    try:
        SRC.append(genesymbol_uniprot_DICT[i])
    except:
        SRC.append(i)

for i in dest:
    try:
        DEST.append(genesymbol_uniprot_DICT[i])
    except:
        DEST.append(i)

HPRD_UNIPROT={'node1':SRC,'node2':DEST}
HPRD_UNIPROT = pd.DataFrame(HPRD_UNIPROT)
HPRD_UNIPROT.to_csv('UNIPROT/HPRD.txt', sep=' ', index=False)

# ----------------------------------------------------------------------

src_set=set(src)
dest_set=set(dest)
nodes_set=src_set.union(dest_set)
nodes_=list(nodes_set)

try:
    mg = mygene.MyGeneInfo()
    out = mg.querymany(nodes_, scopes= 'symbol', fields='entrezgene', species='human', verbose=False)
except:
    pass


NODES_genesymbol=[]
NODES_entrezgene=[]
for i in range(len(out)):
    try:
        NODES_entrezgene.append(out[i]['entrezgene'])
        NODES_genesymbol.append(out[i]['query'])
    except:
        pass

genesymbol_entrezgene_DICT=dict(zip(NODES_genesymbol, NODES_entrezgene))

SRC=[]
DEST=[]

for i in src:
    try:
        SRC.append(genesymbol_entrezgene_DICT[i])
    except:
        SRC.append(i)

for i in dest:
    try:
        DEST.append(genesymbol_entrezgene_DICT[i])
    except:
        DEST.append(i)

HPRD_ENTREZGENE={'node1':SRC,'node2':DEST}
HPRD_ENTREZGENE = pd.DataFrame(HPRD_ENTREZGENE)
HPRD_ENTREZGENE.to_csv('ENTREZ/HPRD.txt', sep=' ', index=False)

# ####################################################################### STRING: #######################################################################

nodes_=[]
node_attributes=[]
for u, v in nx.get_node_attributes(STRING,'represents').items():
    symbols_=v.split('hgnc.symbol:')
    symbol=symbols_[1]
    nodes_.append(u)
    node_attributes.append(symbol)
nodes_attributes_DICT=dict(zip(nodes_, node_attributes))

edges_=STRING.edges
edges_list=list(edges_)
LIST=list(map(list, zip(*edges_list)))

src=[]
dest=[]
for i in LIST[0]:
    src.append(nodes_attributes_DICT[i])
for i in LIST[1]:
    dest.append(nodes_attributes_DICT[i])

STRING_GENE_SYMBOL={'node1':src,'node2':dest}
STRING_GENE_SYMBOL = pd.DataFrame(STRING_GENE_SYMBOL)
STRING_GENE_SYMBOL.to_csv('GENE_SYMBOL/STRING.txt', sep=' ', index=False)

# ----------------------------------------------------------------------

src_set=set(src)
dest_set=set(dest)
nodes_set=src_set.union(dest_set)
nodes_=list(nodes_set)

try:
    mg = mygene.MyGeneInfo()
    # out = mg.querymany(nodes_, scopes= 'symbol', fields='entrezgene', species='human', verbose=False)
    out = mg.querymany(nodes_, scopes= 'symbol', fields='uniprot', species='human', verbose=False)
except:
    pass

nodes_genesymbol=[]
nodes_uniprot=[]
for i in range(len(out)):
    try:
        nodes_uniprot.append(out[i]['uniprot'])
        nodes_genesymbol.append(out[i]['query'])
    except:
        pass

NODES_genesymbol=[]

NODES_uniprot=[]

for i in range(len(nodes_uniprot)):
    try:
        NODES_uniprot.append(nodes_uniprot[i]['Swiss-Prot'])
        NODES_genesymbol.append(nodes_genesymbol[i])
    except:
        pass

genesymbol_uniprot_DICT=dict(zip(NODES_genesymbol, NODES_uniprot))

SRC=[]
DEST=[]

for i in src:
    try:
        SRC.append(genesymbol_uniprot_DICT[i])
    except:
        SRC.append(i)

for i in dest:
    try:
        DEST.append(genesymbol_uniprot_DICT[i])
    except:
        DEST.append(i)

STRING_UNIPROT={'node1':SRC,'node2':DEST}
STRING_UNIPROT = pd.DataFrame(STRING_UNIPROT)
STRING_UNIPROT.to_csv('UNIPROT/STRING.txt', sep=' ', index=False)

# ----------------------------------------------------------------------

src_set=set(src)
dest_set=set(dest)
nodes_set=src_set.union(dest_set)
nodes_=list(nodes_set)

try:
    mg = mygene.MyGeneInfo()
    out = mg.querymany(nodes_, scopes= 'symbol', fields='entrezgene', species='human', verbose=False)
except:
    pass


NODES_genesymbol=[]
NODES_entrezgene=[]
for i in range(len(out)):
    try:
        NODES_entrezgene.append(out[i]['entrezgene'])
        NODES_genesymbol.append(out[i]['query'])
    except:
        pass

genesymbol_entrezgene_DICT=dict(zip(NODES_genesymbol, NODES_entrezgene))

SRC=[]
DEST=[]

for i in src:
    try:
        SRC.append(genesymbol_entrezgene_DICT[i])
    except:
        SRC.append(i)

for i in dest:
    try:
        DEST.append(genesymbol_entrezgene_DICT[i])
    except:
        DEST.append(i)

STRING_ENTREZGENE={'node1':SRC,'node2':DEST}
STRING_ENTREZGENE = pd.DataFrame(STRING_ENTREZGENE)
STRING_ENTREZGENE.to_csv('ENTREZ/STRING.txt', sep=' ', index=False)

############################################################ HIPPIE: ##################################################################################

src=[]
dest=[]
for u, v in nx.get_edge_attributes(HIPPIE,'name').items():
    nodes_=v.split(' (interacts-with) ')
    src.append(nodes_[0])
    dest.append(nodes_[1])
HIPPIE_SYMBOL = {'node1':src,'node2':dest}
HIPPIE_SYMBOL = pd.DataFrame(HIPPIE_SYMBOL)
HIPPIE_SYMBOL.to_csv('GENE_SYMBOL/HIPPIE.txt', sep=' ', index=False)

# ------------------------------------------------------------------------

src_set=set(src)
dest_set=set(dest)
nodes_set=src_set.union(dest_set)
nodes_=list(nodes_set)

try:
    mg = mygene.MyGeneInfo()
    # out = mg.querymany(nodes_, scopes= 'symbol', fields='entrezgene', species='human', verbose=False)
    out = mg.querymany(nodes_, scopes= 'symbol', fields='uniprot', species='human', verbose=False)
except:
    pass

NODES_genesymbol=[]
NODES_uniprot=[]
for i in range(len(out)):
    try:
        NODES_uniprot.append(out[i]['uniprot']['Swiss-Prot'])
        NODES_genesymbol.append(out[i]['query'])
    except:
        pass

genesymbol_uniprot_DICT=dict(zip(NODES_genesymbol, NODES_uniprot))

SRC=[]
DEST=[]

for i in src:
    try:
        SRC.append(genesymbol_uniprot_DICT[i])
    except:
        SRC.append(i)

for i in dest:
    try:
        DEST.append(genesymbol_uniprot_DICT[i])
    except:
        DEST.append(i)

HIPPIE_UNIPROT={'node1':SRC,'node2':DEST}
HIPPIE_UNIPROT = pd.DataFrame(HIPPIE_UNIPROT)
HIPPIE_UNIPROT.to_csv('UNIPROT/HIPPIE.txt', sep=' ', index=False)

# ----------------------------------------------------------------------

src_set=set(src)
dest_set=set(dest)
nodes_set=src_set.union(dest_set)
nodes_=list(nodes_set)

try:
    mg = mygene.MyGeneInfo()
    out = mg.querymany(nodes_, scopes= 'symbol', fields='entrezgene', species='human', verbose=False)
except:
    pass


NODES_genesymbol=[]
NODES_entrezgene=[]
for i in range(len(out)):
    try:
        NODES_entrezgene.append(out[i]['entrezgene'])
        NODES_genesymbol.append(out[i]['query'])
    except:
        pass

genesymbol_entrezgene_DICT=dict(zip(NODES_genesymbol, NODES_entrezgene))

SRC=[]
DEST=[]

for i in src:
    try:
        SRC.append(genesymbol_entrezgene_DICT[i])
    except:
        SRC.append(i)

for i in dest:
    try:
        DEST.append(genesymbol_entrezgene_DICT[i])
    except:
        DEST.append(i)

HIPPIE_ENTREZGENE={'node1':SRC,'node2':DEST}
HIPPIE_ENTREZGENE = pd.DataFrame(HIPPIE_ENTREZGENE)
HIPPIE_ENTREZGENE.to_csv('ENTREZ/HIPPIE.txt', sep=' ', index=False)

################################################################### HuRI: ##################################################

src=[]
dest=[]
for u, v in nx.get_edge_attributes(HuRI,'name').items():
    nodes_=v.split(' (interacts with) ')
    src_ensemble_id=nodes_[0].split('ensembl:')
    src_ensemble_id=src_ensemble_id[1]
    src.append(src_ensemble_id)
    dest_ensemble_id=nodes_[1].split('ensembl:')
    dest_ensemble_id=dest_ensemble_id[1]
    dest.append(dest_ensemble_id)

src_set=set(src)
dest_set=set(dest)
nodes_set=src_set.union(dest_set)
nodes_=list(nodes_set)

try:
    mg = mygene.MyGeneInfo()
    # out = mg.querymany(nodes_, scopes= 'symbol', fields='entrezgene', species='human', verbose=False)
    out = mg.querymany(nodes_, scopes= 'ensembl.gene', fields='symbol', species='human', verbose=False)
except:
    pass

NODES_ensembl=[]
NODES_genesymbol=[]
for i in range(len(out)):
    try:
        NODES_genesymbol.append(out[i]['symbol'])
        NODES_ensembl.append(out[i]['query'])
    except:
        pass

ensembl_genesymbol_DICT=dict(zip(NODES_ensembl, NODES_genesymbol))

SRC=[]
DEST=[]

for i in src:
    try:
        SRC.append(ensembl_genesymbol_DICT[i])
    except:
        SRC.append(i)

for i in dest:
    try:
        DEST.append(ensembl_genesymbol_DICT[i])
    except:
        DEST.append(i)

HuRI_GENE_SYMBOL={'node1':SRC,'node2':DEST}
HuRI_GENE_SYMBOL = pd.DataFrame(HuRI_GENE_SYMBOL)
HuRI_GENE_SYMBOL.to_csv('GENE_SYMBOL/HuRI.txt', sep=' ', index=False)

# # ------------------------------------------------------------------------

src_set=set(src)
dest_set=set(dest)
nodes_set=src_set.union(dest_set)
nodes_=list(nodes_set)

try:
    mg = mygene.MyGeneInfo()
    # out = mg.querymany(nodes_, scopes= 'symbol', fields='entrezgene', species='human', verbose=False)
    out = mg.querymany(nodes_, scopes= 'ensembl.gene', fields='uniprot', species='human', verbose=False)
except:
    pass

NODES_ensembl=[]
NODES_uniprot=[]
for i in range(len(out)):
    try:
        NODES_uniprot.append(out[i]['uniprot']['Swiss-Prot'])
        NODES_ensembl.append(out[i]['query'])
    except:
        pass

ensembl_uniprot_DICT=dict(zip(NODES_ensembl, NODES_uniprot))

SRC=[]
DEST=[]

for i in src:
    try:
        SRC.append(ensembl_uniprot_DICT[i])
    except:
        SRC.append(i)

for i in dest:
    try:
        DEST.append(ensembl_uniprot_DICT[i])
    except:
        DEST.append(i)

HuRI_UNIPROT={'node1':SRC,'node2':DEST}
HuRI_UNIPROT = pd.DataFrame(HuRI_UNIPROT)
HuRI_UNIPROT.to_csv('UNIPROT/HuRI.txt', sep=' ', index=False)

# # ----------------------------------------------------------------------

src_set=set(src)
dest_set=set(dest)
nodes_set=src_set.union(dest_set)
nodes_=list(nodes_set)

try:
    mg = mygene.MyGeneInfo()
    # out = mg.querymany(nodes_, scopes= 'symbol', fields='entrezgene', species='human', verbose=False)
    out = mg.querymany(nodes_, scopes= 'ensembl.gene', fields='entrezgene', species='human', verbose=False)
except:
    pass

NODES_ensembl=[]
NODES_entrezgene=[]
for i in range(len(out)):
    try:
        NODES_entrezgene.append(out[i]['entrezgene'])
        NODES_ensembl.append(out[i]['query'])
    except:
        pass

ensembl_entrezgene_DICT=dict(zip(NODES_ensembl, NODES_entrezgene))

SRC=[]
DEST=[]

for i in src:
    try:
        SRC.append(ensembl_entrezgene_DICT[i])
    except:
        SRC.append(i)

for i in dest:
    try:
        DEST.append(ensembl_entrezgene_DICT[i])
    except:
        DEST.append(i)

HuRI_ENTREZGENE={'node1':SRC,'node2':DEST}
HuRI_ENTREZGENE = pd.DataFrame(HuRI_ENTREZGENE)
HuRI_ENTREZGENE.to_csv('ENTREZ/HuRI.txt', sep=' ', index=False)




