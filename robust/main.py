from .ppi import *
from .steinerdiv import ExpMinMaxDiverseSteinerTreeComputer
import networkx as nx
import os.path
import warnings
import pandas as pd

def run(seeds, network='BioGRID', namespace='GENE_SYMBOL', alpha=0.25, beta=0.9, n=30, tau=0.1, study_bias_scores=None, gamma=1.0, outfile=None):
    is_graphml=0
    # Check the namespace and parse the network.
    namespace = _check_namespace(namespace)
    if network in ['APID', 'BioGRID', 'HPRD', 'STRING']:
        network = f'data/networks/{namespace}/{network}.txt'
    else:
        network,is_graphml=_check_and_preprocess_network(network, namespace)
    network = read_ppi_network(network, is_graphml)

    # Parse the study bias scores if provided and set the edge weights.
    if study_bias_scores is None:
        edge_weights = UnitEdgeWeight()
    else:
        path_to_study_bias_scores = _get_path_to_study_bias_scores(study_bias_scores, namespace)
        add_study_bias_scores_to_network(path_to_study_bias_scores, network)
        gamma = _check_gamma(gamma)
        edge_weights = BiasAwareEdgeWeight(network, gamma)

    # print(path_to_study_bias_scores)
    

    # Kick out seeds not in network.
    terminals = _get_terminals(seeds)
    terminals = list(set(terminals).intersection(set(network.nodes)))

    # Set up instance and engine.
    ppi_instance = PPIInstance(network, terminals, edge_weights)
    engine = ExpMinMaxDiverseSteinerTreeComputer(initial_fraction=alpha, reduction_factor=beta)

    # Compute the module.
    steiner_trees = engine(ppi_instance, n=n)
    module_as_df = steiner_trees.get_occurrences(include_terminals=True)
    module_as_df = module_as_df[module_as_df["%occurrences"] >= tau]

    # module_as_df.to_csv('module_as_df.txt', index=False)

    module_as_subgraph = steiner_trees.get_subgraph(threshold=tau)

    # Add connected component ID to nodes contained in module.
    comp_idx = 0
    for comp in sorted(nx.connected_components(module_as_subgraph), key=len, reverse=True):
        for node in comp:
            module_as_subgraph.nodes[node]['connected_components_id'] = comp_idx
        comp_idx += 1

    # Save module if requested.
    if outfile is not None:
        _save_module(module_as_df, module_as_subgraph, outfile)


    

    # Return module.
    return module_as_df, module_as_subgraph


def _check_namespace(namespace):
    if namespace not in ['GENE_SYMBOL', 'ENTREZ', 'UNIPROT']:
        warnings.warn(f'Illegal value {namespace} for parameter "namespace".\n'
                      f'==> Setting parameter "namespace" to "GENE_SYMBOL"')
        namespace = 'GENE_SYMBOL'
    return namespace


def _get_path_to_study_bias_scores(study_bias_scores, namespace):
    if os.path.exists(study_bias_scores):
        return study_bias_scores
    if study_bias_scores not in ['BAIT_USAGE', 'STUDY_ATTENTION']:
        warnings.warn(f'Illegal value {study_bias_scores} for parameter "study_bias_scores".\n'
                      f'==> Setting parameter "study_bias_scores" to "BAIT_USAGE"')
        study_bias_scores = 'BAIT_USAGE'
    return f'data/study_bias_scores/{namespace}/{study_bias_scores}.csv'


def _check_and_preprocess_network(network, namespace):
    is_graphml=0
    if type(network) is nx.Graph:
        is_graphml=1
    elif network.endswith('.graphml'):
        is_graphml=1
        network=nx.read_graphml(network)
    elif isinstance(network, pd.DataFrame):
        network.to_csv('data/networks/{namespace}_customNetwork.txt')
        network='data/networks/{namespace}_customNetwork.txt'
    else:
        if not (network.endswith('.txt') or network.endswith('.csv') or network.endswith('.tsv')):
            raise ValueError(f'Illegal network type: {network}')
    return network, is_graphml


def _check_gamma(gamma):
    if gamma > 1:
        warnings.warn(f'Illegal value {gamma} > 1 for parameter "gamma".\n'
                      f'==> Setting parameter "gamma" to 1.0')
        gamma = 1.0
    if gamma < 0:
        warnings.warn(f'Illegal value {gamma} < 0 for parameter "gamma".\n'
                      f'==> Setting parameter "gamma" to 0.0')
        gamma = 0.0
    return gamma


def _get_terminals(seeds):
    if os.path.exists(seeds):
       seeds = read_terminals(seeds)
    # if not isinstance(seeds, list):
    #     raise ValueError(f'Illegal value {seeds} for parameter "seeds".\n'
    #                      f'Must be a path or a list of gene/protein identifiers.')
    return seeds


def _save_module(solution_as_df, subgraph, outfile):
    if outfile.endswith(".csv"):
        solution_as_df.to_csv(outfile)
    elif outfile.endswith(".graphml"):
        nx.write_graphml(subgraph, outfile)
    else:
        nx.write_edgelist(subgraph, outfile, data=False)