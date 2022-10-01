# Installation

Install conda environment as follows (there also exists an environment.yml file, but it contains more packages than necessary)
```bash
conda create --name biosteiner python=3.7
conda activate biosteiner
conda install numpy matplotlib pandas networkx pip jupyter
pip install pcst_fast
```

# Running ROBUST

You can simply run ROBUST by calling
```bash
python3 robust/robust.py robust/data/seeds/GENE_SYMBOL/30_ms_seeds.txt
```
The positional arguments are:
```
[1] file containing the list of seed genes/ proteins
[2] path to output file
```

The optional arguments are:
```
--network NETWORK	Description: specify input network to run ROBUST on, type: string, default: 'BioGRID'. Options include: in-built network ('BioGRID', 'APID', 'HPRD', 'STRING') or file_path (.graphml, .txt, .csv, .tsv)
--namespace {'ENTREZ', 'GENE_SYMBOL', 'UNIPROT'}		Description: gene/ protein identifier options for network and study bias data, type=str, default: 'GENE_SYMBOL'
--alpha INITIAL_FRACTION							Description: initial fraction for ROBUST, type=float, expected range=[0,1], default: 0.25
--beta REDUCTION_FACTOR							Description: reduction factor for ROBUST, type=float, expected range=[0,1], default: 0.90
--n NO_OF_STEINER_TREES						Description: # of steiner trees for ROBUST, type=int, expected range=(0,+inf], default: 30
--tau THRESHOLD									Description: threshold value for ROBUST, type=float, expected range=(0,1], default: 0.1
--study-bias-scores							Description: edge weight function used by ROBUST, type=str, default: 'NONE'. Options include: in-built options ('NONE'/'None', 'BAIT_USAGE', 'STUDY_ATTENTION') or file_path (.graphml, .txt, .csv, .tsv)
--gamma										Description: study bias data regulator, type=float, expected range=[0,1], default: 1.00
```

The suffix of the path to the output file you specify, determine the format of the output.
You can either choose
- .graphml: A .graphml file is written that contains the following vertex properties: isSeed, significance, nrOfOccurrences, connected_components_id, trees
- .csv: A .csv file which contains a vertex table with #occurrences, %occurrences, terminal (isSeed) 
- everything else: An edge list

# Evaluating ROBUST-2.00

For a large-scale empirical evaluation of ROBUST-2.00, please follow the instructions given here: https://github.com/bionetslab/robust-eval.

# Citing ROBUST-2.00

Please cite ROBUST-2.00 as follows:
- J. Bernett, D. Krupke, S. Sadegh1, J. Baumbach, S. P. Fekete, T. Kacprowski, M. List1, D. B. Blumenthal: Robust disease module mining via enumeration of diverse prize-collecting Steiner trees, *Bioinformatics* 38(6), pp. 1600-1606, https://doi.org/10.1093/bioinformatics/btab876.
