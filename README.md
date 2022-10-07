# Installation

Install conda environment as follows (there also exists an environment.yml but it contains more packages than necessary)
```bash
conda create --name robust python=3.7
conda activate robust
conda install numpy matplotlib pandas networkx pip jupyter
pip install pcst_fast
```

# Running ROBUST

You can simply run robust by calling
```bash
python robust.py data/data-example1-prec-puberty/BioGRID.txt data/data-example1-prec-puberty/prec-pub-seeds.txt prec_puberty.graphml
```
The positional arguments are:
```
[1] file with a list of seed genes (separator: newline-separated)
[2] path to output file (supported output file types: .graphml, .csv, others) [read more below]


The suffix of the path to the output file you specify, determine the format of the output.
You can either choose
- .graphml: A .graphml file is written that contains the following vertex properties: isSeed, significance, nrOfOccurrences, connected_components_id, trees
- .csv: A .csv file which contains a vertex table with #occurrences, %occurrences, terminal (isSeed) 
- everything else: An edge list
```



The optional arguments are:
```
--initial_fraction INITIAL_FRACTION							Description: initial fraction for ROBUST, type=float, expected range=[0,1], default: 0.25
--reduction_factor REDUCTION_FACTOR							Description: reduction factor for ROBUST, type=float, expected range=[0,1], default: 0.90
--number_of_steiner_trees NO_OF_STEINER_TREES						Description: # of steiner trees for ROBUST, type=int, expected range=(0,+inf], default: 30
--threshold THRESHOLD									Description: threshold value for ROBUST, type=float, expected range=(0,+inf], default: 0.1
--node_namespace {'ENTREZ_GENE_ID', 'GENE_SYMBOL', 'UNIPROT_PROTEIN_ID'}		Description: gene/ protein identifier options for study bias data, type=str, default: 'GENE_SYMBOL'
--edge_cost {'UNIFORM', 'ADDITIVE', 'EXPONENTIAL'}					Description: function for calculating edge costs, type=str, default: UNIFORM
--normalize {'BAIT_USAGE', 'STUDY_ATTENTION', 'CUSTOM'}					Description: study bias data options to be used for normalization, type=str, default: 'BAIT_USAGE'
--lambda										Description: lambda value for ROBUST-version-2.00, type=float, expected range=[0,1], default: 0.50




[1] file providing the network:
	
	Input options:
	- A two-column edgelist. File types and corresponding separators are as follows: 1. '.txt' file should be space-separated 2. '.tsv' file should be tab-separated 3. '.csv' file should be comma-separated. No other file  formats except '.txt', '.csv' and '.tsv' are accepted at the moment.
	- A valid .graphml file
	- In-built network name {'BioGRID', 'APID', 'HPRD', 'STRING'}




```


# Evaluating ROBUST

For a large-scale empirical evaluation of ROBUST, please follow the instructions given here: https://github.com/bionetslab/robust-eval.

# Citing ROBUST

Please cite ROBUST as follows:
- **citation will be added once available**
