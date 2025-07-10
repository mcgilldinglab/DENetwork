# DENetwork
Unveiling regulatory and signaling networks behind differentially-expressed genes

  * [About](#about)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Example](#example)
  * [License](#license)
  * [Credits](#credits)
  * [Contact](#contact)
  * [Citation](#citation)
  * [RNA-seq datasets](#rna-seq-datasets)


## About

Differential gene expression analysis of RNA-sequencing (RNA-seq) data offers crucial insights into biological differences between sample groups. However, the conventional focus on differentially-expressed (DE) genes often omits non-DE regulators, which are an integral part of such differences. Moreover, DE genes frequently serve as passive indicators of transcriptomic variations rather than active influencers, limiting their utility as intervention targets. To address these shortcomings, we have developed DENetwork. This innovative approach deciphers the intricate regulatory and signaling networks driving transcriptomic variations between conditions with distinct phenotypes. Unique in its integration of both DE and critical non-DE genes in a graphical model, DENetwork enhances the capabilities of traditional differential gene expression analysis tools, such as DESeq2. Our application of DENetwork to an array of simulated and real datasets showcases its potential to encapsulate biological differences, as demonstrated by the relevance and statistical significance of enriched gene functional terms. DENetwork offers a robust platform for systematically characterizing the biological mechanisms that underpin phenotypic differences, thereby augmenting our understanding of biological variations and facilitating the formulation of effective intervention strategies.

<p align="center"> 
  <img src="https://github.com/mcgilldinglab/DENetwork/blob/main/images/github_denetwork.svg" />
</p> 


## Installation

Python 3 is required.

1. Clone or download this repository.

2. Navigate to the downloaded DENetwork folder and install required Python packages.

```bash
$ cd DENetwork
$ python3 setup.py install
```
3. We use OmicVerse to perform differential analysis (DESeq2), please follow [OmicVerse documentation](https://omicverse.readthedocs.io/en/latest/index.html).

You can skip this step if you don't need to perform differential gene expression analysis, or are using another method (e.g. edgeR).

## Usage

**1. Obtain the following data files:**

  * RNA-seq gene expression count data for wildtype (healthy) and disease samples (required to run DESeq2)
  * Receptors specific to your disease condition (required to run DENetwork)

  The data should be formatted as shown in the [Example](#example) section. 

**2. Run DESeq2**

Please perform differential analysis by following the OmicVerse tutorial: [Bulk RNA-seq Differential Analysis with DESeq2](https://omicverse.readthedocs.io/en/latest/Tutorials-bulk/t_deseq2/).

You can also choose to use other methods for differential gene expression analysis (e.g. edgeR), instead of DESeq2. If you choose to do so, please make sure that your input files for the next step are as indicated in the [Example](#example) section.

**3. Run DENetwork**

```bash
usage: run_denetwork.py [-h] -n NAME -d DE_FILE -g GENE_FILE -r RECEP_FILE -t
                        {de,tf} -s NUM_SIGNIFICANT_GENES

optional arguments:
  -h, --help            show this help message and exit

  -n NAME, --name NAME, required  
                        name of the DENetwork model

  -d DE_FILE, --de_file DE_FILE, required 
                        tab-delimited file (.tsv) of differentially-expressed
                        genes, along with their baseMean, log2FoldChange,
                        lfcSE, stat, pvalue, and padj values

  -g GENE_FILE, --gene_file GENE_FILE, required 
                        tab-delimited file (.tsv) containing all genes and
                        their |log2FoldChange|

  -r RECEP_FILE, --recep_file RECEP_FILE, required 
                        text file (.txt) of disease-specific receptors

  -t {de,tf}, --targets {de,tf}, required 
                        choose whether the targets are differentially-
                        expressed genes (de) OR transcription factors (tf)

  -s NUM_SIGNIFICANT_GENES, --num_significant_genes NUM_SIGNIFICANT_GENES, required 
                        number of top-ranked genes to select from the (near)
                        optimal network; these selected genes make up the
                        final network (signaling/regulatory pathway) and can
                        be used for GO biological process & reactome pathway
                        analyses

  -G GAMMA, --gamma GAMMA, optional
                        Float, Optional, is 1.0 by default. Used to change how
                        much the log2 fold changes (node attributes) influence
                        the path score.

  -K MAX_NUM_INTERMEDIATE_NODES_IN_SHORTEST_PATHS, --max_num_intermediate_nodes_in_shortest_paths MAX_NUM_INTERMEDIATE_NODES_IN_SHORTEST_PATHS, optional
                        Integer, Optional, is 5 by default. Parameter for
                        finding shortest paths. Is the max number of
                        intermediate nodes allowed when finding shortest paths
                        between source and target pairs.

  -p MAX_P_VALUE_OF_SIGNIFICANT_PATHS, --max_p_value_of_significant_paths MAX_P_VALUE_OF_SIGNIFICANT_PATHS, optional
                        Float, Optional, is 0.01 by default. Is used as a cut-
                        off for finding significant shortest paths.

  -N NUM_TOP_PATHS, --num_top_paths NUM_TOP_PATHS, optional
                        Integer, Optional, is 5 by default. Parameter for
                        number of top significant shortest paths kept between
                        each source and target pair.

  -v PERCENT_IMPROVEMENT_FOR_CONVERGENCE, --percent_improvement_for_convergence PERCENT_IMPROVEMENT_FOR_CONVERGENCE, optional
                        Float, Optional, is 5.0 by default. Percent
                        improvement required to determine convergence of
                        global optimal solution (& to terminate iterations).
                        
  -T NUM_TRYS, --num_trys NUM_TRYS, optional
                        Int, Optional, is 5 by default. Number of trys allowed
                        for continuous no improvements before terminating &
                        determining the converged solution.
```

DENetwork should take a couple minutes to run.
There are many model parameters in DENetwork, so they are not listed as optional arguments in run_denetwork.py. They can be manually changed in the run_denetwork.py file. Otherwise, the default model parameters are used.


## Example

The 'example' folder contains an example of RNA-seq gene expression data (GSE192528_RawCountsAnnotated.xlsx) and disease-specific receptors (receptors_influenza.txt) for an IAV (influenza) dataset.

Your RNA-seq data should have the following columns:
* column of genes, with 'Gene' as the column name
* a column of gene counts for each sample, with the disease/wildtype condition name in the column name (e.g. A22_influenza as the column name for an influenza sample)

```bash
$ python3 run_deseq2.py -n influenza -r example/GSE192528_RawCountsAnnotated.xlsx -w uninfected -d influenza -o example_deseq2_output
```

After running DESeq2, you will end up with the following files, in the 'example_deseq2_output' folder, that you will need to run DENetwork:
* 3 .tsv files containing differentially-expressed (DE) genes starting with 'DE'. 
    1. 'DE_all_influenza.tsv' contains all DE genes
    2. 'DE_pos_influenza.tsv' contains all upregulated DE genes
    3. 'DE_neg_influenza.tsv' contains all downregulated DE genes. 

  You can choose either of the 3 files for input into DENetwork. Ideally, you want to obtain tens to hundreds (~30-200) DE genes for input into DENetwork. For example, DE_pos_influenza.tsv is chosen and it contains 34 upregulated DE genes. The optional input arguments (-b, -t, -p, -l) can be changed to increase or decrease the number of DE genes outputted. 

* 'gene_log2fc_influenza.tsv' which contains a gene column and a log2 fold change (in absolute value) column.

The 'data' folder contains 2 files that are mandatory for running DENetwork:
* 'ppi_ptm_pd.hgnc.txt' contains a list of PPI scores (both protein-protein and protein-dna) 
* 'TFs.txt' contains a list of all TFs that are used to find significant TFs, if TFs are chosen as the target nodes

```bash
$ python3 run_denetwork.py -n influenza -d example_deseq2_output/DE_pos_influenza.tsv -g example_deseq2_output/gene_log2fc_influenza.tsv -r example/receptors_influenza.txt -t de -s 100
```

You can choose either DE genes or transcription factors (TFs) as the target nodes in DENetwork. Here, DE genes are chosen using the '-t de' option. You can also choose to use either all, only upregulated, or only downregulated DE genes. Here, upregulated (pos) DE genes are chosen using the '-d example_deseq2_output/DE_pos_influenza.tsv' option. 

The output files of DENetwork are in the 'figures', 'files', and 'results' folders. Each DENetwork model has its own subdirectory within each of these 3 folders. In this example, the subdirectory is named 'influenza' (using the '-n influenza' option when running DENetwork).


* 'results/influenza' contains the output of DENetwork:
    * 'genes_go.txt' contains all 100 top-ranked genes ('-s 100' option)
    * 'receptors_go.txt' contains all receptors amongst the top genes (targets that are also receptors are listed as receptors)
    * 'targets_go.txt' contains all targets amongst the top genes
    * 'internal_nodes.txt' contains all intermediate signaling/regulatory genes amonst the top genes
    * 'optimal_network.sif' contains the connections between the top-ranked genes
    * 'optimal_network.noa' contains the attributes (type and node size) of the top-ranked genes
    * 'nodes_ranking.tsv' contains the ranking list of the top genes as well as their score changes and node type

    The .sif and .noa files can be imported into Cytoscape to view the final signaling/regulatory network. The genes in the .txt files can be used for GO term enrichment and pathway analyses.

* 'figures/influenza' contains the following folders:
    * 'path_scores_distribution' contains distribution plots of the path scores 
    * 'S_Q' contains S(Q) and score penalty plots (from equation 5 in the manuscript)

    Figure names without a number at the end were plotted when finding a local optimal solution, and those with a number at the end were plotted when finding the (near) optimal solution. 

    If there is more than one improvement during the iterative refinement process, there will be an additional improvement plot named 'step5_improvement_8_5_5.svg' (where l=8, N=5, t=5 in the optimal_graph function in denetwork.py)
  
* 'files/influenza' contains the following files:
    * 'local_op_g.pickle' pickle file containing the local optimal graph object
    * 'optimal_graph_8_5_5.pickle' pickle file containing the (near) optimal graph object
    * 'shortest_paths.pickle' pickle file containing the shortest paths found in the initial fully-connected network
    * 'nodes_ranking_score_differences_8_5_5.txt' and 'nodes_ranking_scores_8_5_5.txt' contain the score differences and scores of node-removals in the (near) optimal network

## License
DENetwork is licensed under the terms of the MIT license.

## Credits
DENetwork was developed by the Ding Lab @ McGill University, and implemented by Ting-Yi Su.

## Contact
Jun Ding jun.ding@mcgill.ca <br />
Ting-Yi Su ting-yi.su@mail.mcgill.ca

## Citation

```
@article {Su2023,
  author = {Ting-Yi Su and Quazi S. Islam and Steven K. Huang and Carolyn J. Baglole and Jun Ding},
  title = {DENetwork: Unveiling Regulatory and Signaling Networks Behind Differentially-Expressed Genes},
  elocation-id = {2023.06.27.546719},
  year = {2023},
  doi = {10.1101/2023.06.27.546719},
  publisher = {Cold Spring Harbor Laboratory},
  URL = {https://www.biorxiv.org/content/early/2023/06/27/2023.06.27.546719},
  eprint = {https://www.biorxiv.org/content/early/2023/06/27/2023.06.27.546719.full.pdf},
  journal = {bioRxiv}
}
```

## RNA-seq datasets

The other 2 RNA-seq datasets used in the manuscript are available under the 'rna_seq_data' folder.
  * 'Alox15_macrophage_rna_seq_data.tsv': data for wildtype and Alox15-/- mouse macrophages
  * 'covid_rna_seq_data.tsv': data for patients tested positive or negative for the SARS-CoV-2 virus
