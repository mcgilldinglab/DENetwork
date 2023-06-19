# DENetwork
Unveiling Regulatory and Signaling Networks Behind Differentially Expressed Genes

  * [About](#about)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Example](#example)
  * [License](#license)
  * [Credits](#credits)
  * [Contact](#contact)
  * [Citation](#citation)


## About

The identification of differentially expressed genes from RNA-seq data is instrumental for understanding the biological differences between sample groups. However, the focus on differential genes can overlook the potential primary drivers of these differences and viable intervention targets. Additionally, key non-differential regulators may be missed by standard RNA-seq differential gene analysis. To mitigate these limitations, we present DENetwork, an innovative approach that deciphers the regulatory and signaling networks responsible for transcriptomic variations between conditions with distinct phenotypes. This methodology augments conventional differential gene analysis tools, such as DESeq2, by embracing both differential and non-differential critical genes. Through applying DENetwork to diverse simulated and real datasets, we exhibit its efficacy in encapsulating biological differences, denoted by the relevance and statistical significance of enriched gene functional terms. As a robust tool, DENetwork systematically characterizes the biological mechanisms driving phenotypic differences, thus deepening our understanding of biological variations between conditions and aiding in the design of potent intervention strategies.

<p align="center"> 
  <img src="https://github.com/mcgilldinglab/DENetwork/blob/main/images/github_denetwork.svg" />
</p> 


## Installation

Python 3 is required.

1. Clone or download this repository.

2. Navigate to the downloaded DENetwork folder and download required Python packages

```bash
$ cd DENetwork
$ python3 setup.py install
```

## Usage

**1. Obtain the following data files:**

  * RNA-seq gene expression data for WT and disease samples (required to run DESeq2)
  * Receptors specific to your disease condition (required to run DENetwork)

  The data should be formatted as shown in the [Example](#example) section. 

**2. Run DESeq2**

```bash
usage: run_deseq2.py [-h] -n NAME -r RAW_COUNTS -w WT_CONDITION -d
                     DISEASE_CONDITION -o OUTPUT_DIR [-b BASE_MEAN_THRESHOLD]
                     [-t {padj,pvalue}] [-p P_VALUE_THRESHOLD]
                     [-l LOG2_FOLD_CHANGE_THRESHOLD]

optional arguments:
  -h, --help            show this help message and exit

  -n NAME, --name NAME, required
                        name of disease condition, used to name the output
                        files

  -r RAW_COUNTS, --raw_counts RAW_COUNTS, required
                        comma- (.csv), tab-delimited (.tsv) or excel (.xlsx)
                        file of raw RNA-seq gene expression data counts that
                        contain samples of the wildtype and disease conditions
                        (if includes samples for other conditions, they will
                        be filtered out)

  -w WT_CONDITION, --wt_condition WT_CONDITION, required
                        wildtype (healthy) condition

  -d DISEASE_CONDITION, --disease_condition DISEASE_CONDITION, required
                        disease condition

  -o OUTPUT_DIR, --output_dir OUTPUT_DIR, required
                        output directory where output files are stored

  -b BASE_MEAN_THRESHOLD, --base_mean_threshold BASE_MEAN_THRESHOLD, optional
                        Integer, Optional, is 50 by default. Filters the
                        DESeq2 result matrix to remove genes with low
                        expression counts.

  -t {padj,pvalue}, --p_value_threshold_type {padj,pvalue}, optional
                        padj/pvalue, Optional, is padj by default. Choose
                        whether to use the unadjusted (pvalue) or FDR adjusted
                        (padj) p-values to filter DESeq2 result matrix.

  -p P_VALUE_THRESHOLD, --p_value_threshold P_VALUE_THRESHOLD, optional
                        Float, Optional, is 0.05 by default. Filters the
                        DESeq2 result matrix to remove insignificant genes.

  -l LOG2_FOLD_CHANGE_THRESHOLD, --log2_fold_change_threshold LOG2_FOLD_CHANGE_THRESHOLD, optional
                        Positive Float, Optional, is 0.6 by default for
                        upregulated DE genes and -0.6 for downregulated DE
                        genes. Filters DE genes.
```
If you already have a list of DE genes and do not need to run DESeq2, please make sure that your input files for the next step are as indicated in the [Example](#example) section.


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
                        lfcSE, stat, pvalue, and padj

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
                        analysis
```

There are many model parameters in DENetwork, so they are not listed as optional arguments in run_denetwork.py. They can be manually changed in the run_denetwork.py file. Otherwise, the default model parameters are used.


## Example

The 'example' folder contains an example of RNA-seq gene expression data (GSE192528_RawCountsAnnotated.xlsx) and disease-specific receptors (receptors_influenza.txt) for an IAV (influenza) dataset.

Your RNA-seq data should have the following columns:
* column of genes, with 'Gene' as the column name
* a column of gene counts for each sample, with the disease/wildtype condition name in the column name (e.g. A22_influenza as the column name for an influenza sample)

```bash
$python3 run_deseq2.py -n influenza -r example/GSE192528_RawCountsAnnotated.xlsx -w uninfected -d influenza -o example_deseq2_output -t pvalue
```

After running DESeq2, you will end up with the following files, in the 'example_deseq2_output' folder, that you will need to run DENetwork:
* 3 .tsv files containing differentially-expressed (DE) genes starting with 'DE'. 
    1. 'DE_all_influenza.tsv' contains all DE genes
    2. 'DE_pos_influenza.tsv' contains all upregulated DE genes
    3. 'DE_neg_influenza.tsv' contains all downregulated DE genes. 

  You can choose either of the 3 files for input into DENetwork. Ideally, you want to obtain a couple hundred (~200-400) DE genes for input into DENetwork. For example, DE_pos_influenza.tsv is chosen and it contains 236 upregulated DE genes. The optional input arguments (-b, -t, -p, -l) can be changed to increase or decrease the number of DE genes outputted. 

* 'gene_log2fc_influenza.tsv' which contains a gene column and a log2 fold change (in absolute value) column.

The 'data' folder contains 2 files that are mandatory for running DENetwork:
* 'ppi_ptm_pd.hgnc.txt' contains a list of PPI scores (both protein-protein and protein-dna) 
* 'TFs.txt' contains a list of all TFs that are used to find significant TFs, if TFs are chosen as the target nodes

```bash
$python3 run_denetwork.py -n influenza -d example_deseq2_output/DE_pos_influenza.tsv -g example_deseq2_output/gene_log2fc_influenza.tsv -r example/receptors_influenza.txt -t de -s 100
```

You can choose either DE genes or transcription factors (TFs) as the target nodes in DENetwork. Here, DE genes are chosen using the '-t de' option. You can also choose to use either all, only upregulated, or only downregulated DE genes. Here, upregulated (pos) DE genes are chosen using the '-d example_deseq2_output/DE_pos_influenza.tsv' option. 

The output files of DENetwork are in the 'figures', 'files', and 'results' folders. Each DENetwork model has its own subdirectory within each of these 3 folders. In this example, the subdirectory is named 'influenza' (using the '-n influenza' option when running DENetwork).


* 'results/influenza' contains the output of DENetwork:
    * 'genes_go.txt' contains all top-ranked genes kept using the '-s' argument
    * 'receptors_go.txt' contains all receptors amongst the top genes kept (targets that are also receptors, are listed as receptors)
    * 'targets_go.txt' contains all targets amongst the top genes kept
    * 'internal_nodes.txt' contains all intermediate signaling/regulatory genes amonst the top genes kept
    * 'optimal_network.sif' contains the connections between the top-ranked genes
    * 'optimal_network.noa' contains the attributes (type and node size) of the top-ranked genes
    * 'nodes_ranking.tsv' contains the ranking list of the top genes as well as their score changes and node type

    The .sif and .noa file can be imported into Cytoscape to view the final signaling/regulatory network. The genes in the .txt files can be used for GO term enrichment and pathway analyses.

* 'figures/influenza' contains the following folders:
    * 'path_scores_distribution' folder containing distribution plots of the path scores 
    * 'S_Q' folder containing S(Q) and score penalty plots

    Figure names without a number at the end were plotted when finding a local optimal solution, and those with a number at the end were plotted when finding the (near) optimal solution. 

    If there is more than one improvement during the iterative refinement process, there will be additional an improvement plot named 'step5_improvement_8_5_5.svg' (where l=8, N=5, t=5 in the optimal_graph function in denetwork.py)
  
* 'files/influenza' contains the following files:
    * 'local_op_g.pickle' pickle file containing the local optimal graph object
    * 'optimal_graph_8_5_5.pickle' pickle file containing the (near) optimal graph object
    * 'shortest_paths.pickle' pickle file containing the shortest paths found
    * 'nodes_ranking_score_differences_8_5_5.txt' and 'nodes_ranking_scores_8_5_5.txt' contain the score differences and scores of node-removals


## License
DENetwork is licensed under the terms of the MIT license.

## Credits
DENetwork was developed by the Ding Lab @ McGill University, and implemented by Ting-Yi Su.

## Contact
Jun Ding jun.ding@mcgill.ca <br />
Ting-Yi Su ting-yi.su@mail.mcgill.ca

## Citation
