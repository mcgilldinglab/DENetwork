# DENetwork
Unveiling Regulatory and Signaling Networks Behind Differentially Expressed Genes

  * [About](#about)
  * [Required Packages](#required-packages)
  * [Installation](#installation)
  * [Usage](#usage)
  * [Results](#results)
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

## Required Packages

Python 3 is required.

## Installation

1. Clone or download this repository.

2. Navigate to the downloaded DENetwork folder and download required Python packages

```bash
$ cd DENetwork
$ python3 setup.py install
```

## Usage (FIX INDENTATION...IS REALLY MESSY RN)

1. Obtain the following data files:
  * RNA-seq gene expression data for WT and disease samples
  * Receptors specific to your disease condition

The data should be formatted as shown in the Example section. 

2. Run DESeq2

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
                        comma- (.csv) or tab-delimited (.tsv) file of raw RNA-
                        seq gene expression data counts that contain samples
                        of the wildtype and disease conditions (if includes
                        other conditions, they will be filtered out)

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
                        whether to use the unadjusted (pvalue) or adjusted
                        (padj) p-values to filter DESeq2 result matrix.

  -p P_VALUE_THRESHOLD, --p_value_threshold P_VALUE_THRESHOLD, optional
                        Float, Optional, is 0.05 by default. Filters the
                        DESeq2 result matrix to remove insignificant genes.

  -l LOG2_FOLD_CHANGE_THRESHOLD, --log2_fold_change_threshold LOG2_FOLD_CHANGE_THRESHOLD, optional
                        Positive Float, Optional, is 0.6 by default for
                        upregulated DE genes and -0.6 for downregulated DE
                        genes. Filters DE genes.
```

Aftering running DESeq2, you should get the following files (to use as input for the next step):
  1. 3 files containing differentially-expressed (DE) genes
    * DE_pos_*.tsv containing all upregulated DE genes
    * DE_neg_*.tsv containing all downregulated DE genes
    * DE_all_*.tsv containing all DE genes
  2. a file containing all genes and their log2 fold changes (gene_log2fc.tsv)

3. Run DENetwork

```bash
usage: run_denetwork.py [-h] -n NAME -d DE_FILE -g GENE_FILE -r RECEP_FILE -t
                        {de,tf}

optional arguments:
  -h, --help            show this help message and exit

  -n NAME, --name NAME, required
                        name of the DENetwork model

  -d DE_FILE, --de_file DE_FILE, required
                        file of differentially-expressed genes

  -g GENE_FILE, --gene_file GENE_FILE, required
                        file containing all genes, their (output result matrix
                        from DESeq2)

  -r RECEP_FILE, --recep_file RECEP_FILE, required
                        file of disease-specific receptors

  -t {de,tf}, --targets {de,tf}, required
                        choose whether the targets are differentially-
                        expressed genes (de) OR transcription factors (tf)
```

There are many model parameters in DENetwork, so they are not listed as optional arguments when running run_denetwork.py. They can be manually changed in the run_denetwork.py file. Otherwise, the default model parameters are used.

## Results

## Example

## License
DENetwork is licensed under the terms of the MIT license.

## Credits
DENetwork was developed by the Ding Lab @ McGill University, and implemented by Ting-Yi Su.

## Contact
Jun Ding jun.ding@mcgill.ca <br />
Ting-Yi Su ting-yi.su@mail.mcgill.ca

## Citation
