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

## Usage

1. Obtain the following data files:
  * RNA-seq gene expression data for WT and disease samples
  * Receptors specific to your disease condition

The data should be formatted as shown in the Example section. 

2. Run DESeq2 (FIX TABBING)

```bash
$ python3 run_deseq2.py -c COUNTDATA -m METADATA (NOT DONE YET)
```

Aftering running DESeq2, you should get the following files (to use as input for the next step):
  1. 3 files containing differentially-expressed (DE) genes
    * DE_pos_*.tsv containing all upregulated DE genes
    * DE_neg_*.tsv containing all downregulated DE genes
    * DE_all_*.tsv containing all DE genes
  2. a file containing all genes and their log2 fold changes (gene_log2fc.tsv)

3. Run DENetwork

```bash
$ python3 run_deseq2.py -n NAME -d DEFILE -g GENEFILE -r RECEPFILE -t {de,tf}
```

```bash
usage: run_denetwork.py [-h] -n NAME -d DEFILE -g GENEFILE -r RECEPFILE -t {de,tf}

optional arguments:
  -h, --help            show this help message and exit
  -n NAME, --name NAME  name of the DENetwork model
  -d DEFILE, --defile DEFILE
                        file of differentially-expressed genes
  -g GENEFILE, --genefile GENEFILE
                        file containing all genes, their (output result matrix
                        from DESeq2)
  -r RECEPFILE, --recepfile RECEPFILE
                        file of disease-specific receptors
  -t {de,tf}, --targets {de,tf}
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
