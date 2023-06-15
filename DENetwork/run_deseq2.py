import pdb,sys,os
from multiprocessing import Pool
from File import *
import numpy as np
from BioUtils import BioList
import pandas as pd
import csv
import math
import pandas as pd
import statistics as stats
import os.path as osp

def runDESeq2(datMatrixFile,metaMatrixFile,A,B):
    from diffexpr.py_deseq import py_DESeq2
    datMat=TabFile(datMatrixFile).read("\t")
    cells=datMat[0][1:]
    genes=[item[0] for item in datMat[1:]]
    dat=[item[1:] for item in datMat[1:]]
    pdata=pd.DataFrame(data=np.array(dat,dtype='float'))
    pdata.columns=cells
    pdata['Gene']=genes
    
    metaMat=TabFile(metaMatrixFile).read("\t")
    pmeta=pd.DataFrame(data=metaMat[1:],columns=metaMat[0])
    pmeta.index=[item[0] for item in metaMat[1:]]
    
    dds = py_DESeq2(count_matrix = pdata,
               design_matrix = pmeta,
               design_formula = '~ condition',
               gene_column = 'Gene') # <- telling DESeq2 this should be the gene ID column
    dds.run_deseq()
    nMatrix=dds.normalized_count()
    dds.get_deseq_result(contrast=['condition',A,B])
    res=dds.deseq_result
    return [res,nMatrix]

# get gene_loc2fc here!
def find_log2fc(res_matrix_file, gene_log2fc_file):
    # get log2fc
    gene_log2fc = {}
    max_log2fc = 0.0
    with open(res_matrix_file, 'r') as res, open(gene_log2fc_file, 'w') as gene_w:
        next(res)
        for line in res:
            items = line.split('\t')
            gene = items[0].upper()
            log2fc = float(items[2])
            if log2fc < 0.0: log2fc = abs(log2fc)
            gene_log2fc[gene] = log2fc # may have multiple log2fc for same gene (if gene repeated in scRNA seq counts), get the last entry
            if log2fc > max_log2fc:
                max_log2fc = log2fc
            
        for gene in gene_log2fc:
            gene_w.write(gene + '\t' + str(gene_log2fc[gene]) + '\n')
    print('max log2fc:', max_log2fc)

def get_all_de_genes(raw_counts, name, wt_condition, disease_condition, wt_condition_name, disease_condition_name):
    print('Finding all DE genes for:', name)
    df = pd.read_csv(raw_counts, index_col=None)

    sample_names = list(df)[1:] # sample names
    column_names = []

    # files to write to & use
    meta_data_file = osp.join('data', name + '_meta_data.tsv')
    count_data_file = osp.join('data', name + '_count_data.tsv')
    res_matrix_file = osp.join('data', 'res_matrix_' + name + '.tsv')
    gene_log2fc_file = osp.join('data', 'gene_log2fc_' + name + '.tsv')
    de_genes_file = osp.join('data', 'DE_all_' + name + '.tsv')
    normalized_matrix_file = osp.join('data', name + '_normalized_matrix.tsv')

    # create meta data file
    with open(meta_data_file, 'w') as f:
        f.write('\t'.join(['sample', 'condition']) + '\n')
        column_names = ['Gene']
        for sample_name in sample_names: # remove gene name column
            condition = ''
            if wt_condition in sample_name:
                condition = wt_condition_name
                f.write('\t'.join([sample_name, condition]) + '\n')
                column_names.append(sample_name)
            elif disease_condition in sample_name:
                condition = disease_condition_name
                f.write('\t'.join([sample_name, condition]) + '\n')
                column_names.append(sample_name)
            else:
                continue

    # write disease_conditon & wt_condition specific conditions to file
    df_selected = df[column_names]
    df_selected.to_csv(count_data_file,sep="\t",quoting=csv.QUOTE_NONE, index=False)
    
    # run DESeq2
    [res,nMatrix]=runDESeq2(count_data_file, meta_data_file, disease_condition_name, wt_condition_name)

    # get normalized count data for log2fc
    nMatrix = nMatrix.set_axis(nMatrix['Gene']) # set rowname to gene column
    nMatrix = nMatrix.drop(['Gene'], axis=1) # then drop gene column
    # drop all rows with 0s in all columns
    nMatrix_filtered = nMatrix.loc[~(nMatrix==0).all(axis=1)]
    nMatrix_filtered.to_csv(normalized_matrix_file,sep="\t",quoting=csv.QUOTE_NONE)

    # process result matrix
    res = res.dropna() # drop all N/As
    res = res.set_axis(res['Gene'])
    res = res.drop(['Gene'], axis=1) # drop 'Gene' column
    res.to_csv(res_matrix_file,sep="\t",quoting=csv.QUOTE_NONE)

    # get log2fc for all genes
    find_log2fc(res_matrix_file, gene_log2fc_file)

    # further process result matrix to get DE genes
    res = res.loc[(res['baseMean'] > 500)] 
    res = res.loc[(res['padj'] < 0.05)] 
    # get all DE genes
    de = res[(res['log2FoldChange'] > 1.5) | (res['log2FoldChange'] < -1.5)]
    de.to_csv(de_genes_file,sep="\t",quoting=csv.QUOTE_NONE)


# NEED TO FINISH THE FOLLOWING FUNCTION LATER FOR THE ORIGINAL 3 DATASETS (INFLUENZA, COVID, MACROPHAGE)
# def get_pos_and_neg_de_genes(raw_counts, name, disease_condition, wt_condition, disease_condition_name, wt_condition_name):

def main():

    # IPF_LF comparisons
    get_all_de_genes(osp.join('data', 'experimental_raw_counts.csv'), '5_IPF_SiC_Un_vs._IPF_SiC_TG', 'IPF_SiC_Un', 'IPF_SiC_TG', 'Control', 'TGF-B_treated')
    get_all_de_genes(osp.join('data', 'experimental_raw_counts.csv'), '6_IPF_SiC_Un_vs._IPF_SiH_Un', 'IPF_SiC_Un', 'IPF_SiH_Un', 'Control', 'HuR_KO')
    get_all_de_genes(osp.join('data', 'experimental_raw_counts.csv'), '7_IPF_SiH_Un_vs._IPF_SiH_TG', 'IPF_SiH_Un', 'IPF_SiH_TG', 'Control', 'TGF-B_treated')
    get_all_de_genes(osp.join('data', 'experimental_raw_counts.csv'), '8_IPF_SiC_TG_vs._IPF_SiH_TG', 'IPF_SiC_TG', 'IPF_SiH_TG', 'Control', 'HuR_KO')

    # NLF & IPF-LF comparisons
    get_all_de_genes(osp.join('data', 'experimental_raw_counts.csv'), '9_NLF_SiC_Un_vs._IPF_SiC_Un', 'NLF_SiC_Un', 'IPF_SiC_Un', 'Control', 'IPF')
    get_all_de_genes(osp.join('data', 'experimental_raw_counts.csv'), '10_NLF_SiC_TG_vs._IPF_SiC_TG', 'NLF_SiC_TG', 'IPF_SiC_TG', 'Control', 'IPF')
    get_all_de_genes(osp.join('data', 'experimental_raw_counts.csv'), '11_NLF_SiH_Un_vs._IPF_SiH_Un', 'NLF_SiH_Un', 'IPF_SiH_Un', 'Control', 'IPF')
    get_all_de_genes(osp.join('data', 'experimental_raw_counts.csv'), '12_NLF_SiH_TG_vs._IPF_SiH_TG', 'NLF_SiH_TG', 'IPF_SiH_TG', 'Control', 'IPF')


if __name__=='__main__':
    main()