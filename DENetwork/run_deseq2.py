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

class RunDESeq2:
    def __init__(self, name, raw_counts, wt_condition, disease_condition, output_dir):
        self.name = name
        self.raw_counts = raw_counts
        self.wt_condition = wt_condition
        self.disease_condition = disease_condition
        self.output_dir = output_dir

    def runDESeq2(self, datMatrixFile,metaMatrixFile,A,B):
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

    # get gene_loc2fc here
    def find_log2fc(self, res_matrix_file, gene_log2fc_file):
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

    def get_all_de_genes(self, args.base_mean_threshold, args.p_value_threshold_type, args.p_value_threshold, args.log2_fold_change_threshold):
        print('Finding all DE genes for:', name)
        df = pd.read_csv(raw_counts, index_col=None)

        sample_names = list(df)[1:] # sample names
        column_names = []

        # files to write to & use
        meta_data_file = osp.join(self.output_dir, self.name + '_meta_data.tsv')
        count_data_file = osp.join(self.output_dir, self.name + '_count_data.tsv')
        res_matrix_file = osp.join(self.output_dir, 'res_matrix_' + self.name + '.tsv')
        gene_log2fc_file = osp.join(self.output_dir, 'gene_log2fc_' + self.name + '.tsv')
        de_all_genes_file = osp.join(self.output_dir, 'DE_all_' + self.name + '.tsv')
        de_pos_genes_file = osp.join(self.output_dir, 'DE_pos_' + self.name + '.tsv')
        de_neg_genes_file = osp.join(self.output_dir, 'DE_neg_' + self.name + '.tsv')
        normalized_matrix_file = osp.join(self.output_dir, self.name + '_normalized_matrix.tsv')

        # create meta data file
        with open(meta_data_file, 'w') as f:
            f.write('\t'.join(['sample', 'condition']) + '\n')
            column_names = ['Gene']
            for sample_name in sample_names: # remove gene name column
                condition = ''
                if self.wt_condition in sample_name:
                    f.write('\t'.join([sample_name, self.wt_condition]) + '\n')
                    column_names.append(sample_name)
                elif self.disease_condition in sample_name:
                    f.write('\t'.join([sample_name, self.disease_condition]) + '\n')
                    column_names.append(sample_name)
                else:
                    continue

        # write disease_condition & wt_condition specific conditions to file
        df_selected = df[column_names]
        df_selected.to_csv(count_data_file,sep="\t",quoting=csv.QUOTE_NONE, index=False)
        
        # run DESeq2
        [res,nMatrix] = self.runDESeq2(count_data_file, meta_data_file, self.disease_condition, self.wt_condition)

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
        self.find_log2fc(res_matrix_file, gene_log2fc_file)

        # further process result matrix to get DE genes
        res = res.loc[(res['baseMean'] > base_mean_threshold)] 
        res = res.loc[(res[p_value_threshold_type] < p_value_threshold)] 

        # get all DE genes
        de_all = res[(res['log2FoldChange'] > log2_fold_change_threshold) | (res['log2FoldChange'] < -log2_fold_change_threshold)]

        # get upregulated DE genes
        de_pos = de_all.loc[(de_all['log2FoldChange'] > 0)] 
        # get downregulated DE genes
        de_neg = de_all.loc[(de_all['log2FoldChange'] < 0)] 

        # write de_all, de_pos, and de_neg to file
        de_all.to_csv(de_all_genes_file,sep="\t",quoting=csv.QUOTE_NONE)
        de_pos.to_csv(de_pos_genes_file,sep="\t",quoting=csv.QUOTE_NONE)
        de_neg.to_csv(de_neg_genes_file,sep="\t",quoting=csv.QUOTE_NONE)


def main():

    # reguired arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--name', required=True, help='name of disease condition, used to name the output files')  
    parser.add_argument('-r', '--raw_counts', required=True, help='comma- (.csv) or tab-delimited (.tsv) file of raw RNA-seq gene expression data counts that contain samples of the wildtype and disease conditions (if includes other conditions, they will be filtered out)')
    parser.add_argument('-w', '--wt_condition', required=True, help='wildtype (healthy) condition')
    parser.add_argument('-d', '--disease_condition', required=True, help='disease condition')
    parser.add_argument('-o', '--output_dir', required=True, help='output directory where output files are stored')

    # optional arguments
    parser.add_argument('-b','--base_mean_threshold',required=False, default=50, help='Integer, Optional, is 50 by default. Filters the DESeq2 result matrix to remove genes with low expression counts.')
    parser.add_argument('-t','--p_value_threshold_type',required=False, default='padj', help='padj/pvalue, Optional, is padj by default. Choose whether to use the unadjusted (pvalue) or adjusted (padj) p-values to filter DESeq2 result matrix.', choices=['padj', 'pvalue'])
    parser.add_argument('-p','--p_value_threshold',required=False, default=0.05, help='Float, Optional, is 0.05 by default. Filters the DESeq2 result matrix to remove insignificant genes.')
    parser.add_argument('-l','--log2_fold_change_threshold',required=False, default=0.6, help='Positive Float, Optional, is 0.6 by default for upregulated DE genes and -0.6 for downregulated DE genes. Filters DE genes.')

    args = parser.parse_args()

    r = RunDESeq2(args.name, args.raw_counts, args.wt_condition, args.disease_condition, args.output_dir)
    r.get_all_de_genes(args.base_mean_threshold, args.p_value_threshold_type, args.p_value_threshold, args.log2_fold_change_threshold)


if __name__=='__main__':
    main()