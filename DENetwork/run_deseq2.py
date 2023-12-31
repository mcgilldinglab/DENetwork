'''

This script contains the RunDESeq2 class for running DESeq2.
Author: Ting-Yi Su ting-yi.su@mail.mcgill.ca

'''
import numpy as np
import pandas as pd
import csv
import math
import pandas as pd
import os.path as osp
import argparse
from simple_tools import check_create_dir, get_list_from_delimited_file
from diffexpr.py_deseq import py_DESeq2

class RunDESeq2:
	def __init__(self, name, raw_counts, wt_condition, disease_condition, output_dir):
		self.name = name
		self.raw_counts = raw_counts
		self.wt_condition = wt_condition
		self.disease_condition = disease_condition
		self.output_dir = output_dir
		check_create_dir(self.output_dir)

	def runDESeq2(self, data_matrix_file, meta_matrix_file, A, B):
		data_matrix = get_list_from_delimited_file(data_matrix_file, '\t')
		cells = data_matrix[0][1:]
		genes = [item[0] for item in data_matrix[1:]]
		data = [item[1:] for item in data_matrix[1:]]
		pdata = pd.DataFrame(data =np.array(data,dtype=float))
		pdata.columns = cells
		pdata['Gene'] = genes
		
		meta_matrix = get_list_from_delimited_file(meta_matrix_file, '\t')
		pmeta = pd.DataFrame(data=meta_matrix[1:],columns=meta_matrix[0])
		pmeta.index=[item[0] for item in meta_matrix[1:]]
		
		dds = py_DESeq2(count_matrix = pdata,
				   design_matrix = pmeta,
				   design_formula = '~ condition',
				   gene_column = 'Gene') # <- telling DESeq2 this should be the gene ID column
		dds.run_deseq()
		normalized_matrix = dds.normalized_count()
		dds.get_deseq_result(contrast=['condition',A,B])
		res = dds.deseq_result
		return [res, normalized_matrix]

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

	def get_all_de_genes(self, base_mean_threshold, p_value_threshold_type, p_value_threshold, log2_fold_change_threshold):
		print('DESeq2 hyperparameters:')
		print('base mean threshold:', base_mean_threshold)
		print('p-value threshold type:', p_value_threshold_type)
		print('p-value threshold:', p_value_threshold)
		print('log2fc threshold:', '<' + str(-log2_fold_change_threshold), 'OR', '>' + str(log2_fold_change_threshold))
		print('Finding all DE genes for:', self.name)
		
		df = []
		if '.xlsx' in self.raw_counts:
			df = pd.read_excel(self.raw_counts, index_col=None)
		elif '.csv' in self.raw_counts:
			df = pd.read_csv(self.raw_counts, index_col=None) # default delimiter is a comma
		elif '.tsv' in self.raw_counts:
			df = pd.read_csv(self.raw_counts, index_col=None, sep='\t')
		else:
			print('Error! File type of RNA-seq data needs to be either .xlsx (Excel), .csv (comma-delimited), or .tsv (tab-delimited)!')
			return

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
		[res, normalized_matrix] = self.runDESeq2(count_data_file, meta_data_file, self.disease_condition, self.wt_condition)

		# get normalized count data for log2fc
		normalized_matrix = normalized_matrix.set_axis(normalized_matrix['Gene']) # set rowname to gene column
		normalized_matrix = normalized_matrix.drop(['Gene'], axis=1) # then drop gene column
		# drop all rows with 0s in all columns
		nMatrix_filtered = normalized_matrix.loc[~(normalized_matrix==0).all(axis=1)]
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

		print('Number of DE genes:', len(de_all.index))
		print('Number of upregulated DE genes:', len(de_pos.index))
		print('Number of downregulated DE genes:', len(de_neg.index))


def main():

	# reguired arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('-n', '--name', required=True, help='name of disease condition; used to name the output files')  
	parser.add_argument('-r', '--raw_counts', required=True, help='comma- (.csv), tab-delimited (.tsv) or excel (.xlsx) file of raw RNA-seq gene expression data counts for samples of wildtype and disease conditions (samples of irrelevant conditions will be filtered out)')
	parser.add_argument('-w', '--wt_condition', required=True, help='wildtype (healthy) condition')
	parser.add_argument('-d', '--disease_condition', required=True, help='disease condition')
	parser.add_argument('-o', '--output_dir', required=True, help='output directory where output files are stored')

	# optional arguments
	parser.add_argument('-b','--base_mean_threshold',required=False, default=50, type=float, help='Integer, Optional, is 50 by default. Filters the DESeq2 result matrix to remove genes with low expression counts.')
	parser.add_argument('-t','--p_value_threshold_type',required=False, default='padj', help='padj/pvalue, Optional, is padj by default. Choose whether to use unadjusted (pvalue) or FDR adjusted (padj) p-values to filter the DESeq2 result matrix.', choices=['padj', 'pvalue'])
	parser.add_argument('-p','--p_value_threshold',required=False, default=0.05, type=float, help='Float, Optional, is 0.05 by default. Filters the DESeq2 result matrix to remove insignificant genes.')
	parser.add_argument('-l','--log2_fold_change_threshold',required=False, default=0.6, type=float, help='Positive Float, Optional, is 0.6 by default for upregulated DE genes and -0.6 for downregulated DE genes. Filters DE genes.')

	args = parser.parse_args()

	r = RunDESeq2(args.name, args.raw_counts, args.wt_condition, args.disease_condition, args.output_dir)
	r.get_all_de_genes(args.base_mean_threshold, args.p_value_threshold_type, args.p_value_threshold, args.log2_fold_change_threshold)


if __name__=='__main__':
	main()