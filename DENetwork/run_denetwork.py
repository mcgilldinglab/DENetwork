import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import statistics as stats
import scipy.stats
import math
from kneed import KneeLocator
import random
import numpy as np
from numpy.linalg import norm
from multiprocessing import Pool 
from functools import partial
import timeit
import math
import warnings
import copy
import pickle
import os
import os.path as osp
from ast import literal_eval
from model import Data, Graph, SignificantGenes
import random
import argparse
from simple_tools import check_create_dir

class RunModel:
	def __init__(self, model_type_list, data_dir, recepfile):
		'''
		Types of models listed as 'target_type'_'DE_gene_type'_'disease_type'
			- types of targets include 'TF' (transcription factors) or 'DE' (differentially expressed genes)
			- types of DE genes include 'pos' (upregulated DE genes), 'neg' (downregulated DE genes), or 'all' (all DE genes)
		'''
		self.model_type_list = model_type_list # lists of types of models ('TF_pos_*', 'TF_neg_*', 'DE_pos_*', 'DE_neg_*', 'TF_all_*', 'DE_all_*'), * = disease
		self.data_dir = data_dir
		self.recepfile = recepfile

	def run_all_models(self): # for experimental comparisons
		# run each model type
		for model_type in self.model_type_list:
			print('RUNNING MODEL:', model_type, '----------------')
			# initialize parameters
			disease_name = model_type[7:]
			DEfile = osp.join(self.data_dir, 'DE' + model_type[2:] + '.tsv')
			genefile = osp.join(self.data_dir, '_'.join(['gene_log2fc', disease_name + '.tsv']))
			TFfile = osp.join(self.data_dir, 'TFs.txt') # same for every disease (b/c are looking for enriched TFs)
			recepfile = ''
			if self.recepfile == '':
				recepfile = osp.join(self.data_dir, '_'.join(['receptors', disease_name + '.txt']))
			else:
				recepfile = self.recepfile
			PPIfile = osp.join(self.data_dir, 'ppi_ptm_pd_hgnc.txt') # also same for every disease
			print('DE file:', DEfile)
			print('gene file:', genefile)
			print('receptor file:', recepfile)
			# run model
			if 'TF_' in model_type: # run with TFs as targets
				self.run_model(model_type, DEfile, genefile, TFfile, recepfile, PPIfile)
			elif 'DE_' in model_type: # run with DE genes as targets, put '' in the TFfile parameter
				self.run_model(model_type, DEfile, genefile, '', recepfile, PPIfile)

	def run_model(self, model_type, DEfile, genefile, TFfile, recepfile, PPIfile):

		# initialize graph
		g = Graph(model_type, DEfile, genefile, TFfile, recepfile, PPIfile)
		targets = g.init_graph(True, 'FDR') # default is adjust p-values w/ FDR
		
		# no model can be created if have no targets present
		if targets == '':
			print('No targets were found. Terminating program...')
			return

		# filter paths
		g.filter_paths(5, True, 0.02)

		# get top N paths
		g.top_N_paths(5, True)

		# get local optimum
		g.best_signaling_network(50, -1, "concave")

		# print info
		g.print_info()

		# get global optimum
		# k, l, N, i_local, t, curve, direction, seed
		g.optimal_graph(5, 8, 5, 50, 5, "concave", 0, 0.001) #default seed 0

		# rank nodes
		g.rank_nodes()

		# get top 100 significant genes
		sg = SignificantGenes(g, 100) # get 100 top siginificant genes
		sg.get_ranking()
		sg.get_go()
		sg.get_noa_sif()


def main():

	script_dir = osp.dirname(__file__)
	data_dir = osp.join(script_dir, 'data')
	files_dir = osp.join(script_dir, 'files')
	figures_dir = osp.join(script_dir, 'figures')
	results_dir = osp.join(script_dir, 'results')

	# create data, files, figures, and results directories if they don't exist already
	check_create_dir(data_dir)
	check_create_dir(files_dir)
	check_create_dir(figures_dir)
	check_create_dir(results_dir)

	# input arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--data_dir')  #/home/tingyisu/scratch/modeller/files
	parser.add_argument('-p', '--pdb_download_dir') # /home/tingyisu/scratch/pdb_cif
	args = parser.parse_args()



