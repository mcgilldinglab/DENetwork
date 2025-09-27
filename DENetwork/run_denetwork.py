'''

This script contains the RunModel class for running DENetwork.
Author: Ting-Yi Su ting-yi.su@mail.mcgill.ca

'''
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
from denetwork import Data, Graph, SignificantGenes
import random
import argparse
from simple_tools import check_create_dir
import tracemalloc

class RunModel:
	def __init__(self, script_dir, name, DEfile, genefile, TFfile, recepfile, PPIfile, num_significant_genes, gamma, max_num_intermediate_nodes_in_shortest_paths, max_p_value_of_significant_paths, num_top_paths, percent_improvement_for_convergence, num_trys):
		self.script_dir = script_dir
		self.name = name
		self.DEfile = DEfile
		self.genefile = genefile
		self.TFfile = TFfile # TFfile = '' if using DE genes as the targets
		self.recepfile = recepfile
		self.PPIfile = PPIfile
		self.num_significant_genes = num_significant_genes
		self.gamma = gamma
		self.K = max_num_intermediate_nodes_in_shortest_paths
		self.p = max_p_value_of_significant_paths
		self.N = num_top_paths
		self.v = percent_improvement_for_convergence
		self.t = num_trys
		# print hyperparameters
		print('DENetwork hyperparameters:')
		print('gamma:', self.gamma)
		print('K:', self.K)
		print('p:', self.p)
		print('N:', self.N)
		print('v:', self.v)
		print('t:', self.t)

	def run_model(self):

		# initialize graph

		g = Graph(self.script_dir, self.name, self.DEfile, self.genefile, self.TFfile, self.recepfile, self.PPIfile, self.gamma)

		
		# --------------------STEP 1--------------------
		''' init_graph
			adj = True/False, adjust p-values using FDR (only for when using TFs as targets)
			adj_type = 'FDR'/'BON' FDR or Bonferroni; multiple testing correction for p-values (only for when using TFs as targets)
			num_tfs = 'all'/int # of enriched TFs to select (only for when using TFs as targets)
			returns the targets of the network

			default parameters: adj=True, adj_type='FDR', num_tfs='all'
		'''
		targets = g.init_graph() 
		
		# no model can be created if have no targets present
		if targets == '':
			print('No targets were found. Terminating program...')
			return


		# --------------------STEP 2--------------------
		''' filter paths
			K = # of edges in path (only paths with k # of edges are kept) between source-target pairs
			r = True/False, whether to run shortest_path, if 'False' will assume that shortest_paths.pickle exists within the 'files' directory
			p_val = p-value cutoff for finding significant shortest paths
			returns a new graph w/ candidate paths of p-value < p_val
		'''
		g.filter_paths(K=self.K, r=True, p_val=self.p)
		# g.filter_paths(k=5, r=True, p_val=0.02)	


		# --------------------STEP 3--------------------
		''' top_N_paths
			N = number of top paths with the maximal path scores to keep
			p = print graph info True/False
		'''
		g.top_N_paths(N=self.N, p=True)
		# g.top_N_paths(N=5, p=True)
		
		# --------------------STEP 4--------------------
		''' best_signaling_network
			i = total number of paths to remove one by one
			j = current number of times function is called (should ALWAYS be -1 here when finding the local optimum); is used to distinguish whether the local or global optimum is being found
			curve = shape of the graph, either 'concave' or 'convex'
		'''
		# get local optimum
		g.best_signaling_network(i=50, j=-1, curve="concave")

		# print info
		g.print_info()

		# NOTE: the local optimal network will not change no matter how many times run_denetwork.py is run
		# only the (near) global optimal network will be different at each run (due to random edges being added)

		# --------------------STEP 5--------------------
		''' optimal_graph
			K = # of edges allowed when finding the shortest path between source-target pairs
			l = # of paths to keep for each source-target pair in filter_paths2
			N = # of top paths to keep
			i_local = total number of paths to remove one by one when finding the best local optimum solution
			t = # of trys allowed if no improvement is seen at each continuous iteration
			curve = 'concave'/'convex'; shape of plot in best_signaling_network
			seed = size of random seed used to shuffle list of random edges
			p = fraction of random edges to draw to total list of available edges
			v = %improvement required to determine convergence of global optimal solution (& to terminate program)
		'''
		g.optimal_graph(K=self.K, l=self.K+3, N=self.N, i_local=50, t=self.t, curve="concave", seed=0, p=0.001, v=self.v)
		# g.optimal_graph(k=5, l=8, N=5, i_local=50, t=5, curve="concave", seed=0, p=0.001, v=5)

		# --------------------STEP 6--------------------
		''' rank_nodes
			ranks the nodes in the (near) optimal solution from most important to least important, and writes them to file
		'''
		g.rank_nodes()


		# --------------------STEP 7--------------------
		''' SignificantGenes
			g = object of (near) global optimual solution
			n = number of significant genes to keep
		'''
		sg = SignificantGenes(self.script_dir, g=g, n=self.num_significant_genes) 
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
	parser.add_argument('-n', '--name', required=True, help='name of the DENetwork model')  
	parser.add_argument('-d', '--de_file', required=True, help='tab-delimited file (.tsv) of differentially-expressed genes, along with their baseMean, log2FoldChange, lfcSE, stat, pvalue, and padj values')
	parser.add_argument('-g', '--gene_file', required=True, help='tab-delimited file (.tsv) containing all genes and their |log2FoldChange|')
	parser.add_argument('-r', '--recep_file', required=True, help='text file (.txt) of disease-specific receptors')
	parser.add_argument('-t', '--targets', required=True, help='choose whether the targets are differentially-expressed genes (de) OR transcription factors (tf)', choices=['de', 'tf'])
	parser.add_argument('-s', '--num_significant_genes', required=True, type=int, help='number of top-ranked genes to select from the (near) optimal network; these selected genes make up the final network (signaling/regulatory pathway) and can be used for GO biological process & reactome pathway analyses')
	
	# optional arguments
	parser.add_argument('-G','--gamma', required=False, default=1.0, type=float, help='Float, Optional, is 1.0 by default. Used to change how much the log2 fold changes (node attributes) influence the path score.')
	parser.add_argument('-K', '--max_num_intermediate_nodes_in_shortest_paths', required=False, default=5, type=int, help='Integer, Optional, is 5 by default. Parameter for finding shortest paths. Is the max number of intermediate nodes allowed when finding shortest paths between source and target pairs.')
	parser.add_argument('-p', '--max_p_value_of_significant_paths', required=False, default=0.01, type=float, help='Float, Optional, is 0.01 by default. Is used as a cut-off for finding significant shortest paths.')
	parser.add_argument('-N', '--num_top_paths', required=False, default=5, type=int, help='Integer, Optional, is 5 by default. Parameter for number of top significant shortest paths kept between each source and target pair.')
	parser.add_argument('-v', '--percent_improvement_for_convergence', required=False, default=5, type=float, help='Float, Optional, is 5.0 by default. Percent improvement required to determine convergence of global optimal solution (& to terminate iterations).')
	parser.add_argument('-T', '--num_trys', required=False, default=5, type=int, help='Int, Optional, is 5 by default. Number of trys allowed for continuous no improvements before terminating & determining the converged solution.')
	args = parser.parse_args()

	# TFfile and PPIfile are the same for all models
	tf_file = '' # DE genes as targets
	if args.targets == 'tf': # TFs as targets
		tffile = osp.join(data_dir, 'TFs.txt')
	ppi_file = osp.join(data_dir, 'ppi_ptm_pd_hgnc.txt')


	# record running time and total memory usage
	tracemalloc.start()
	start = timeit.default_timer()

	# script_dir, name, DEfile, genefile, TFfile, recepfile, PPIfile
	r = RunModel(script_dir, args.name, args.de_file, args.gene_file, tf_file, args.recep_file, ppi_file, args.num_significant_genes, args.gamma, args.max_num_intermediate_nodes_in_shortest_paths, args.max_p_value_of_significant_paths, args.num_top_paths, args.percent_improvement_for_convergence, args.num_trys)
	r.run_model()

	stop = timeit.default_timer()

	# print total running time (in seconds) & max memory usage (in bytes)
	print('TOTAL MEMORY USAGE (BYTES):', tracemalloc.get_traced_memory())
	print('TOTAL RUNNING TIME (SECONDS): ', stop - start)

	tracemalloc.stop()


if __name__=='__main__':
	main()



