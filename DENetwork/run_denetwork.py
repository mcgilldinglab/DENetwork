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
from simple_tools import check_create_dir, pickle_load

class RunModel:
	def __init__(self, script_dir, name, DEfile, genefile, TFfile, recepfile, PPIfile):
		self.script_dir = script_dir
		self.name = name
		self.DEfile = DEfile
		self.genefile = genefile
		self.TFfile = TFfile # TFfile = '' if using DE genes as the targets
		self.recepfile = recepfile
		self.PPIfile = PPIfile

	def run_model(self):

		# initialize graph

		g = Graph(self.script_dir, self.name, self.DEfile, self.genefile, self.TFfile, self.recepfile, self.PPIfile)

		
		# --------------------STEP 1--------------------
		''' init_graph
			adj = True/False, adjust p-values using FDR
			adj_type = 'FDR'/'BON' (FDR or Bonferroni; multiple testing correction for p-values)
			num_tfs = 'all'/int # of enriched TFs to select 
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
			k = # of edges in path (only paths with k # of edges are kept) between source-target pairs
			r = True/False, whether to run shortest_path, if 'False' will assume that shortest_paths.pickle exists within the 'files' directory
			p_val = p-value cutoff for finding significant shortest paths
			returns a new graph w/ candidate paths of p-value < p_val
		'''
		g.filter_paths(k=5, r=True, p_val=0.02)


		# --------------------STEP 3--------------------
		''' top_N_paths
			N = number of top paths with the maximal path scores to keep
			p = print graph info True/False
		'''
		g.top_N_paths(N=5, p=True)

		
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
			k = # of edges allowed when finding the shortest path between source-target pairs
			l = # of paths to keep for each source-target pair in filter_paths2
			N = # of top paths to keep
			i_local = total number of paths to remove one by one when finding the best local optimum solution
			t = # of trys allowed if no improvement is seen at each continuous iteration
			curve = 'concave'/'convex'; shape of plot in best_signaling_network
			seed = size of random seed used to shuffle list of random edges
			p = fraction of random edges to draw to total list of available edges
		'''
		g.optimal_graph(k=5, l=8, N=5, i_local=50, t=5, curve="concave", seed=0, p=0.001)

		# --------------------STEP 6--------------------
		''' rank_nodes
			ranks the nodes in the (near) optimal solution from most important to least important, and writes them to file
		'''
		g.rank_nodes()


		# --------------------STEP 7--------------------
		''' SignificantGenes
			g = object of (near) global optimum solution
			n = number of significant genes to keep
		'''

		sg = SignificantGenes(self.script_dir, g=g, n=100) 
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
	# name, DEfile, genefile, recepfile
	parser = argparse.ArgumentParser()
	parser.add_argument('-n', '--name', required=True, help='name of the DENetwork model')  
	parser.add_argument('-d', '--de_file', required=True, help='file of differentially-expressed genes')
	parser.add_argument('-g', '--gene_file', required=True, help='file containing all genes, their (output result matrix from DESeq2)')
	parser.add_argument('-r', '--recep_file', required=True, help='file of disease-specific receptors')
	parser.add_argument('-t', '--targets', required=True, help='choose whether the targets are differentially-expressed genes (de) OR transcription factors (tf)', choices=['de', 'tf'])
	args = parser.parse_args()

	# TFfile and PPIfile are the same for all models
	tf_file = '' # DE genes as targets
	if args.targets == 'tf': # TFs as targets
		tffile = osp.join(data_dir, 'TFs.txt')
	ppi_file = osp.join(data_dir, 'ppi_ptm_pd_hgnc.txt')

	# script_dir, name, DEfile, genefile, TFfile, recepfile, PPIfile
	r = RunModel(script_dir, args.name, args.de_file, args.gene_file, tf_file, args.recep_file, ppi_file)
	r.run_model()


if __name__=='__main__':
	main()



