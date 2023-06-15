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
from statsmodels.stats.multitest import fdrcorrection
from simple_tools import check_create_dir, pickle_load, pickle_dump, write_list_to_file_convert_to_str

class Data:
	def __init__(self, DEfile, genefile, TFfile, recepfile, PPIfile):
		self.DEfile = DEfile
		self.genefile = genefile
		self.TFfile = TFfile # TFfile = '' if using DE genes as the targets
		self.recepfile = recepfile
		self.PPIfile = PPIfile
		self.genes_log2fc = {}
		self.ppis = {} # dict (key = tuple of interacting pair, value = log2 ppi score)
		self.ppis_dict = {} # dict of dict of interacting pairs, value = log2 ppi score

	def get_TFs(self, DEgenes, num_genes, adj, adj_type, num_tfs):
		'''
		TFfile = file of TFs and their target genes
		DEgenes = list of DE genes
		num_genes = total number of genes
		adj = T/F, adjust p-values
		adj_type = 'FDR'/'BON' (FDR or Bonferroni correction)
		num_tfs = 'all'/int # of enriched TFs to select 
		returns a list of TFs w/ p-value < 0.05
		'''
		# create a dictionary of TFs; key = TF, value = list of gene targets
		tfs_targets = {}
		with open(self.TFfile, 'r') as tf_f:
			next(tf_f)
			for line in tf_f: 
				items = line.strip().split('\t')
				if items[0] not in tfs_targets:
					tfs_targets[items[0].upper()] = [items[1].upper()]
				else: 
					tfs_targets[items[0].upper()].append(items[1].upper())


		# get p-values for each TF
		pbg = len(DEgenes)/num_genes # background probability of DE genes
		print('Background probability:', pbg)

		tf_name_list = list(tfs_targets.keys())
		tf_list = []

		# store all p-vals
		p_val_list = []
		
		for tf in tfs_targets:
			n = len(tfs_targets[tf]) # total number of target genes for this TF
			x = len([target for target in tfs_targets[tf] if target in DEgenes]) # number of target genes that are differential for this TF
			p_val = 1 - scipy.stats.binom.cdf(x-1, n, pbg)
			p_val_list.append(p_val)

		final_pval_list = p_val_list.copy()
		p_val_cutoff = 0.05
		# change final_pval_list or p_val_cutoff based on adj_type
		if adj and adj_type == 'FDR':
			rejected, final_pval_list = fdrcorrection(p_val_list)
			print('Finding enriched TFs using FDR correction')
		if adj and adj_type == 'BON':
			p_val_cutoff = 0.05/len(p_val_list)
			print('Finding enriched TFs using Bonferroni correction')

		if num_tfs == 'all':
			# get tfs with pvals < 0.05
			for i in range(len(tf_name_list)):
				if final_pval_list[i] < p_val_cutoff:
					tf_list.append(tf_name_list[i])
			print(str(len(tf_list)) + '/' + str(len(tfs_targets)) + ' TFs are enriched.')
		else:
			tf_pval_dict = {} # key = TF, value = final pval
			for i in range(len(tf_name_list)):
				tf_pval_dict[tf_name_list[i]] = final_pval_list[i]
			# sort 
			sorted_tf_pval_dict = {k: v for k, v in sorted(tf_pval_dict.items(), key=lambda item: item[1])} # sort in ascending order
			# get top tfs
			tf_list = list(sorted_tf_pval_dict.keys())[:num_tfs]
			print(str(len(tf_list)) + '/' + str(len(tfs_targets)) + ' enriched TFs were selected.')

		return tf_list
	
	# get log2fc from file
	def get_log2fc(self):
		with open(self.genefile, 'r') as gene_f:
			for line in gene_f:
				items = line.strip().split('\t')
				self.genes_log2fc[items[0].upper()] = float(items[1])
	
	# creates a PPI dict of dict; for adding edge PPI scores to new graphs
	def get_ppi_dict(self): 
		with open(self.PPIfile, 'r') as ppis_f:
			# get ppis
			ppis_all = {} 
			# add all ppi scores to a dict; key = interaction pair (tuple), value = list of ppi scores
			for line in ppis_f:
				items = line.split('\t')
				node1, node2 = items[0].upper(), items[2].upper()
				ppi_score = float(items[3].strip())
				# make sure interaction pair exists and isn't a self-interaction
				if ppi_score != 0.0 and node1 in self.genes_log2fc and node2 in self.genes_log2fc and node1 != node2:
					# add all scores, either pp or pd (any order)
					if (node1, node2) in ppis_all: 
						ppis_all[(node1, node2)].append(ppi_score)
					elif (node2, node1) in ppis_all: 
						ppis_all[(node2, node1)].append(ppi_score)
					else: 
						ppis_all[(node1, node2)] = [ppi_score]

		self.ppis = {k: math.log2(stats.mean(v)) for k, v in ppis_all.items()} # store the mean of all PPIs between pp and pd

		# create a dict of dict of normalized log2 ppi scores
		for key in self.ppis:
			p1, p2 = key
			if p1 not in self.ppis_dict: self.ppis_dict[p1] = {}
			self.ppis_dict[p1][p2] = self.ppis[key]


class Graph:
	def __init__(self, script_dir, name, DEfile, genefile, TFfile, recepfile, PPIfile):
		'''
		script_dir = directory where scripts are located
		name = name of graph
		DEfile = file of DE genes w/ padj < 0.01
		genefile = file of total genes
		wt_indicies = list of column indices of WT sample counts
		ko_indicies = list of column indices of KO sample counts
		TFfile = file of TFs
		recepfile = file of receptors
		PPIfile = file of interacting proteins and their corresponding PPI scores
		'''	
		self.name = name
		self.G = nx.Graph()
		self.nodes = []
		self.edges = []
		self.paths = {}
		self.paths_list = []
		self.num_paths = 0
		self.targets = []
		self.num_targets = 0
		self.recs = []
		self.num_recs = 0
		self.d = Data(DEfile, genefile, TFfile, recepfile, PPIfile) # data
		self.pos_g = nx.Graph() # reference graph for drawing
		self.local_lambda_pen = 0.0 # found in STEP4, use to compare in STEP5
		self.global_lambda_pen = 0.0  # found in STEP5, use to rank nodes
		self.local_SQ = 0.0
		self.global_SQ = 0.0
		self.improv_list = {}
		self.improv_iterations = {}
		# the following are attributes for finding the final network 
		self.l = -1
		self.N = -1
		self.t = -1
		# create self.name specific directories within files and figures folders if they don't exist already
		self.files_dir = osp.join(script_dir, 'files', self.name)
		check_create_dir(self.files_dir)
		self.figures_dir = osp.join(script_dir, 'figures', self.name)
		check_create_dir(self.figures_dir)
		self.results_dir = osp.join(script_dir, 'results', self.name)
		check_create_dir(self.results_dir)
		self.path_scores_fig_dir = osp.join(self.figures_dir, 'path_scores_distribution')
		check_create_dir(self.path_scores_fig_dir)
		self.sq_fig_dir = osp.join(self.figures_dir, 'S_Q')
		check_create_dir(self.sq_fig_dir)
		
		# if not osp.exists('files/' + self.name):
		# 	os.makedirs('files/' + self.name)
		# if not osp.exists('figures/' + self.name + '/path_scores_distribution'):
		# 	os.makedirs('figures/' + self.name + '/path_scores_distribution')
		# if not osp.exists('figures/' + self.name + '/S_Q'):
		# 	os.makedirs('figures/' + self.name + '/S_Q')

	
	def update_all(self, prev_g, prev_paths, prev_paths_list):
		# update node, edge, and path attributes
		self.G = prev_g
		self.nodes = prev_g.nodes
		self.edges = prev_g.edges
		self.num_nodes = prev_g.number_of_nodes()
		self.num_edges = prev_g.number_of_edges() * 2
		self.paths = prev_paths
		self.paths_list = prev_paths_list 
		self.num_paths = len(prev_paths_list)
	
	#------------------------------------------print graph---------------------------------
	
	def print_graph(self):
		for node in self.nodes:
			neighbors = self.G.neighbors(node)
			edges = [str(neighbor) + '(' + str(self.nodes[neighbor]['log2fc']) + '): ' + self.edges[node, neighbor]['ppi'] for neighbor in neighbors]
			print(node + '(' + str(self.nodes[node]['log2fc']) + '):', edges)
	
	# prints #nodes, edges, recs, targets
	def print_info(self):
		print('Graph Details:')
		print('# of nodes =', self.num_nodes)
		print('# of edges =', self.num_edges) # undirected
		print("# of receptors (sources):", self.num_recs)
		if self.d.TFfile == '':
			print('# of DE genes (targets):', self.num_targets)
		else:
			print('# of TFs (targets):', self.num_targets)
		print("# of paths:", self.num_paths)
	
	def draw(self):
		nx.draw(self.G, pos=nx.circular_layout(self.pos_g))

	# for step 5: finding the global optimal graph
	def plot_improvement(self):
		print(self.improv_list, self.improv_iterations)
		plt.plot(self.improv_iterations, self.improv_list)
		plt.xlabel("Iteration number")
		plt.ylabel("% improvement")
		fig_name = osp.join(self.figures_dir, 'step5_improvement_', '_'.join([self.l, str(self.N), str(self.t) + '.svg']))
		plt.savefig(fig_name)

		# plt.savefig('figures/' + self.name + '/step5_improvement_' + self.l + '_' + str(self.N) + '_' + str(self.t) + '.jpg')
	
	# -------------------------------------------helper functions-----------------------------------------
	
	def list_to_dict(self, paths_list):
		'''
		input: a nested list of paths
			paths_list: [[path1], [path2], [path3]]
		output: a nested dictionary of paths with sources and targets as keys
			path_dict: {source1: {target1: [[path1], [path2]], target2: ...}, source2: ...}
		'''
		path_dict = {}
		for path in paths_list:
			if path == []: continue
			if path[0] not in path_dict:
				path_dict[path[0]] = {}
			if path[len(path)-1] not in path_dict[path[0]]:
				path_dict[path[0]][path[len(path)-1]] = [path] # add first
			else:
				temp = path_dict[path[0]][path[len(path)-1]].copy()
				temp.append(path)
				path_dict[path[0]][path[len(path)-1]] = temp
		return path_dict
	
	def get_pairs(self):
		'''
		returns a list of (receptor, target) pairs
		'''
		pairs = []
		for receptor in self.recs:
			for target in self.targets:
				pairs.append((receptor, target))
		return pairs

	def shortest_paths(self, k, rec, target):
		try:
			return [p for p in nx.all_shortest_paths(self.G, source=rec, target=target) if len(p) <= k + 1]
		except nx.exception.NetworkXNoPath:
			return []

	def get_shortest_paths(self, k):
		'''
		k = # of intermediate nodes in path
		returns dict of dict of shortest paths
		'''
		start = timeit.default_timer()
		receptor_target_pairs = self.get_pairs()
		p = Pool()
		func = partial(self.shortest_paths, k)
		paths = p.starmap(func, receptor_target_pairs)
		self.paths_list = sum(paths, [])
		p.close()
		p.join()
		stop = timeit.default_timer()
		print('Time needed to find the shortest paths: ', stop - start)
		return self.list_to_dict(self.paths_list)
	
	def calculate_path_score(self):
		'''
		returns a dict of dict of path scores and a list of scores
		'''
		scores = {}
		scores_list = []
		for receptor in self.paths:
			scores[receptor] = {}
			for target in self.paths[receptor]:
				curr_paths = self.paths[receptor][target] # more than one shortest path
				# score = product of fc of all nodes x product of ppi score of all edges
				# since is in log2, can add them up
				scores[receptor][target] = []
				for path in curr_paths:
					score = 0.0
					for i in range(len(path)-1):
						score += float(self.nodes[path[i]]['log2fc']) + float(self.edges[path[i], path[i+1]]['ppi']) 
					score +=  float(self.nodes[path[len(path) - 1]]['log2fc'])
					scores[receptor][target].append(score)
					scores_list.append(score)
			if not scores[receptor]: # empty
				del scores[receptor]
		return scores, scores_list
	
	def create_graph(self, candidate_paths, remaining_paths): # creates a new graph from the paths found
		'''
		candidate_paths = dict of dict of paths
		remaining_paths = list of paths to add to the new graph
		updates attributes
		'''
		# create new graph from remaining_paths
		self.G = nx.Graph()
		self.paths = candidate_paths
		self.paths_list = remaining_paths
		self.num_paths = len(remaining_paths)
		for path in remaining_paths:
			nx.add_path(self.G, path)	
		# add in node log2fc and edge ppi score
		for node in self.G.nodes:
			self.G.add_node(node, log2fc=self.d.genes_log2fc[node])
		for edge in self.G.edges: # doesn't matter in this case
			node1, node2 = edge
			if node1 in self.d.ppis_dict and node2 in self.d.ppis_dict[node1]: self.G.add_edge(node1, node2, ppi=self.d.ppis_dict[node1][node2])
			elif node2 in self.d.ppis_dict and node1 in self.d.ppis_dict[node2]: self.G.add_edge(node2, node1, ppi=self.d.ppis_dict[node2][node1])   
		# update attributes 
		self.nodes = self.G.nodes
		self.edges = self.G.edges
		self.num_nodes = self.G.number_of_nodes()
		self.num_edges = self.G.number_of_edges() * 2
		self.recs = list(self.paths.keys())
		remaining_targets = [list(v.keys()) for k, v in self.paths.items()]
		self.targets = set(sum(remaining_targets, []))
		self.num_targets = len(self.targets)
		self.recs = [rec for rec in self.recs if rec in self.G.nodes]
		self.num_recs = len(self.recs)
	
	def temp_graph(self, remaining_paths):
		temp_G = nx.Graph()
		for path in remaining_paths:
			nx.add_path(temp_G, path)	
		# add in node log2fc and edge ppi score
		for node in temp_G.nodes:
			temp_G.add_node(node, log2fc=self.d.genes_log2fc[node])
		for edge in temp_G.edges: # doesn't matter in this case
			node1, node2 = edge
			if node1 in self.d.ppis_dict and node2 in self.d.ppis_dict[node1]: temp_G.add_edge(node1, node2, ppi=self.d.ppis_dict[node1][node2])
			elif node2 in self.d.ppis_dict and node1 in self.d.ppis_dict[node2]: temp_G.add_edge(node2, node1, ppi=self.d.ppis_dict[node2][node1])   
		return temp_G
	
	# ----------------------------------------------------MAIN STEPS IN METHOD----------------------------------------------------
	
	# STEP1: initializes a NetworkX undirected graph
	def init_graph(self, adj=True, adj_type='FDR', num_tfs='all'): 
		'''
		adj = T/F, adjust p-values using FDR
		adj_type = 'FDR'/'BON' (FDR or Bonferroni; multiple testing correction for p-values)
		num_tfs = 'all'/int # of enriched TFs to select 
		'''
		
		# get DE genes
		DEgenes = []
		with open(self.d.DEfile, 'r') as DE_f:
			next(DE_f)
			DEgenes = [line.strip().split('\t')[0].upper() for line in DE_f.readlines()]

		print('Number of DE genes:', len(DEgenes))
		# print(DEgenes)

		# get all genes
		self.d.get_log2fc()

		# get targets
		# either TFs or DE genes - positive (upregulated) OR negative (downregulated) OR both
		targets = []
		if self.d.TFfile != '': 
			# get enriched TFs
			targets = self.d.get_TFs(DEgenes, len(self.d.genes_log2fc), adj, adj_type, num_tfs)
			if targets == []:
				return '' # no targets, so can terminate model
		else: # no TF file, then targets are DE genes
			targets = DEgenes

		# get receptors
		recs = []
		with open(self.d.recepfile, 'r') as rec_f:
			rec = [r.strip().upper() for r in rec_f.readlines()]
			recs = list(filter(None, rec))

		# get ppis
		self.d.get_ppi_dict()

		print('-------------------------Initializing graph------------------------')
	
		# add edges
		# self loops (when a gene/protein interacts with itself) are excluded
		for node1 in self.d.ppis_dict:
			for node2 in self.d.ppis_dict[node1]:
				self.G.add_edge(node1, node2, ppi=self.d.ppis_dict[node1][node2]) # add ppi score attribute to edge
				self.G.add_node(node1, log2fc=self.d.genes_log2fc[node1]) # add log2fc attribute to nodes
				self.G.add_node(node2, log2fc=self.d.genes_log2fc[node2])
		
		# update attributes 
		self.nodes = self.G.nodes
		self.edges = self.G.edges
		self.num_nodes = self.G.number_of_nodes()
		self.num_edges = self.G.number_of_edges() * 2
		
		# find number of receptors and targets in the graph
		recs_in = list(set(self.nodes) & set(recs))
		targets_in = list(set(self.nodes) & set(targets))
		
		# update
		self.recs = recs_in
		self.num_recs = len(recs_in)
		self.targets = targets_in
		self.num_targets = len(targets_in)
		
		self.print_info()
		self.pos_g = self.G	

		print('Overlap between sources & targets:', len(set(self.recs) & set(self.targets)))
		
	# STEP2: filter out edges
	def filter_paths(self, k, r, p_val): 
		'''
		k = # of edges in path
		r = True/False, whether to run shortest_path, if 'False' will assume that f is name of file with shortest paths
		p_val = p value cutoff
		returns a new graph w/ candidate paths of p-value < p_val
		'''

		shortest_paths_pickle_file = osp.join(self.files_dir, 'shortest_paths.pickle')
		# fname = 'files/' + self.name + '/shortest_path.pickle'

		if r: 
			self.paths = self.get_shortest_paths(k)
			pickle_dump(self.paths, shortest_paths_pickle_file) # write shortest paths to file

			# with open(fname, 'wb') as f:
			# 	pickle.dump(self.paths, f)
		else:
			self.paths = pickle_load(shortest_paths_pickle_file) # read in shortest paths

			# with open(fname, 'rb') as f:
			# 	self.paths = pickle.load(f)

		print('-------------------------Filtering paths------------------------')

		# calculate score distribution
		scores, scores_list = self.calculate_path_score()
		print("# of paths found with", k, "intermediate nodes =", len(scores_list))
		print("Make sure that path scores are normally distributed:")

		# plot path score distribution (should be roughly Gaussian) & save it to file
		plt.figure()
		plt.hist(scores_list, bins=20)
		path_scores_figname = osp.join(self.path_scores_fig_dir, 'path_scores_distribution.svg')
		plt.savefig(path_scores_figname)
		plt.close()

		# plt.savefig('figures/' + self.name + '/path_scores_distribution/path_scores_distribution.png')

		# calculate p-values
		mu = stats.mean(scores_list)
		sigma = stats.stdev(scores_list)
		deleted = 0
		remaining_paths = []
		candidate_paths = {}
		for receptor in scores:
			candidate_paths[receptor] = {}
			for target in scores[receptor]:
				curr_list = []
				curr_scores = scores[receptor][target]
				for i in range(len(curr_scores)): 
					pval = 1 - 0.5 * (1 + math.erf((curr_scores[i] - mu) / (sigma * math.sqrt(2))))
					if pval < p_val: # get only significant p-values
						curr_list.append(self.paths[receptor][target][i])
						remaining_paths.append(self.paths[receptor][target][i])
					else: # delete insignificant paths
						deleted += 1
				if curr_list != []:
					candidate_paths[receptor][target] = curr_list
		print(deleted, 'paths were removed;', len(scores_list) - deleted, ' paths remain.')
		candidate_paths = {k: v for k, v in candidate_paths.items() if v}

		# create new graph after filtering out insignificant paths
		self.create_graph(candidate_paths, remaining_paths)   
		# print graph info
		self.print_info()

	
	# STEP2: filter out edges
	# different from filter_paths; only use for finding the (near) global optimal graph
	def filter_paths2(self, k, j, l): 
		'''
		k = # of intermediate nodes in path
		j = # of times filter_paths has been called currently (starts from 0)
		l = min # of paths to keep for each receptor-target pair
		returns a new graph w/ candidate paths of p-value < p_val
		'''
		paths = self.get_shortest_paths(k)
		self.paths = paths
		
		# calculate score distribution
		scores, scores_list = self.calculate_path_score()
		print("# of paths found with", k, "intermediate nodes =", len(scores_list))
		print("Make sure that path scores are normally distributed:")

		# plot path score distribution (should be roughly Gaussian) & save it to file
		plt.figure()
		plt.hist(scores_list, bins=20)
		path_scores_figname = osp.join(self.path_scores_fig_dir, 'path_scores_distribution' + str(j) +'.svg')
		plt.savefig(path_scores_figname)
		# plt.savefig('figures/' +  self.name + '/path_scores_distribution/' + 'path_scores_distribution' + str(j) +'.png')
		plt.close()


		# calculate p-values
		mu = stats.mean(scores_list)
		sigma = stats.stdev(scores_list)
		deleted = 0
		remaining_paths = []
		candidate_paths = {}
		for receptor in scores:
			candidate_paths[receptor] = {}
			for target in scores[receptor]:
				curr_list = paths[receptor][target]
				curr_scores = scores[receptor][target]
				curr_pvals = []
				for i in range(len(curr_scores)): 
					pval = 1 - 0.5 * (1 + math.erf((curr_scores[i] - mu) / (sigma * math.sqrt(2))))
					curr_pvals.append(pval)
				# sort curr_list
				curr_pvals.sort(reverse=False)
				if len(curr_pvals) < l:
					remaining_paths.append(curr_list)
					candidate_paths[receptor][target] = paths[receptor][target]
				else:
					keep = sorted(range(len(curr_pvals)), key=lambda x: curr_pvals[x])[:l] # want the minimal p-vals
					to_keep = [curr_list[i] for i in keep]
					remaining_paths.append(to_keep)
					deleted += len(curr_pvals) - len(keep)
					candidate_paths[receptor][target] = to_keep
		print(deleted, 'paths were removed;', len(scores_list) - deleted, ' paths remain.')
		remaining_paths = sum(remaining_paths, [])
		candidate_paths = {k: v for k, v in candidate_paths.items() if v}

		# create new graph after filtering out insignificant paths
		self.create_graph(candidate_paths, remaining_paths)   
		# print graph info & make sure no recs and targets are lost
		self.print_info() 
	
	# STEP3: find the top N paths with the maximal path scores
	def top_N_paths(self, N, p):
		'''
		N = number of top paths with the maximal path scores to keep
		p = print graph info True/False
		'''
		original = self.num_paths
		# calculate path score
		top_paths = {}
		remaining_paths = []
		for receptor in self.paths:
			top_paths[receptor] = {}
			for target in self.paths[receptor]:
				curr_paths = self.paths[receptor][target] # more than one shortest path
				if len(curr_paths) <= N: 
					top_paths[receptor][target] = curr_paths
					remaining_paths.append(curr_paths)
				else:
					scores = []
					for path in curr_paths:
						score = 0.0
						for i in range(len(path)-1):
							score -= float(self.edges[path[i], path[i+1]]['ppi']) + float(self.nodes[path[i]]['log2fc']) + float(self.nodes[path[i+1]]['log2fc'])
						scores.append(score)
					top_N = sorted(range(len(scores)), key=lambda x: scores[x])[:N] # want the minimal path scores
					top = [curr_paths[i] for i in top_N]
					top_paths[receptor][target] = top
					remaining_paths.append(top)
		remaining_paths = sum(remaining_paths, [])
		# create new graph for top N paths
		self.create_graph(top_paths, remaining_paths)
		remaining = self.num_paths
		if p:
			print('-------------------------Finding the top '+ str(N) + ' paths------------------------')
			print(original - remaining, 'paths were removed;', remaining, ' paths remain.')
			self.print_info()

		
	# STEP4: build overall best signaling network (local optimal solution)
	def best_signaling_network(self, i, j, curve):
		'''
		i = total number of paths to remove one by one
		j = current number of times function is called
		curve = shape of the graph, either 'concave' or 'convex'
		'''

		test_g = self.G.copy()
		test_paths = self.paths.copy()
		test_remaining_paths = self.paths_list.copy()

		# add first S(Q)
		_, scores = self.calculate_path_score()
		scores_list = scores.copy()
		total_score = [sum(scores_list)]
		penalty = [sum(d for _, d in test_g.degree)/test_g.number_of_nodes()]

		# remove one path at a time
		for k in range(i):

			# find index of path w/ lowest path score to remove
			min_index = scores_list.index(min(scores_list)) # if more than one min score, will pick first one
			del scores_list[min_index]
			del test_remaining_paths[min_index]
			test_g = self.temp_graph(test_remaining_paths)
			total_score.append(sum(scores_list))
			penalty.append((sum(d for _, d in test_g.degree)/test_g.number_of_nodes()))
		
		# calculate lambda needed to make the score and penalty close to each other
		lambda_pen = stats.mean(total_score)/stats.mean(penalty)
		print('STEP4 lambda penalty: ', lambda_pen)
		actual_penalty = [lambda_pen * i for i in penalty]
		S_Q = [total_score[i] - actual_penalty[i] for i in range(len(total_score))]

		# plot score & penalty separately
		plt.figure()
		plt.plot(total_score)
		plt.plot(actual_penalty)
		plt.xlabel("Number of paths removed")
		plt.legend(["Score", "Penalty"])
		if j == -1:
			plt.savefig(osp.join(self.sq_fig_dir, 'score_penalty.svg'))
			# plt.savefig('figures/' + self.name + '/S_Q/score_penalty.png')
		else: 
			plt.savefig(osp.join(self.sq_fig_dir, 'score_penalty' + str(j+1) + '.svg'))
			# plt.savefig('figures/' + self.name + '/S_Q/score_penalty' + str(j+1) + '.png')
		plt.close()

		# plot S(Q) with line x (which connects the first & last points)
		# find the "elbow"/"knee" using the kneed package
		plt.figure()
		plt.plot(S_Q)
		# plot line from start to end 
		x1 = [0, i]
		y1 = [S_Q[0], S_Q[len(S_Q)-1]]
		plt.plot(x1, y1)
		plt.xlabel("Number of paths removed")
		plt.legend(['S(Q)', 'line X'])
		# x = number of paths removed
		x = [k for k in range(i+1)]

		# find the best "elbow"/"knee" by iterating over different sensitivites
		print('looking for the best "elbow":')
		best_elbow = 0
		best_d = -math.inf
		p1 = np.array([0, S_Q[0]])
		p2 = np.array([i, S_Q[len(S_Q)-1]])
		sensitivity = [1, 3, 5, 10, 100, 200, 400]
		direction = ""
		# get direction ('increasing'/'decreasing') of the graph
		if S_Q[0] > S_Q[len(S_Q)-1]:
			direction = 'decreasing'
		else:
			direction = 'increasing'
		for s in sensitivity:
			with warnings.catch_warnings():
				warnings.filterwarnings('error')
				try:
					kneedle = KneeLocator(x, S_Q, curve=curve, direction=direction, S=s)
					if kneedle.knee_y is not None:
						num = kneedle.elbow
						p3 = np.array([num, kneedle.knee_y])
						# find the best elbow/knee that hars the largest Euclidean distance to line X
						d = norm(np.cross(p2-p1, p1-p3))/norm(p2-p1)
						if d > best_d: 
							best_elbow = num
							best_d = d
						print('elbow:', num, kneedle.knee_y, '; sensitivity:', s) # num paths to remove for best signaling network
					else: print('no elbow found for sensitivity:', s)
				except Warning as e:
					print('no elbow found for sensitivity:', s)
					continue

		print('best elbow:', best_elbow)
		plt.axvline(x=best_elbow, c = 'r')
		# save figure
		if j == -1:
			plt.savefig(self.sq_fig_dir, 'S_Q.svg')
			# plt.savefig('figures/' + self.name + '/S_Q/S_Q.png')
		else: 
			plt.savefig(self.sq_fig_dir, 'S_Q' + str(j+1) + '.svg')
			# plt.savefig('figures/' + self.name + '/S_Q/S_Q' + str(j+1) + '.png')
		plt.close()

		# remove the paths with the lowest scores (# of paths removed = best_elbow)
		to_remove = sorted(range(len(scores)), key=lambda x: scores[x])[:best_elbow] # indices of paths to remove
		local_op_paths = [self.paths_list[i] for i in range(self.num_paths) if i not in to_remove]
		self.create_graph(self.list_to_dict(local_op_paths), local_op_paths)

		# only save initial lambda_pen from first run of STEP4
		if j == -1: 
			self.local_SQ = S_Q[best_elbow]
			self.local_lambda_pen = lambda_pen
			pickle_dump(self, osp.join(self.files_dir, 'local_op_g.pickle'))
			# with open('files/' + self.name + '/local_op_g.pickle', 'wb') as fname:
			# 	pickle.dump(self, fname)

		if j != -1: return lambda_pen

	# STEP5: find (near) global optimal solution
	def optimal_graph(self, k, l, N, i_local, t, curve, seed, p):
		'''
		k = # of edges in shortest path
		l = # of paths to keep for each receptor-target pair in filter_paths2
		N = # of top paths to keep
		i_local = # of iterations used to find the best local signaling network
		t = # of trys if no improvement
		curve = concave/convex; shape of plot in best_signaling_network
		seed = size of random seed used to shuffle list of random edges
		p = fraction of random edges to draw
		'''

		# time amount of time it takes
		start = timeit.default_timer()
		print('local S(Q):',  self.local_SQ)
		print('local lambda_pen:', self.local_lambda_pen)
		_, scores = self.calculate_path_score()
		best_SQ = self.local_SQ
		best_scores = scores
		improv_list = []
		improv_iterations = []
		num_improv = 0
		lambda_pen = 0.0  # lambda_pen in STEP4 (will be updated in each iteration)
		self.l = l
		self.N = N
		self.t = t

		# randomly add edges to local optimal graph
		ppis_set = set(self.d.ppis.keys())

		best_edges = []
		# need to get the complement (e.g. ('B', 'A') for ('A', 'B'))
		# can't get complement when iterating through self.edges, manually add into best_edges
		for edge in list(self.edges):
			node1, node2 = edge
			best_edges.append(edge)
			if (node2, node1) in list(self.edges): print('already exists!') # check
			best_edges.append((node2, node1))

		# get edges from PPI list that aren't in best_g already
		best_set = set(best_edges)
		edges = ppis_set.difference(best_set)
		M = int(len(edges) * p) # randomly draw p number of the edges each time
		print('Number of edges to draw each time:', M)

		# create list to randomly draw edges
		to_draw = []
		num_0 = 0
		num_not_0 = 0
		for edge in edges:
			num = 0
			if (self.d.genes_log2fc[node1] != 0) and (self.d.genes_log2fc[node2] != 0):
				num = int((2**(self.d.ppis[edge] + self.d.genes_log2fc[node1] + self.d.genes_log2fc[node2])) * 100) # num times this edge appears
			curr = [edge] * num
			to_draw += curr # add to list of edges to draw from

		if len(to_draw) == 0:
			print('All the edges with the highest scores exist in the local optimal already!')
			print('The local optimal network is the global optimal. Global lambda penalty = local lambda penalty.')
			return

		# use a seed to shuffle the same way each time
		random.Random(seed).shuffle(to_draw)

		print('Number of unique edges to draw from:', len(set(to_draw)))

		total_num_no_improv = 0
		num_no_improv = 0

		# get recs and targets in graph
		# recs and targets need to stay constant
		recs = self.recs 
		targets = self.targets 

		i = int(len(edges)/M) # max number of iterations
		
		# graph iterative refinement
		# iterate until converges
		for j in range(i): 
			print('---------------------------------ITERATION', j+1, '---------------------------------' )

			# save previous
			prev_g = copy.deepcopy(self.G)
			prev_paths = copy.deepcopy(self.paths)
			prev_paths_list = copy.deepcopy(self.paths_list)
			prev_SQ = best_SQ
			prev_scores = best_scores

			ran_edges = to_draw[j*M:(j+1)*M]
			print('Number of repeated random edges:', len(ran_edges) - len(set(ran_edges)))

			print("prev g:", prev_g.number_of_nodes(), prev_g.number_of_edges() * 2)
			
			# add random edges
			for edge in set(ran_edges):
				node1, node2 = edge
				self.G.add_edge(node1, node2, ppi=self.d.ppis[edge]) # add ppi score attribute to edge
				self.G.add_node(node1, log2fc=self.d.genes_log2fc[node1]) # add log2fc attribute to nodes
				self.G.add_node(node2, log2fc=self.d.genes_log2fc[node2])

			# update nodes and edge attributes
			print("added:", self.G.number_of_nodes() ,"nodes,", self.G.number_of_edges() * 2, " edges", '(' + str((j+1)*M - j*M) + ')')
			
			# redo STEP2 (get new paths) then do STEP3
			self.filter_paths2(k, j, l) # keeps l number of paths for each (receptor, target) pair

			# redo STEP3 
			self.top_N_paths(N, False)

			# redo STEP4
			lambda_pen = self.best_signaling_network(i_local, j, curve)


			# check whether S(Q) has improved
			_, G1_scores = self.calculate_path_score()
			G1_SQ = (sum(G1_scores) - self.local_lambda_pen * (sum(d for _, d in self.G.degree)/self.G.number_of_nodes()))
			
			print('local lambda_pen:', self.local_lambda_pen)
			print("prev g:", prev_g.number_of_nodes(), "nodes,", prev_g.number_of_edges() * 2, "edges,", prev_SQ, "SQ", len(prev_scores), "paths")
			print("current g:", self.num_nodes, "nodes,", self.num_edges, "edges,", G1_SQ, "SQ", len(G1_scores), "paths")

			# repeat until no improvement
			# only add these random edges if it increases the score
			# otherwise continue
			if G1_SQ <= prev_SQ: 
				# update self to prev g
				self.update_all(prev_g, prev_paths, prev_paths_list)
				best_SQ = prev_SQ
				best_scores = prev_scores
				total_num_no_improv += 1
				num_no_improv += 1
				print('NO IMPROVEMENT; current best S(Q) (calculated using local lambda penalty):', best_SQ, ", num_no_improv:", num_no_improv, ', num_improv: ', num_improv)
				if improv_list != []:
					if num_no_improv == t or stats.mean(improv_list) < 5: # have converges if have less than 5% improvement on average or have t continuous no improvements
						# plot a curve of G vs G1 S(Q) (percentage improvement)
						self.global_SQ = best_SQ
						print(improv_list, "current # no improvement:", num_no_improv, "; total # no improvement:", total_num_no_improv)
						# save optimal graph object
						pickle_dump(self, osp.join(self.files_dir, '_'.join(['optimal_graph', str(l), str(N), str(t) + '.pickle'])))
						# with open('files/' + self.name + '/optimal_graph_' + str(l) + '_' + str(N) + '_' + str(t) + '.pickle', 'wb') as f:
						# 	pickle.dump(self, f)
						stop = timeit.default_timer()
						print('Time needed to find the optimal graph: ', stop - start)
						self.improv_list = improv_list
						self.improv_iterations = improv_iterations
						if len(improv_list) > 1:
							# plot improvement
							self.plot_improvement()
						return
						
			else: 
				# find percent improvement
				improv = (G1_SQ - prev_SQ)/abs(prev_SQ) * 100
				improv_list.append(improv)
				improv_iterations.append(j+1)
				best_SQ = G1_SQ
				best_scores = G1_scores
				num_no_improv = 0 # reset
				num_improv += 1
				self.global_lambda_pen = lambda_pen # update global_lambda pen when there is an improvement
				print('IMPROVEMENT:', improv, ' current best S(Q) (calculated using local lambda penalty):', best_SQ, '; num_improv:', num_improv)   
		
		# if reach max # of iterations
		# plot a curve of G vs G1 S(Q) (percentage improvement)
		self.global_SQ = best_SQ
		print(improv_list, "current # no improvement:", num_no_improv, "total # no improvement:", total_num_no_improv)
		# save optimal graph object
		pickle_dump(self, osp.join(self.files_dir, '_'.join(['optimal_graph', str(l), str(N), str(t) + '.pickle'])))
		# with open('files/' + self.name + '/optimal_graph_' + str(l) + '_' + str(N) + '_' + str(t) + '.pickle', 'wb') as f:
		# 	pickle.dump(self, f)
		stop = timeit.default_timer()
		print('Time needed to find the optimal graph: ', stop - start)
		self.improv_list = improv_list
		self.improv_iterations = improv_iterations
		return

	def calculate_path_score2(self, G, paths): # input parameters instead of attributes
		'''
		returns a dict of dict of path scores and a list of scores
		'''
		scores = {}
		scores_list = []
		for receptor in paths:
			scores[receptor] = {}
			for target in paths[receptor]:
				curr_paths = paths[receptor][target] # more than one shortest path
				# score = product of log2fc of all nodes x product of ppi score of all edges
				scores[receptor][target] = []
				for path in curr_paths:
					score = 0.0
					for i in range(len(path)-1):
						score += float(G.nodes[path[i]]['log2fc']) + float(G.edges[path[i], path[i+1]]['ppi']) 
					score +=  float(G.nodes[path[len(path) - 1]]['log2fc'])
					scores[receptor][target].append(score)
					scores_list.append(score)
			if not scores[receptor]: # empty
				del cores[receptor]
		return scores, scores_list

	# use for ranking nodes
	# removes a node from the graph
	def remove_node(self, node):
		original = copy.deepcopy(self.G)
		original.remove_node(node)
		new_paths_list = [path for path in self.paths_list if node not in path]
		new_paths = self.list_to_dict(new_paths_list)
		_, scores = self.calculate_path_score2(original, new_paths)
		SQ = (sum(scores) - self.global_lambda_pen * (sum(d for _, d in original.degree)/original.number_of_nodes()))
		return (node, SQ)

	# final step: rank nodes in the (near) global optimal graph
	def rank_nodes(self):
		'''
		returns a ranking of nodes from most important to least important
		'''
		_, scores = self.calculate_path_score()
		# case where didn't run optimal_graph
		if self.global_lambda_pen == 0.0:
			self.global_lambda_pen = self.local_lambda_pen
		global_SQ = (sum(scores) - self.global_lambda_pen * (sum(d for _, d in self.G.degree)/self.G.number_of_nodes()))
		print('Global S(Q) (calculated using global lambda penalty):', global_SQ)

		# use Pool for parallel processing
		# time amount of time it takes
		start = timeit.default_timer()
		p = Pool()
		SQ_scores = p.map(self.remove_node, self.nodes)
		p.close()
		p.join()
		stop = timeit.default_timer()
		print('Time needed to rank nodes: ', stop - start)

		SQ_score_differences = []
		for (gene, score) in SQ_scores:
			SQ_score_differences.append((gene, score - global_SQ))
		SQ_score_differences.sort(key = lambda x: x[1], reverse=False)
		SQ_scores.sort(key = lambda x: x[1], reverse=False)

		# write to files folder
		score_differences_f = osp.join(self.files_dir, '_'.join(['nodes_ranking_score_differences', str(self.l), str(self.N), str(self.t) + '.txt']))
		scores_f = osp.join(self.files_dir, '_'.join(['nodes_ranking_scores', str(self.l), str(self.N), str(self.t) + '.txt']))

		# score_differences_f = 'files/' + self.name + '/nodes_ranking_score_differences_' + str(self.l) + '_' + str(self.N) + '_' + str(self.t) + '.txt'
		# scores_f = 'files/' + self.name + '/nodes_ranking_scores_' + str(self.l) + '_' + str(self.N) + '_' + str(self.t) + '.txt'
		
		write_list_to_file_convert_to_str(SQ_scores, scores_f)
		write_list_to_file_convert_to_str(SQ_score_differences, score_differences_f)

		# with open(score_differences_f, 'w') as f_diff, open(scores_f, 'w') as f_scores:
		# 	for tup in SQ_scores:
		# 		f_scores.write(str(tup) + '\n')
		# 	for tup in SQ_score_differences:
		# 		f_diff.write(str(tup) + '\n')
		return (SQ_score_differences)

class SignificantGenes:
	def __init__(self, script_dir, g, n):
		self.g = g # (near) optimal graph
		self.n = n # number of significant nodes to keep
		self.nodes = []
		self.node_score = {} # key = nodes, value = score
		self.recs = []
		self.targets = []
		self.internal_nodes = []
		self.ranking = [] # top self.n ranking list
		self.ranking_dict = {} # key = gene, value = rank in list (1-100)
		self.files_dir = osp.join(self.script_dir, 'files', self.g.name)
		self.results_dir = osp.join(self.script_dir, 'results', self.g.name)
		
	# read in ranking from file and take top n
	def get_ranking(self):
		ranking_file = osp.join(self.files_dir, '_'.join(['nodes_ranking_score_differences', str(self.g.l), str(self.g.N), str(self.g.t) + '.txt']))
		# ranking_file = 'files/' + self.g.name + '/nodes_ranking_score_differences_' + str(self.g.l) + '_' + str(self.g.N) + '_' + str(self.g.t) + '.txt'
		i = 0
		num_added = 0
		rank = 1
		with open(ranking_file, 'r') as fname:
			for line in fname:
				line = literal_eval(line) # read in as a tuple
				if i < self.n:
					if line[1] > 0.0:
						print('Not enough nodes with negative score changes. So only', num_added, 'top significant genes were kept.')
						break
					# update self.ranking and self.nodes
					self.nodes.append(line[0])
					self.node_score[line[0]] = str(line[1])
					self.ranking.append(line)
					self.ranking_dict[line[0]] = rank
					rank += 1
					num_added += 1
				else:
					print(len(self.ranking), 'top significant genes were kept.')
					break
				i += 1
	
	# get files for GO term enrichment/pathway analyses
	# returns list of siginificant nodes, targets, recs, and internal nodes
	def get_go(self):
		# create directories and filenames
		output_genes_file = osp.join(self.results_dir, 'genes_go.txt') 
		output_receptors_file = osp.join(self.results_dir, 'receptors_go.txt')
		output_targets_file = osp.join(self.results_dir, 'targets_go.txt')
		output_internal_file = osp.join(self.results_dir, 'internal_nodes_go.txt')
		output_ranking_file = osp.join(self.results_dir, 'nodes_ranking.tsv')

		# output_genes_file = dirname + '/genes_go.txt'
		# output_receptors_file = dirname + '/receptors_go.txt'
		# output_targets_file = dirname + '/targets_go.txt'
		# output_internal_file = dirname + '/internal_nodes_go.txt'
		# output_ranking_file = dirname + '/nodes_ranking.tsv'

		# list of targets that are also source nodes
		source_target_overlap = []

		# write genes, receptors, TFs, and internal nodes to corresponding files
		with open(output_genes_file, 'w') as gfname, open(output_receptors_file, 'w') as rfname, open(output_targets_file, 'w') as tname, open(output_internal_file, 'w') as inname, open(output_ranking_file, 'w') as rname:
			rname.write('\t'.join(['Gene/Protein', 'Score_Change', 'Type']) + '\n')
			# get list of nodes and write to file
			for node in self.nodes:
				gfname.write(node + '\n')
				# get list of recs and write to file
				if node in self.g.recs:
					if node in self.g.targets:
						source_target_overlap.append(node)
					self.recs.append(node)
					rfname.write(node + '\n')
					rname.write('\t'.join([node, self.node_score[node], 'receptor']) + '\n')
				elif node in self.g.targets:
					tname.write(node + '\n')
					self.targets.append(node)
					rname.write('\t'.join([node, self.node_score[node], 'target']) + '\n')
				else: # internal node
					inname.write(node + '\n')
					self.internal_nodes.append(node)
					rname.write('\t'.join([node, self.node_score[node], 'internal node']) + '\n')
		print('Top significant genes, receptors, targets, and internal nodes have been saved for GO molecular function & Reactome Pathways analysis.')
		print(len(source_target_overlap), 'nodes that are both receptors (sources) and DE genes/TFs (targets) are designated as source nodes:', ', '.join(source_target_overlap))

	# add node size to .noa file for viewing on Cytoscape
	# nodes with higher rankings have larger node sizes
	# size 5 is the largest size
	# nodes ranked 1-20 are of size 5, 21-40 (size 4), 41-60 (size 3), 61-80 (size 2), 81-100 (size 1)
	def get_node_size(self, node):
		node_ranking_categories = [21, 41, 61, 81, 101]
		size = 5
		for node_ranking_category in node_ranking_categories:
			if self.ranking_dict[node] < node_ranking_category:
				return size
			size -= 1

	# get files for viewing (near) optimal graph details on Cytoscape			 
	def get_noa_sif(self):
		dirname = 'results/' + self.g.name
		noa_file = dirname + '/optimal_network.noa'
		sif_file = dirname + '/optimal_network.sif'

		# get list of edges
		with open(sif_file, 'w') as sif_fname:
			for node in self.nodes:
				curr_list = [n for n in self.g.G.neighbors(node) if n in self.nodes]
				curr_list = [node, 'ppi'] + curr_list
				sif_fname.write('\t'.join(curr_list) + '\n')
				
		# get list of attributes
		with open(noa_file, 'w') as noa_fname:
			noa_fname.write('\t'.join(['Name', 'Type', 'Node_Size']) + '\n')
			for rec in self.recs:
				size = str(self.get_node_size(rec))
				noa_fname.write('\t'.join([rec, 'source', size]) + '\n')
			for target in self.targets:
				size = str(self.get_node_size(target))
				noa_fname.write('\t'.join([target, 'target', size]) + '\n')
			# get internal nodes
			internal_nodes = list((set(self.g.nodes) - set(self.g.recs)) - set(self.g.targets))
			for node in self.internal_nodes:
				size = str(self.get_node_size(node))
				noa_fname.write('\t'.join([node, 'internal', size]) + '\n')
		print('.noa and .sif files have been saved for viewing optimal graph details on Cytoscape.')

