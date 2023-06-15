# simple functions
import pickle
import os.path as osp
import os
from Bio.PDB import PDBList

def pickle_load(pickle_file):
	with open(pickle_file, 'rb') as f:
		contents = pickle.load(f)
		return contents

def pickle_dump(contents, pickle_file):
	with open(pickle_file, 'wb') as f:
		pickle.dump(contents, f)

# check whether a file/directory exists or not
def check_exist(dirname):
	return osp.exists(dirname)

# creates a directory if it doesn't exist already
def check_create_dir(dirname):
	if not osp.exists(dirname):
		os.mkdir(dirname)

# remove directory named dirname if it exists
def check_remove(dirname):
	if osp.exists(dirname):
		os.remove(dirname)

# def write_list_of_lists_to_tsv(contents, file):
# 	with open(file, 'w') as f:
# 		for line in contents:
# 			line = [str(item) for item in line]
# 			f.write('\t'.join(line) + '\n')

# reads in items from newline delimited file
# and returns a list of the items
def get_list_from_file(newline_delimited_file):
	items = []
	with open(newline_delimited_file, 'r') as f:
		items = [line.strip() for line in f]
	return items

def write_list_to_file(list_to_write, file):
	with open(file, 'w') as f:
		for item in list_to_write:
			f.write(item + '\n')

# create dictionary with keys & lists as values
def initialize_dict_of_lists(keys):
	initialized_dict = {}
	for key in keys:
		initialized_dict[key] = []
	return initialized_dict