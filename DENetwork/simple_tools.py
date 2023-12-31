'''

This script contains simple functions for use in other scripts.
Author: Ting-Yi Su ting-yi.su@mail.mcgill.ca

'''
import pickle
import os.path as osp
import os

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

# reads in items from newline delimited file
# and returns a list of the items
def get_list_from_file(newline_delimited_file):
	items = []
	with open(newline_delimited_file, 'r') as f:
		items = [line.strip() for line in f]
	return items

# for writing a list of strings to file
def write_list_to_file(list_to_write, file):
	with open(file, 'w') as f:
		for item in list_to_write:
			f.write(item + '\n')

# for writing a list of tuples, integers, etc to file
def write_list_to_file_convert_to_str(list_to_write, file):
	with open(file, 'w') as f:
		for item in list_to_write:
			f.write(str(item) + '\n')

# create dictionary with keys & lists as values
def initialize_dict_of_lists(keys):
	initialized_dict = {}
	for key in keys:
		initialized_dict[key] = []
	return initialized_dict

# reads in items from file with any type of delimiter
# and returns a list of the items
def get_list_from_delimited_file(file, sep):
	with open(file, 'r') as f:
		lines = f.readlines()
		list_of_items = [item.strip().split(sep) for item in lines]
	return list_of_items