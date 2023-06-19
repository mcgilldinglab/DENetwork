'''

This script contains functions to process tab-delimited files
Author: Jun Ding jun.ding@mcgill.ca

'''	
class File:
	def __init__(self,fname):
		self.name=fname
	def read(self):
		raise NotImplementedError('Please implement this method')
	def Linux2Win(self):
		f=open(self.name,'r')
		lf=f.readlines()
		f.close()
		lf=[item.strip()+'\r\n' for item in lf]
		lf=''.join(lf)
		f=open(self.name+'.win.txt','w')
		f.write(lf)
		f.close()

class TabFile(File):
	def read(self,sep):
		f=open(self.name,'r')
		lf=f.readlines()
		f.close()
		lf=[item.strip().split(sep) for item in lf]
		return lf
	
	# read the file and build a dictionary (col1->col2)
	def read2dict(self,sep):
		lf=self.read("\t")
		dictf={}
		for i in lf:
			dictf[i[0]]=i[1]
		return dictf
		
	# replace sep1 with sep2
	def convert(self,sep1, sep2):
		infile=self.read(sep1)
		outfile=[]
		for i in infile:
			outfile.append(sep2.join(i))
		return outfile
	
		
	# filter columns using the index list provided in colList
	def filterColRead(self,sep,colList):
		infile=self.read(sep)
		if type(colList)==type([]):
			outfile=[[item[k] for k in colList] for item in  infile]
		elif type(colList)==type(1):
			outfile=[item[colList:] for item in infile]
		else:
			print("input colList must be a integer or  a list of integer")
			return None
		return outfile
	
	# a read generator that reads the file with specified chunk size
	def readchunk(self,sep,chunksize=999,header=True):
		f=open(self.name,'r')
		header=f.readline()
		
		ct=0
		chunkLines=[]
		for line in f:
			line=line.strip().split(sep)
			if ct<chunksize:
				chunkLines.append(line)
				ct+=1
			else:
				chunkLines.append(line)
				yield chunkLines
				chunkLines=[]
				ct=0
		yield chunkLines
				




		
		
		
