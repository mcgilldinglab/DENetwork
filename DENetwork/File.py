import abc
import h5py
import numpy as np
import pdb
import pandas as pd
import anndata
from scipy.sparse import csr_matrix,coo_matrix,lil_matrix
import scanpy as sc


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


class RFile(File):
	def __init__(self,fname):
		self.name=fname
	def read(self):
		import rpy2
		import rpy2.robjects as robjects
		from rpy2.robjects import pandas2ri
		readRDS=robjects.r['readRDS']
		df=readRDS(self.name)
		df = pandas2ri.rpy2py_dataframe(df)
		return df
		
class HugeMatrixFile:
	def __init__(self,fname):
		self.name=fname
	# read and save it to sparse matrix
	def read(self,sep,sIndex):
		f=open(self.name,'r')
		var=f.readline()
		var=var.strip().split(sep)[sIndex:]
		obs=[]
		
		row=[]
		col=[]
		data=[]
		
		line_ct=0
		
		for line in f:
			line=line.strip().split(sep)
			fc=line[:sIndex]
			obs.append(fc)
			ei=line[sIndex:]
			for j in range(len(ei)):
				if float(ei[j])!=0:
					row.append(line_ct)
					col.append(j)
					data.append(float(ei[j]))
			line_ct+=1
			if line_ct%1000==0:
				print(line_ct)
		M=len(obs)
		N=len(var)
		GM=coo_matrix((data,(row,col)),shape=(M,N))
		return [GM,obs,var]
		
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
				
	# conver the tab file to h5 file
	def toH5(self,sep,outname,cellLabels,chunksize=999):
		# sep: delimiter 
		# outname: output file name prefix
		# cellLabels: a list of attributes (string) for each of the cells
		# output file format: 
		# X: data matrix
		# var: gene annotations
		# obs: cell annotations
		
		si=len(cellLabels)
		f=open(self.name,'r')
		
		genes=f.readline().strip().split(sep)
		genes=genes[si:]
		nCols=len(genes)
		f.close()
		
		chunks=self.readchunk(sep,chunksize,True)
		hf=h5py.File(outname+'.h5','w')
		dset=hf.create_dataset('X', (chunksize,nCols),maxshape=(None,nCols),dtype='float',compression="gzip",compression_opts=9)
		
		row_count=0
		
		cells=[]
		chunk_ct=0
		
		for chunk in chunks:
			chunk_X=np.array([item[si:] for item in chunk],dtype='float')
			chunk_cells=[item[:si] for item in chunk]
			cells+=chunk_cells
			dset.resize(row_count+len(chunk),axis=0)
			dset[row_count:]=chunk_X 
			row_count+=len(chunk)
			chunk_ct+=1
			print("%s cells"%(chunk_ct*(chunksize+1)))
			
		
		dt = h5py.string_dtype()
		
		# writing var (genes~)
		genes=np.array(genes,dtype=h5py.string_dtype(encoding='utf-8'))
		gVars=hf.create_group('var')
		gVars.create_dataset('index',data=genes,dtype=dt)
		
		# writing obs (cells~)
		gObs=hf.create_group('obs')
		for i in range(len(cellLabels)):
			li=cellLabels[i]
			ci=np.array([item[i] for item in cells],dtype=h5py.string_dtype(encoding='utf-8'))
			gObs.create_dataset(li,data=ci,dtype=dt)
		
		hf.close()
		
		# read file to anndata object (sparse_matrix)


			
class H5File(File):
	def read(self):
		hfdata=h5py.File(self.name,'r')
		return hfdata
	
	# convert the full matrix to sparse matrix and store it to h5ad file 
	def toSparseAnnData(self,outname,BLOCK=5000):
		hfdata=h5py.File(self.name,'r')
		X=hfdata['X']
		row=np.empty(0,dtype=np.intc)
		col=np.empty(0,dtype=np.intc)
		dat=np.empty(0,dtype=np.float64)
		(M,N)=X.shape
		
		p=0
		while (p+BLOCK<M):
			XP=X[p:p+BLOCK]
			pcm=coo_matrix(XP)
			prow=[p+item for item in pcm.row]
			dat=np.append(dat,pcm.data)
			col=np.append(col,pcm.col)
			row=np.append(row,prow)
			p+=BLOCK
			print(p)
		XP=X[p:]
		pcm=coo_matrix(XP)
		prow=[p+item for item in pcm.row]
		dat=np.append(dat,pcm.data)
		col=np.append(col,pcm.col)
		row=np.append(row,prow)
		
		spm=csr_matrix((dat,(row,col)),shape=(M,N))
		xobs=pd.DataFrame(data=dict(hfdata['obs']))
		xobs=xobs.set_index(['index'])
		
		xvar=pd.DataFrame(data=dict(hfdata['var']))
		xvar=xvar.set_index(['index'])
		
		andata=anndata.AnnData(X=spm,obs=xobs,var=xvar)
		andata.write_h5ad(outname,compression=9)
		return andata
		
		
	def toAnnData(self):
		hfdata=self.read()
		xobs=pd.DataFrame(index=np.array(hfdata['obs']['index']))
		xvar=pd.DataFrame(index=np.array(hfdata['var']['index']))
		andata=anndata.AnnData(X=np.array(hfdata['X']),obs=xobs,var=xvar)
		return andata

class SC10xFile(File):
	def __init__(self,cname,gname,mname):
		self.cname=cname
		self.gname=gname
		self.nname=mname 
	
	
class H5ADFile(File):
	def __init__(self,fname):
		self.d1=anndata.read_h5ad(fname)
		
	def preprocess(self,PCT_MT=5):
		datControl=self.d1
		sc.pp.filter_cells(datControl, min_genes=200)
		sc.pp.filter_genes(datControl, min_cells=3) 
		datControl.var['mt'] = datControl.var_names.str.startswith('MT-')
		sc.pp.calculate_qc_metrics(datControl, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
		sc.pl.violin(datControl, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
				 jitter=0.4, multi_panel=True)
		datControl = datControl[datControl.obs.pct_counts_mt < PCT_MT, :]
		sc.pp.normalize_total(datControl, target_sum=1e4)
		sc.pp.log1p(datControl)
		#sc.pp.highly_variable_genes(datControl, min_mean=0.00005, max_mean=7, min_disp=0.05)
		#datControl = datControl[:, datControl.var.highly_variable]
		self.d1=datControl 
		
	def cluster(self):
		datControl=self.d1
		sc.tl.pca(datControl, svd_solver='arpack')
		sc.pp.neighbors(datControl, n_neighbors=10, n_pcs=40)
		sc.tl.diffmap(datControl)
		sc.tl.leiden(datControl)
		sc.tl.paga(datControl)
		sc.pl.paga(datControl)
		sc.tl.umap(datControl,init_pos='paga')
		sc.pl.umap(datControl,color='leiden',legend_loc="on data")
		self.d1=datControl
		
	def getDE(self):
		datControl=self.d1
		sc.tl.rank_genes_groups(datControl, 'leiden', method='wilcoxon')
		sc.pl.rank_genes_groups(datControl, n_genes=25, sharey=False)	
		deg=pd.DataFrame(datControl.uns['rank_genes_groups']['names'])
		deg.to_csv(f"{self.name}_topDEGenes.txt",sep="\t")
		datControl.layers['scaled'] = sc.pp.scale(datControl, copy=True).X
		sc.pl.rank_genes_groups_matrixplot(datControl, n_genes=3, use_raw=False, vmin=-3, vmax=3, cmap='bwr', layer='scaled',save="DE_heatmap")
		

class FastaFile(File):
	def read(self):
		f=open(self.name,'r')
		lf=f.readlines()
		f.close()
		lf="".join(lf)
		lf=lf.split('>')[1:]
		for i in range(len(lf)):
			lf[i]=lf[i].split('\n')
			lf[i]=[lf[i][0],''.join(lf[i][1:])]
		return lf

class TFBSFile(File):
	def read(self):
		f=open(self.name,'r')
		lf=f.readlines()
		f.close()
		lf="".join(lf)
		lf=lf.split('>')[1:]
		for i in range(len(lf)):
			lf[i]=lf[i].split('\n')
			lf[i]=[item for item in lf[i] if item!='']
		return lf
	
	
class LineFile(File):
	def read(self):
		f=open(self.name,'r')
		lf=f.readlines()
		f.close()
		lf=[item.strip() for item in lf]
		return lf
			

class MotifFile(File):
	def read(self):
		f=open(self.name,'r')
		lf=f.readlines()
		f.close()
		lf="".join(lf)
		lf=lf.split('>')[1:]
		for i in range(len(lf)):
			lf[i]=lf[i].split('\n')
			lf[i]=[item for item in lf[i] if item!='']
			lf[i]=[item.split() for item in lf[i]]
			lf[i][0]=lf[i][0][0]
		return lf

class MicroArrayFile(File):
	def __init__(self,FileList):
		self.FileList=FileList
		self.E=self.read()
		
	def read(self):
		out=[]
		for i in self.FileList:
			fi=TabFile(i).read('\t')
			ID=[item[0] for item in fi]
			fi=[item[1] for item in fi]
			out.append(fi)
		L=len(ID)
		ME=[]
		for i in range(L):
			ei=[item[i] for item in out]
			ME.append([ID[i]]+ei)
		return ME
		
	def toGeneSymbol(self,geneMap):
		gm=TabFile(geneMap).read('\t')
		dgm={}
		for i in gm:
			if len(i)>9:
				mi=i[9].split('//')
				if len(mi)>1:
					mi=mi[1].strip()
					dgm[i[0]]=mi
		for i in range(len(self.E)):
			if self.E[i][0] in dgm:
				self.E[i][0]=dgm[self.E[i][0]]
			else:
				self.E[i]=[]
		self.E=[item for item in self.E if item!=[]]
		return self.E
		
		
		
		
