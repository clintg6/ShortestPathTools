import networkit
import numpy as np
import networkx as nx
import nibabel as nib
import multiprocessing as mp
from scipy.interpolate import UnivariateSpline
from utils import poolCN, pathFinder, Utils

numProcs = mp.cpu_count()
pool = mp.Pool(processes=numProcs)

class SPT(Utils):
	def __init__(self, mitns_graph="", parc_fn="", edgelist_fn=""):
		self.mitns_graph = mitns_graph
		self.g = mitns_graph.graph
		self.lut = mitns_graph.coordinate_lut
		self.voxel_coords = mitns_graph.voxel_coords
		if parc_fn != "": self.parc = nib.load(parc_fn).get_data()[::-1,::-1,:]
		self.edgelist_fn = edgelist_fn

	def computeConnectome(self):
		if not hasattr(self, 'paths'): 
			print("Run findPaths first.")
		else:
			pathSplit = np.array_split(self.paths, numProcs)

			data = [(self.edgelist_fn,self.parc,paths,self.lut) for paths in pathSplit]
			results = pool.map_async(poolCN, data)
			ret = results.get()	
			cn = np.stack(ret,2).sum(2)

			self.cn = cn

	def computeDisconnectome(self, lesion_fn=""):
		if not hasattr(self, 'cn'): 
			print("Run computeConnectome first.")
		else:
			g = self.g

			parc = self.parc
			paths = self.paths
			inlesion_cn = np.zeros((parc.max()+1,parc.max()+1))

			lut = self.lut

			lesion = nib.load(lesion_fn).get_data()[::-1,::-1,:]
			lesion_coords = np.array(np.unravel_index(np.flatnonzero(lesion.flatten(order="F")), lesion.shape, order="F")).T
			lesion_coords = np.stack([coord for coord in lesion_coords if tuple(coord) in lut])
			lesion_inds = np.stack([lut[tuple(coord)] for coord in lesion_coords if tuple(coord) in lut])
			lesion_inds = set(lesion_inds)

			for path in paths:
				path_inds = [lut[path[step][0],path[step][1],path[step][2]] for step in range(len(path))]
	
				if not lesion_inds.isdisjoint(path_inds): 
					sL = path[0]
					tL = path[-1]
					labels = parc[sL[0],sL[1],sL[2]],parc[tL[0],tL[1],tL[2]]
					w = [g.weight(path_inds[step], path_inds[step+1]) for step in range(len(path_inds)-1)]
					prob = np.exp(-np.sum(w)/len(w))
					inlesion_cn[labels[0],labels[1]] = inlesion_cn[labels[0],labels[1]] + prob
					inlesion_cn[labels[1],labels[0]] = inlesion_cn[labels[0],labels[1]]
			self.inlesion_cn = inlesion_cn
			self.loss_cn = np.nan_to_num(inlesion_cn/self.cn)

	def findPaths(self):
		if not hasattr(self, 'subSamples'): 
			print("Run getsubSamples first.")
		else:
			subSamplesSplit = np.array_split(self.subSamples, numProcs)
			data = [(self.edgelist_fn,self.voxel_coords,samples) for samples in subSamplesSplit]
			results = pool.map_async(pathFinder, data)
			ret = results.get()	
			self.paths = [p for paths in ret for p in paths]

	def getsubSamples(self):
		parc = self.parc
		labels2check = list(combinations(range(1,parc.max()+1),2))
		numRegions = parc.max()
		label_lut = {}
		for jj in range(1,numRegions+1):
			k0 = (parc == jj)
			coords = np.array(np.unravel_index(np.flatnonzero(k0.flatten(order="F")), parc.shape, order="F")).T
			inds = [mitns.coordinate_lut[tuple(c)] for c in coords]
			label_lut[jj] = inds
		subSamples = []
		for pair in labels2check:
			numInds = min([len(label_lut[pair[0]]), len(label_lut[pair[1]])])
			p0 = np.random.choice(label_lut[pair[0]],(numInds,),False)
			p1 = np.random.choice(label_lut[pair[1]],(numInds,),False)
			subP = np.vstack((p0,p1)).T
			subSamples.append(subP.tolist())
		self.subSamples = np.array(sorted([p for s in subSamples for p in s]))

	def findMDS(self, seed_edge=None, smooth=0.1):
		if not hasattr(self, 'loss_cn'):
			print("Run computeDisconnectome first.")
		else: 
			k = self.parc.max()
			x = np.arange(k-1)
			dprofile, subG = self.FindStopNodes(k, seed_edge)
			spline = UnivariateSpline(x, dprofile, s=smooth)
			self.dprofile = spline(x)
			koptimal = np.argmax(spline(x)) + 2
			self.koptimal = koptimal
			dprofile, subG = self.FindStopNodes(koptimal, seed_edge)
			self.MDS = subG

