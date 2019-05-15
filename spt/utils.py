import os
import networkit
import numpy as np
import nibabel as nib
from scipy.stats import mode
from itertools import combinations
from scipy.ndimage.morphology import binary_dilation

def pathFinder(data):
	edgelist_fn = data[0]
	voxel_coords = data[1]
	vpairs = data[2]
	reader = networkit.graphio.EdgeListReader(' ',0,directed=True)
	g = reader.read(edgelist_fn)
	sourceChange = np.append(np.where(np.roll(vpairs[:,0],1)!=vpairs[:,0])[0],vpairs.shape[0])
	sourceNodes = vpairs[sourceChange[:-1],0]
	vpairs = vpairs[:,1]
	pathStore = []

	for jj in range(sourceNodes.shape[0]):
		i = sourceNodes[jj]
		d = networkit.distance.Dijkstra(g, i)
		d.run()

		for ii in range(sourceChange[jj],sourceChange[jj+1],1):
			j = vpairs[ii]
			pathStore.append(voxel_coords[np.array(d.getPath(j))].astype(np.uint8))
                           
	return pathStore

def poolCN(data):
	edgelist_fn = data[0]
	parc = data[1]
	paths = data[2] 
	lut = data[3]
	reader = networkit.graphio.EdgeListReader(' ',0,directed=True)
	g = reader.read(edgelist_fn)
	cn = np.zeros((parc.max()+1,parc.max()+1))

	for path in paths:
		sL = path[0]
		tL = path[-1]				
		labels = parc[sL[0],sL[1],sL[2]],parc[tL[0],tL[1],tL[2]]
		w = [g.weight(lut[path[step][0],path[step][1],path[step][2]], lut[path[step+1][0],path[step+1][1],path[step+1][2]]) for step in range(len(path)-1)]
		prob = np.exp(-np.sum(w)/len(w))
		cn[labels[0],labels[1]] = cn[labels[0],labels[1]] + prob
		cn[labels[1],labels[0]] = cn[labels[0],labels[1]]
	return cn

class Utils(object):

	def writeEdgelist(self, edgelist_fn):
		# edgelist_fn: output filename for the edgelist
		newG = networkit.graphio.EdgeListWriter(' ',0,True)
		newG.write(self.mitns_graph.graph, edgelist_fn)
		self.edgelist_fn = edgelist_fn

	def parcellateWM(self, wm_path, gm_parc_path):
		# wm_path: path to layer 1 white matter mask
		lps_neighbor_shifts = {
		 'a': np.array([ 0, -1,  0]),
		 'ai': np.array([ 0, -1, -1]),
		 'as': np.array([ 0, -1,  1]),
		 'i': np.array([ 0,  0, -1]),
		 'l': np.array([1, 0, 0]),
		 'la': np.array([ 1, -1,  0]),
		 'lai': np.array([ 1, -1, -1]),
		 'las': np.array([ 1, -1,  1]),
		 'li': np.array([ 1,  0, -1]),
		 'lp': np.array([1, 1, 0]),
		 'lpi': np.array([ 1,  1, -1]),
		 'lps': np.array([1, 1, 1]),
		 'ls': np.array([1, 0, 1]),
		 'p': np.array([0, 1, 0]),
		 'pi': np.array([ 0,  1, -1]),
		 'ps': np.array([0, 1, 1]),
		 'r': np.array([-1,  0,  0]),
		 'ra': np.array([-1, -1,  0]),
		 'rai': np.array([-1, -1, -1]),
		 'ras': np.array([-1, -1,  1]),
		 'ri': np.array([-1,  0, -1]),
		 'rp': np.array([-1,  1,  0]),
		 'rpi': np.array([-1,  1, -1]),
		 'rps': np.array([-1,  1,  1]),
		 'rs': np.array([-1,  0,  1]),
		 's': np.array([0, 0, 1])}

		wm = nib.load(wm_path).get_data()[::-1,::-1,:]
		wm_coords = np.array(np.unravel_index(np.flatnonzero(wm.flatten(order="F")), wm.shape, order="F")).T  
		parc = nib.load(gm_parc_path).get_data()[::-1,::-1,:]
		parc[wm==1] = 0
		wm_parc = np.zeros(parc.shape, dtype=np.uint16)

		for j, starting_voxel in enumerate(wm_coords):
			nvals = []
			for i, name in enumerate(lps_neighbor_shifts):
				coord = tuple(starting_voxel + lps_neighbor_shifts[name])
				if parc[coord]: nvals.append(parc[coord])
			if len(nvals) > 0: wm_parc[tuple(starting_voxel)] = mode(nvals)[0][0] 
		dir_path = os.path.dirname(os.path.realpath(wm_path))
		new_path = dir_path + "/wm_layer1_parcellated.nii.gz"
		aff = nib.load(wm_path).affine
		new_nib = nib.Nifti1Image(wm_parc.astype(np.int16)[::-1,::-1,:], aff)
		nib.save(new_nib, new_path)
		self.parc = wm_parc.astype(np.int16) 

	def get_layer1_wm(self, gm_path, wm_path):
		gm = nib.load(gm_path).get_data()
		wm = nib.load(wm_path).get_data()
		aff = nib.load(wm_path).affine

		gm_dilated = binary_dilation(gm).astype(gm.dtype)
		new = (gm_dilated == 1) & (wm == 1)
		new[gm == 1] = 0
		new_nib = nib.Nifti1Image(new.astype(np.int16), aff)
		dir_path = os.path.dirname(os.path.realpath(wm_path))
		new_path = dir_path + "/wm_layer1.nii.gz"
		nib.save(new_nib, new_path)

	def FindStopNodes(self, k, seed_edge=None):
		g = nx.Graph()
		loss_cn = self.loss_cn
		full_dis = np.flatnonzero(np.triu(loss_cn)>=1.0)
	
		if seed_edge == None:
			if len(full_dis) > 1:
				edges = [(loss_cn[divmod(p, loss_cn.shape[1]),:].sum(),divmod(p, loss_cn.shape[1])) for p in full_dis]
				edges = sorted(edges)[::-1]
				first_edge = edges[-1][1]
		
			else: first_edge = divmod(loss_cn.argmax(), loss_cn.shape[1])
		else: first_edge = seed_edge
		dprofile = [loss_cn[first_edge]]
		g.add_edge(first_edge[0],first_edge[1], weight = loss_cn[first_edge])
		for ii in range(3, k+1):
			nodes = np.array(g.nodes())
			region_pool = loss_cn[nodes,:].sum(axis=0)
			region_pool[nodes] = 0
			new_node = region_pool.argmax()
			pre_weight = g.size(weight='weight')
			for n in nodes:
				g.add_edge(n, new_node, weight = loss_cn[n,new_node])

			post_weight = g.size(weight='weight')
			rat = (post_weight-pre_weight)
			dprofile.append(rat)
		return dprofile, g

	def savePaths(self, paths_fn):
		if self.paths: np.save(paths_fn, self.paths)
	def loadPaths(self, paths_fn):
		self.paths = np.load(paths_fn).tolist()
	def saveConnectome(self, connectome_fn):
		if self.cn: np.save(connectome_fn, self.cn)
	def loadConnectome(self, cn_fn):
		self.paths = np.load(cn_fn)
	def saveDisconnectome(self, disconnectome_fn):
		if self.loss_cn: np.save(disconnectome_fn, self.loss_cn)
	def loadDisconnectome(self, loss_cn_fn):
		self.paths = np.load(loss_cn_fn)
