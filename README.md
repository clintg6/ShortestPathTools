# Shortest Path Brain Connectivity Tools 

This package contains fast, accurate, and robust tools for studying white matter connectivity. Tools include methods for rapidly generating shortest paths for the whole brain in white matter voxel graphs, parcellating white matter at the gray-white interface with a cortical parcellation, constructing shortest-path probability weighted structural connectomes and disconnectomes, and finding maximally disconnected subnetworks in disconnectomes. If you find any of our tools useful in your research, please cite us: [Coming soon.](https://www.mitpressjournals.org/doi/abs/10.1162/NETN_a_00035)

## Installation
Shortest path tools requires MITTENS, networkx, networkit, numpy, nibabel, scipy. MITTENS can be installed by following the instructions at: [https://github.com/mattcieslak/MITTENS/](https://github.com/mattcieslak/MITTENS/)The others can be easily installed using: pip install networkx networkit numpy, nibabel scipy. To install our shortest path tools, download the current version, enter its directory, and run pip to install it:

```bash
$ git clone https://github.com/clintg6/ShortestPathTools.git
$ cd ShortestPathTools
$ pip install -e .
```
To begin using our shortest path tools, a MITTENS Voxel Graph is needed.
For further instruction on using MITTENS to construct a Voxel Graph [click here](https://github.com/mattcieslak/MITTENS): 

## Building a Voxel Graph with MITTENS
```
$ python
from mittens import MITTENS
nifti_mitns = MITTENS(nifti_prefix="YOUR_PREFIX_HERE")
voxel_graph = nifti_mitns.build_graph(doubleODF=True, weighting_scheme='negative_log_p')

```
This loads the doubleODF transition probability nifti files output by MITTENS and builds a voxel graph with the -ln(transition probability) as the edge weights.

## Parcellating a single layer of white matter
A single layer of white matter at the gray-white interface needs to be parcellated to properly define the start and end points for Dijkstra's algorithm to find the shortest path. To create this kind of a parcellation a couple things are necessary before running the following commands: a gray matter and white matter mask and a gray matter parcellation. Once you have these files ready, run the following commands to create the white matter parcellation.

```
from spt import SPT
spt = SPT(voxel_graph)
spt.get_layer1_wm(gm_path,wm_path)
spt.parcellateWM(wm_path, gm_parc_path)
```
If you have already created a single layer white matter parcellation, you can enter the following command to prepare for the shortest path computation.

```
from spt import SPT
parc_fn = "full path to parcellation file"
spt = SPT(voxel_graph, parc_fn)
```

## Save voxel graph edgelist
To enable the shortest paths to be computed in a parallel manner, you must first save the voxel graph out in edgelist format so that it can be loaded by the node/pool workers because the voxel graph cannot be pickled. To do this simply run the following command:
```
edgelist_fn = "full path to where edge list is to be saved"
spt.writeEdgelist(edgelist_fn)
```

If you have already created an edgelist, you can enter the following command to prepare for the shortest path computation.

```
from spt import SPT
spt = SPT(voxel_graph, parc_fn, edgelist_fn)
```
## Finding Shortest Paths
Now that the Voxel Graph and its edgelist and the white matter parcellation have been created, they can be input into our toolbox to generate whole brain shortest paths and construct shortest path probability weighted connectomes and disconnectomes.

```
spt.getsubSamples() # generates the unique subsampled voxel source target pairs for every possible combination of cortical region pairs
spt.findPaths() # find the shortest paths between the source target pairs 
```

## Constructing a connectome
To construct a shortest-path probability weighted connectome run the following command:
```
spt.computeConnectome()
```
## Constructing a disconnectome
To construct a disconnectome, run the following command and specify the path to a patient's lesion mask.
```
spt.computeDisconnectome()
```

## Finding Maximally Disconnected Subnetworks
To find the maximally  disconnected subgraph that contains the set of cortical regions that share the greatest disconnectivity due to the lesion, run the following command.
```
spt.findMDS() 
```
Cite
========
[Coming soon.](https://www.mitpressjournals.org/doi/abs/10.1162/NETN_a_00035)

Credits
========
This source code was sponsored by a grant from the GE/NFL head health challenge. 
The content of the information does not necessarily reflect the position or
the policy of these companies, and no official endorsement should be inferred.

Authors
-------
 * Clint Greene


