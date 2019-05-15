import os
from setuptools import setup, find_packages

setup(
    name = "ShortestPathTools",
    version = "1.0",
    author = "Clint Greene",
    author_email = "clint@ece.ucsb.edu",
    description = ("Graph based tools for studying brain connectivity."),
    keywords = "graphs, connectomes, brain, dMRI, tractography, brain networks, structural networks, shortest paths, diffusion MRI",
    url = "",
      install_requires=[],
#          "scipy",
#          "numpy",
#          "nibabel",
#          "networkx",
#          "networkit",
#          "nilearn",
#          "multiprocessing"],
    packages=find_packages(),
)
