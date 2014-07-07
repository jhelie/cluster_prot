Please cite: [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.10504.png)](http://dx.doi.org/10.5281/zenodo.10504)

cluster_prot
============

Python utility to detect protein clusters within MD simulations outputs (structure and/or trajectories)

usage
-----
python cluster_prot.py --help

requirements
------------
The following Python modules are needed:
- MDAnalysis
- scikit-learn (for DBSCAN implementation)

features
-----------
- user defined or automatic identification of proteins
- time evolution of protein clusters
- definition of size groups 
- statistics
- 1D and 2D graphs of proteins clustering
- output files for visualisation in VMD

To do
-----
- calculate_statistics
- plotting
- gro case


