################################################################################################################################################
# IMPORT MODULES
################################################################################################################################################

#import general python tools
import argparse
import itertools
import operator
from operator import itemgetter
import sys, os, shutil
import os.path
import math

#import python extensions/packages to manipulate arrays
import numpy 				#to manipulate arrays
import scipy 				#mathematical tools and recipesimport MDAnalysis

#import graph building module
import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
import matplotlib.cm as cm			#colours library
import matplotlib.colors as mcolors
import matplotlib.ticker
from matplotlib.ticker import MaxNLocator
from matplotlib.font_manager import FontProperties
fontP=FontProperties()

#import clustering algorithms
from sklearn.cluster import DBSCAN
import networkx as nx

#import MDAnalysis
import MDAnalysis
from MDAnalysis import *
import MDAnalysis.analysis
import MDAnalysis.analysis.leaflet
import MDAnalysis.analysis.distances

#set MDAnalysis to use periodic boundary conditions
MDAnalysis.core.flags['use_periodic_selections'] = True
MDAnalysis.core.flags['use_KDTree_routines'] = False

################################################################################################################################################
# RETRIEVE USER INPUTS
################################################################################################################################################

#create parser
#=============
version_nb="0.1.2"
parser = argparse.ArgumentParser(prog='cluster_prot', usage='', add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, description=\
'''
*******************************************
v''' + version_nb + '''
author: Jean Helie
git: https://github.com/jhelie/cluster_prot
*******************************************
	
[ Description ]

This script identifies protein clusters in a gro file or throughout a trajectory
using either a connectivity based or a density based (DBSCAN) algorithm on the
x,y,z coordinates of lipids headgroups.
	
It produces 3 types of outputs (only the 1st can be obtained without specifying groups, see note 5):
 - 2D plots: time evolution of the cluster size each protein is involved in
 - 1D plots: time evolution of the % of protein represented by each size group
 - stability statistics: maximun nb of consecutive frames each size group existed for

[ Requirements ]

The following python modules are needed:
 - MDAnalysis
 - scikit learn (python-sklearn), for density based clustering

[ Notes ]

1. It's a good idea to trjconv the xtc first and only outputs the proteins,  as the 
   script will run MUCH faster. Also, use the -pbc mol option.	

2. Two clustering algorithms can be used, the most appropriate parameters depends on the protein
   of interest:
   -connectivity: based on networkX, a protein is considered in a cluster if its within a distance
                  less than --cutoff from another protein. This means that a single protein can 
                  act as a connectory between two otherwise disconnected protein clusters. This
				  algorithm can be ran using either the distance between the center of geometry of the
				  proteins ('cog') or the minimum distante between proteins ('min').
				  The 'min' option scales as the square of the number of proteins in the system and
				  is thus very slow for large systems.
   -density: based on the DBSCAN algorithm implemented in scikit, a protein is considered in a cluster
             if is surrounded by at least --neighbours other proteins within a radius of --radius.
             This density based approach is usually less suited to the detection of protein
	         clusters but as a general rule the more compact the clusters, the smaller --radius and
	         the higher --neighbours can be - for details on this algorithm see its online
	         documentation. This algorithm is selected by setting the --algorithm option to 'density'.

3. Clusters statistics can be binned by defining size groups (-g).
   Size groups should be specified in a file where each line has the format:
       'lower_group_size,upper_group_size, colour'
   and respect the following rules:
    - to specify an open ended group use 'max', e.g. '3,max,color'
    - groups should be ordered by increasing size and their boundaries should not overlap
    - boundaries are inclusive so you can specify one size groups with 'size,size,color'
    - colours MUST be specified, either as single letter code, hex code  or colormap name (see note 6)
    - in case a colormap is used its name must be specified as the color of each cluster
    - any cluster size not fallig within the specified size groups will be labeled as 'other' and
      coloured in grey (#C0C0C0).
	
4. Colours of individual cluster sizes use the matplotlib 'jet' color scheme and cannot be
   modified. To use something else just specify a group file! (see above note).

5. If you don't want to use custom colours for your groups, you can use the name of a standard
   matplotlib color map - don't forget it should be specified on each line of your group file.
   Type 'cluster_prot --colour_maps' to see a list of their names.

6. Proteins are detected automatically but you can specify an input file to define your own selection
   with the -p option. In this case this file should contain on each line a protein selection string
   that can be passed as the argument of the MDAnalysis 'selectAtoms()' routine (e.g. 'bynum 1:344').

7. The size (or size group) of the cluster each protein is detected to be involved can be visualised
   with VMD. This can be done either with pdb files (output frequency controled via -w flag) or with 
   the xtc trajectory.
     - pdb file: the clustering info for each protein is stored in the beta factor column. Just open
                 the pdb with VMD and choose Draw Style > Coloring Method > Beta 
     - xtc file: the clustering info is stored in a .txt file in /4_VMD/ and you can load it into the
                 user field in the xtc by sourcing the script 'set_user_fields.tcl' and running the 
                 procedure 'set_cluster_prot'

[ Usage ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: structure file [.gro]
-x			: trajectory file [.xtc] (optional)
-o			: name of output folder
-g			: cluster groups definition file, see note 3
-p			: protein selection file, (optional, see note 6)
-b			: beginning time (ns)
-e			: ending time (ns)	
-t 		10	: process every t-frames
-w			: write annotated pdbs every [w] processed frames (optional, see note 7)
--smooth		: nb of points to use for data smoothing (optional)
--algorithm	cog	: 'cog','min' or 'density', see note 2

Algorithm options (see note 3)
-----------------------------------------------------
--cutoff 	8	: networkX cutoff distance for lipid-lipid contact (Angtrom)
--radius 	20	: DBSCAN search radius (Angtrom)
--neighbours 	3	: DBSCAN minimum number of neighbours within a circle of radius --radius	
 
Other options
-----------------------------------------------------
--colour_maps		: show list of standard colour maps, see note 5
--version		: show version number and exit
-h, --help		: show this menu and exit
  
''')

#data options
parser.add_argument('-f', nargs=1, dest='grofilename', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-x', nargs=1, dest='xtcfilename', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-o', nargs=1, dest='output_folder', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-g', nargs=1, dest='cluster_groups_file', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-p', nargs=1, dest='selection_file_prot', default=['auto'], help=argparse.SUPPRESS)
parser.add_argument('-b', nargs=1, dest='t_start', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-e', nargs=1, dest='t_end', default=[10000000000000], type=int, help=argparse.SUPPRESS)
parser.add_argument('-t', nargs=1, dest='frames_dt', default=[10], type=int, help=argparse.SUPPRESS)
parser.add_argument('-w', nargs=1, dest='frames_write_dt', default=[1000000000000000], type=int, help=argparse.SUPPRESS)
parser.add_argument('--algorithm', dest='m_algorithm', choices=['cog','min','density'], default='cog', help=argparse.SUPPRESS)
parser.add_argument('--smooth', nargs=1, dest='nb_smoothing', default=[0], type=int, help=argparse.SUPPRESS)
#algorithm options
parser.add_argument('--cutoff', nargs=1, dest='cutoff_connect', default=[8], type=float, help=argparse.SUPPRESS)
parser.add_argument('--radius', nargs=1, dest='dbscan_dist', default=[20], type=float, help=argparse.SUPPRESS)
parser.add_argument('--neighbours', nargs=1, dest='dbscan_nb', default=[3], type=int, help=argparse.SUPPRESS)
#other options
parser.add_argument('--colour_maps', dest='show_colour_map', action='store_true', help=argparse.SUPPRESS)
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#store inputs
#============
args=parser.parse_args()
args.grofilename=args.grofilename[0]
args.xtcfilename=args.xtcfilename[0]
args.output_folder=args.output_folder[0]
args.cluster_groups_file=args.cluster_groups_file[0]
args.selection_file_prot=args.selection_file_prot[0]
args.t_start=args.t_start[0]
args.t_end=args.t_end[0]
args.frames_dt=args.frames_dt[0]
args.frames_write_dt=args.frames_write_dt[0]
args.nb_smoothing=args.nb_smoothing[0]
args.dbscan_dist=args.dbscan_dist[0]
args.dbscan_nb=args.dbscan_nb[0]
args.cutoff_connect=args.cutoff_connect[0]

#show colour maps
#----------------
if args.show_colour_map:
	print ""
	print "The following standard matplotlib color maps can be used:"
	print ""
	print "Spectral, summer, coolwarm, pink_r, Set1, Set2, Set3, brg_r, Dark2, hot, PuOr_r, afmhot_r, terrain_r,"
	print "PuBuGn_r, RdPu, gist_ncar_r, gist_yarg_r, Dark2_r, YlGnBu, RdYlBu, hot_r, gist_rainbow_r, gist_stern, "
	print "gnuplot_r, cool_r, cool, gray, copper_r, Greens_r, GnBu, gist_ncar, spring_r, gist_rainbow, RdYlBu_r, "
	print "gist_heat_r, OrRd_r, CMRmap, bone, gist_stern_r, RdYlGn, Pastel2_r, spring, terrain, YlOrRd_r, Set2_r, "
	print "winter_r, PuBu, RdGy_r, spectral, flag_r, jet_r, RdPu_r, Purples_r, gist_yarg, BuGn, Paired_r, hsv_r, "
	print "bwr, cubehelix, YlOrRd, Greens, PRGn, gist_heat, spectral_r, Paired, hsv, Oranges_r, prism_r, Pastel2, "
	print "Pastel1_r, Pastel1, gray_r, PuRd_r, Spectral_r, gnuplot2_r, BuPu, YlGnBu_r, copper, gist_earth_r, "
	print "Set3_r, OrRd, PuBu_r, ocean_r, brg, gnuplot2, jet, bone_r, gist_earth, Oranges, RdYlGn_r, PiYG,"
	print "CMRmap_r, YlGn, binary_r, gist_gray_r, Accent, BuPu_r, gist_gray, flag, seismic_r, RdBu_r, BrBG, Reds,"
	print "BuGn_r, summer_r, GnBu_r, BrBG_r, Reds_r, RdGy, PuRd, Accent_r, Blues, Greys, autumn, cubehelix_r, "
	print "nipy_spectral_r, PRGn_r, Greys_r, pink, binary, winter, gnuplot, RdBu, prism, YlOrBr, coolwarm_r,"
	print "rainbow_r, rainbow, PiYG_r, YlGn_r, Blues_r, YlOrBr_r, seismic, Purples, bwr_r, autumn_r, ocean,"
	print "Set1_r, PuOr, PuBuGn, nipy_spectral, afmhot."
	print ""
	sys.exit(0)

#sanity check
#============
if not os.path.isfile(args.grofilename):
	print "Error: file " + str(args.grofilename) + " not found."
	sys.exit(1)
if args.cluster_groups_file!="no" and not os.path.isfile(args.cluster_groups_file):
	print "Error: file " + str(args.cluster_groups_file) + " not found."
	sys.exit(1)
if args.xtcfilename=="no":
	if '-t' in sys.argv:
		print "Error: -t option specified but no xtc file specified."
		sys.exit(1)
	elif '-b' in sys.argv:
		print "Error: -b option specified but no xtc file specified."
		sys.exit(1)
	elif '-e' in sys.argv:
		print "Error: -e option specified but no xtc file specified."
		sys.exit(1)
	elif '-w' in sys.argv:
		print "Error: -w option specified but no xtc file specified."
		sys.exit(1)
	elif '--smooth' in sys.argv:
		print "Error: --smooth option specified but no xtc file specified."
		sys.exit(1)
elif not os.path.isfile(args.xtcfilename):
	print "Error: file " + str(args.xtcfilename) + " not found."
	sys.exit(1)
if args.m_algorithm!="density":
	if '--radius' in sys.argv:
		print "Error: --radius option specified but --algorithm option set to '" + str(args.m_algorithm) + "'."
		sys.exit(1)
	elif '--neighbours' in sys.argv:
		print "Error: --neighbours option specified but --algorithm option set to '" + str(args.m_algorithm) + "'."
		sys.exit(1)
else:
	if '--cutoff' in sys.argv:
		print "Error: --cutoff option specified but --algorithm option set to 'density'."
		sys.exit(1)

#create folders and log file
#===========================
if args.output_folder=="no":
	if args.xtcfilename=="no":
		args.output_folder="cluster_prot_" + args.grofilename[:-4]
	else:
		args.output_folder="cluster_prot_" + args.xtcfilename[:-4]
if os.path.isdir(args.output_folder):
	print "Error: folder " + str(args.output_folder) + " already exists, choose a different output name via -o."
	sys.exit(1)
else:
	#create folders
	#--------------
	os.mkdir(args.output_folder)
	#1 sizes
	os.mkdir(args.output_folder + "/1_sizes")
	if args.xtcfilename!="no":
		os.mkdir(args.output_folder + "/1_sizes/1_1_plots_2D")
		os.mkdir(args.output_folder + "/1_sizes/1_1_plots_2D/png")
		os.mkdir(args.output_folder + "/1_sizes/1_2_plots_1D")
		os.mkdir(args.output_folder + "/1_sizes/1_2_plots_1D/png")
		os.mkdir(args.output_folder + "/1_sizes/1_2_plots_1D/xvg")
		os.mkdir(args.output_folder + "/1_sizes/1_3_biggest")
		os.mkdir(args.output_folder + "/1_sizes/1_3_biggest/png")
		os.mkdir(args.output_folder + "/1_sizes/1_3_biggest/xvg")
		if args.nb_smoothing>1:
			os.mkdir(args.output_folder + "/1_sizes/1_4_plots_1D_smoothed")
			os.mkdir(args.output_folder + "/1_sizes/1_4_plots_1D_smoothed/png")
			os.mkdir(args.output_folder + "/1_sizes/1_4_plots_1D_smoothed/xvg")
	#2 groups
	if args.cluster_groups_file!="no":
		os.mkdir(args.output_folder + "/2_groups")
		if args.xtcfilename!="no":
			os.mkdir(args.output_folder + "/2_groups/2_1_plots_2D")
			os.mkdir(args.output_folder + "/2_groups/2_1_plots_2D/png")
			os.mkdir(args.output_folder + "/2_groups/2_2_plots_1D")
			os.mkdir(args.output_folder + "/2_groups/2_2_plots_1D/png")
			os.mkdir(args.output_folder + "/2_groups/2_2_plots_1D/xvg")
			if args.nb_smoothing>1:
				os.mkdir(args.output_folder + "/2_groups/2_3_plots_1D_smoothed")
				os.mkdir(args.output_folder + "/2_groups/2_3_plots_1D_smoothed/png")
				os.mkdir(args.output_folder + "/2_groups/2_3_plots_1D_smoothed/xvg")
	#3 snapshots
	os.mkdir(args.output_folder + "/3_snapshots")
	os.mkdir(args.output_folder + "/3_snapshots/sizes")
	if args.cluster_groups_file!="no":
		os.mkdir(args.output_folder + "/3_snapshots/groups")	
	#4 VMD
	if args.xtcfilename!="no":
		os.mkdir(args.output_folder + "/4_VMD")	
	
	#create log
	#----------
	filename_log=os.getcwd() + '/' + str(args.output_folder) + '/cluster_prot.log'
	output_log=open(filename_log, 'w')		
	output_log.write("[cluster_prot v" + str(version_nb) + "]\n")
	output_log.write("\nThis folder and its content were created using the following command:\n\n")
	tmp_log="python cluster_prot.py"
	for c in sys.argv[1:]:
		tmp_log+=" " + c
	output_log.write(tmp_log + "\n")
	output_log.close()
	#copy input files
	#----------------
	if args.cluster_groups_file!="no":
		shutil.copy2(args.cluster_groups_file,args.output_folder + "/")

################################################################################################################################################
# DICTIONARIES
################################################################################################################################################

#color maps dictionaries
colormaps_possible=['Spectral', 'summer', 'coolwarm', 'pink_r', 'Set1', 'Set2', 'Set3', 'brg_r', 'Dark2', 'hot', 'PuOr_r', 'afmhot_r', 'terrain_r', 'PuBuGn_r', 'RdPu', 'gist_ncar_r', 'gist_yarg_r', 'Dark2_r', 'YlGnBu', 'RdYlBu', 'hot_r', 'gist_rainbow_r', 'gist_stern', 'gnuplot_r', 'cool_r', 'cool', 'gray', 'copper_r', 'Greens_r', 'GnBu', 'gist_ncar', 'spring_r', 'gist_rainbow', 'RdYlBu_r', 'gist_heat_r', 'OrRd_r', 'CMRmap', 'bone', 'gist_stern_r', 'RdYlGn', 'Pastel2_r', 'spring', 'terrain', 'YlOrRd_r', 'Set2_r', 'winter_r', 'PuBu', 'RdGy_r', 'spectral', 'flag_r', 'jet_r', 'RdPu_r', 'Purples_r', 'gist_yarg', 'BuGn', 'Paired_r', 'hsv_r', 'bwr', 'cubehelix', 'YlOrRd', 'Greens', 'PRGn', 'gist_heat', 'spectral_r', 'Paired', 'hsv', 'Oranges_r', 'prism_r', 'Pastel2', 'Pastel1_r', 'Pastel1', 'gray_r', 'PuRd_r', 'Spectral_r', 'gnuplot2_r', 'BuPu', 'YlGnBu_r', 'copper', 'gist_earth_r', 'Set3_r', 'OrRd', 'PuBu_r', 'ocean_r', 'brg', 'gnuplot2', 'jet', 'bone_r', 'gist_earth', 'Oranges', 'RdYlGn_r', 'PiYG', 'CMRmap_r', 'YlGn', 'binary_r', 'gist_gray_r', 'Accent', 'BuPu_r', 'gist_gray', 'flag', 'seismic_r', 'RdBu_r', 'BrBG', 'Reds', 'BuGn_r', 'summer_r', 'GnBu_r', 'BrBG_r', 'Reds_r', 'RdGy', 'PuRd', 'Accent_r', 'Blues', 'Greys', 'autumn', 'cubehelix_r', 'nipy_spectral_r', 'PRGn_r', 'Greys_r', 'pink', 'binary', 'winter', 'gnuplot', 'RdBu', 'prism', 'YlOrBr', 'coolwarm_r', 'rainbow_r', 'rainbow', 'PiYG_r', 'YlGn_r', 'Blues_r', 'YlOrBr_r', 'seismic', 'Purples', 'bwr_r', 'autumn_r', 'ocean', 'Set1_r', 'PuOr', 'PuBuGn', 'nipy_spectral', 'afmhot']

################################################################################################################################################
# DATA LOADING
################################################################################################################################################

# Create size groups
#===================

groups_number=0
groups_sizes_dict={}
groups_boundaries=[]
groups_colors_nb=0
groups_colors_dict={}
groups_colors_list=[]
groups_colors_map="custom"

if args.cluster_groups_file!="no":
	#read group definition file
	print "\nReading cluster groups definition file..."
	with open(args.cluster_groups_file) as f:
		lines = f.readlines()
	groups_number=len(lines)
	for g_index in range(0,groups_number):
		l_content=lines[g_index].split(',')
		tmp_beg=int(l_content[0])
		tmp_end=l_content[1]
		groups_colors_dict[g_index]=l_content[2][:-1]					#[:-1] to get rid of the final '\n' character
		if tmp_end=="max":
			tmp_end=100000												#put a stupidly big size to cap the open ended group (might beed increasing for super-large systems...)
		else:
			tmp_end=int(tmp_end)
		groups_boundaries.append([tmp_beg,tmp_end])
		
	#check for boundaries overlapping
	prev_beg=groups_boundaries[0][0]
	prev_end=groups_boundaries[0][1]
	if prev_end<prev_beg:
		print "Error: the upper boundary is smaller than the lower boundary for specified cluster groups" + str(g) + "."
		sys.exit(1)
	for g in groups_boundaries[1:]:
		if g[1]<g[0]:
			print "Error: the upper boundary is smaller than the lower boundary for group " + str(g) + "."
			sys.exit(1)
		if g[0]<=prev_end:
			print "Error: specified cluster groups [" + str(prev_beg) + "," + str(prev_end) + "] and " + str(g) + " overlap or are not in increasing order."
			sys.exit(1)
		prev_beg=g[0]
		prev_end=g[1]
	
	#check if a custom color map has been specified or not
	if groups_number>1 and len(numpy.unique(groups_colors_dict.values()))==1:
		if numpy.unique(groups_colors_dict.values())[0] in colormaps_possible:
			groups_colors_map=numpy.unique(groups_colors_dict.values())[0]
		else:
			print "Error: either the same color was specified for all groups or the color map '" + str(numpy.unique(groups_colors_dict.values())[0]) + "' is not valid."
			sys.exit(1)

	#create equivalency table between groups and sizes
	for g_index in range(0,groups_number):
		bb=groups_boundaries[g_index]
		tmp_beg=bb[0]
		tmp_end=bb[1]
		for tmp_size in range(tmp_beg, tmp_end+1):
			groups_sizes_dict[tmp_size]=g_index
	for tmp_size in list(set(range(1,max(groups_sizes_dict.keys())))-set(groups_sizes_dict.keys())): 	#this handles potentially unaccounted for sizes up to the maximum specified by the user
		groups_sizes_dict[tmp_size]=groups_number
	if max(groups_sizes_dict.keys())!=100000:															#this handles potentially unaccounted for sizes above the maximum specified by the user (in case it's not an open group)
		for tmp_size in range(max(groups_sizes_dict.keys())+1,100001):
			groups_sizes_dict[tmp_size]=groups_number

	#display results
	print " -found " + str(groups_number) + " cluster groups:"
	for g_index in range(0,groups_number):
		if groups_boundaries[g_index][1]==100000:
			print "   g" + str(g_index) + "=" + str(groups_boundaries[g_index][0]) + "+, " + str(groups_colors_dict[g_index])
		else:
			print "   g" + str(g_index) + "=" + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + ", " + str(groups_colors_dict[g_index])
		
# Load universe
#==============
if args.xtcfilename=="no":
	print "\nLoading file..."
	U=Universe(args.grofilename)
	all_atoms=U.selectAtoms("all")
	nb_atoms=all_atoms.numberOfAtoms()
	nb_frames_xtc=1
	nb_frames_processed=1
else:
	print "\nLoading trajectory..."
	U=Universe(args.grofilename, args.xtcfilename)
	all_atoms=U.selectAtoms("all")
	nb_atoms=all_atoms.numberOfAtoms()
	nb_frames_xtc=U.trajectory.numframes
	nb_frames_processed=0
	U.trajectory.rewind()

# Identify proteins
#==================
proteins_nb=0
proteins_sele={}
proteins_sele_string={}
proteins_sele_string_VMD={}
#case: selection file provided
#-----------------------------
if args.selection_file_prot!="auto":
	print "\nReading protein selection file..."
	with open(args.selection_file_prot) as f:
		lines = f.readlines()
	proteins_nb=len(lines)
	proteins_sele["all"]=MDAnalysis.core.AtomGroup.AtomGroup([])
	for p_index in range(0,proteins_nb):
		progress='\r -creating proteins selections: ' + str(p_index+1) + '/' + str(proteins_nb) + '        '
		sys.stdout.flush()
		sys.stdout.write(progress)
		try:
			print " p["+str(p)+"]=U.selectAtoms("+lines[p_index][0:-1]+")"
			proteins_sele[p_index]=U.selectAtoms(lines[p_index][1:-2])
			proteins_sele["all"]+=proteins_sele[p_index]
			proteins_boundaries[p_index]=[proteins_sele[p_index].indices()[0]+1,proteins_sele[p_index].indices()[proteins_sele[p_index].numberOfAtoms()]+1]
			proteins_sele_string[p_index]="bynum " + str(proteins_boundaries[p_index][0]) + ":" + str(proteins_boundaries[p_index][1])
			proteins_sele_string_VMD[p_index]="serial " + str(proteins_boundaries[p_index][0]) + " to " + str(proteins_boundaries[p_index][1])
#			if args.m_algorithm=="min":
#				proteins_sele[p_index].set_segid(str(p_index))
		except:
			print "Error: invalid selection string."
			sys.exit(1)
	proteins_nb_atoms=proteins_sele["all"].numberOfAtoms()
#case: automatic detection
#-------------------------
else:
	print "\nIdentifying proteins..."
	proteins_sele["all"] = U.selectAtoms("protein")
	proteins_nb_atoms=proteins_sele["all"].numberOfAtoms()
	if proteins_nb_atoms==0:
		print "Error: no protein detected."
		sys.exit(1)
	proteins_ca_nb={}
	proteins_ca_nmax=0
	proteins_ca_group={}
	proteins_boundaries={}

	#retrieve coord of each proteins based on non sequential residues number (works if resnum are non sequential)
	#------------------------------------------------------------------------------------------------------------
	#retrieve 1st atom info
	prec_resnum=proteins_sele["all"][0].resnum
	prec_segid=proteins_sele["all"][0].segid
	prec_atnum=proteins_sele["all"][0].number+1
	prev_atnum=proteins_sele["all"][0].number+1	#atom corresponding to the beginning of the current protein
	#browse following atoms
	for a in proteins_sele["all"][1:]:
		delta_res=a.resnum-prec_resnum
		delta_atm=a.number+1-prec_atnum
		if delta_res<0 or a.segid!=prec_segid or delta_atm>1:
			proteins_boundaries[proteins_nb]=[prev_atnum,prec_atnum]
			proteins_nb+=1
			prev_atnum=a.number+1
		prec_resnum=a.resnum    	
		prec_atnum=a.number+1
		prec_segid=a.segid		
	#add last protein section
	if prev_atnum<proteins_sele["all"][proteins_nb_atoms-1].number:
		proteins_boundaries[proteins_nb]=[prev_atnum,proteins_sele["all"][proteins_nb_atoms-1].number+1]
		proteins_nb+=1
	#display results
	print " -protein found:", proteins_nb
	print " -protein boundaries (atom numbers): see protein.sele file"
	#create protein selections and save into a txt file
	filename=os.getcwd() + '/' + str(args.output_folder) + '/proteins.sele'
	output_stat = open(filename, 'w')	
	output_stat.write("This file was generated by the script cluster_prot v" + str(version_nb) +"\n")
	output_stat.write("\n")
	output_stat.write("****************************\n")
	output_stat.write("Protein MDAnalysis selection\n")
	output_stat.write("****************************\n")
	output_stat.write("\n")	
	for p_index in range(0, proteins_nb):
		progress='\r -creating proteins selections: ' + str(p_index+1) + '/' + str(proteins_nb) + '        '
		sys.stdout.flush()
		sys.stdout.write(progress)
		proteins_sele_string[p_index]="bynum " + str(proteins_boundaries[p_index][0]) + ":" + str(proteins_boundaries[p_index][1])
		proteins_sele_string_VMD[p_index]="serial " + str(proteins_boundaries[p_index][0]) + " to " + str(proteins_boundaries[p_index][1])
		proteins_sele[p_index]=U.selectAtoms(proteins_sele_string[p_index])
#		if args.m_algorithm=="min":
#			proteins_sele[p_index].set_segid(str(p_index))
		output_stat.write("U.selectAtoms(\"" + proteins_sele_string[p_index] + "\")\n")
	output_stat.close()
	print ""

################################################################################################################################################
# FUNCTIONS: algorithm
################################################################################################################################################

def get_distances(box_dim):
	
	#method: use minimum distance between proteins
	#---------------------------------------------
	if args.m_algorithm=="min":
		dist_matrix=100000*numpy.ones((proteins_nb,proteins_nb))
		for n in range(proteins_nb,1,-1):
			p_index=proteins_nb-n
			for pp in range(p_index+1,proteins_nb):
				dist_matrix[p_index,pp]=numpy.min(MDAnalysis.analysis.distances.distance_array(numpy.float32(proteins_sele[p_index].coordinates()),numpy.float32(proteins_sele[pp].coordinates()),box_dim))	
				dist_matrix[pp,p_index]=dist_matrix[p_index,pp]

#			tmp_neighbours=proteins_sele["all"].selectAtoms("around " + str(args.cutoff_connect) + " (" + str(proteins_sele_string[p_index]) + ")")			
#			if tmp_neighbours.numberOfAtoms()>0:
#				for pp in set(tmp_neighb.segids()):
#					#debug
#					print " -", pp
#					dist_matrix[p_index,pp]=numpy.min(MDAnalysis.analysis.distances.distance_array(numpy.float32(proteins_sele[p_index].coordinates()),numpy.float32(tmp_neighbours.selectAtoms("segid " + str(pp)).coordinates()),box_dim))	
#				dist_matrix[pp,p_index]=dist_matrix[p_index,pp]
											
	#method: use distance between cog
	#--------------------------------
	else:
		#store cog coord
		tmp_proteins_coords=numpy.zeros((proteins_nb,3))
		for p_index in range(0,proteins_nb):
			tmp_proteins_coords[p_index,:]=proteins_sele[p_index].centerOfGeometry()
		#calculate distances
		dist_matrix=MDAnalysis.analysis.distances.distance_array(numpy.float32(tmp_proteins_coords), numpy.float32(tmp_proteins_coords), box_dim)

	return dist_matrix
def detect_clusters_connectivity(dist, box_dim):
	
	#use networkx algorithm
	connected=(dist<args.cutoff_connect)
	network=nx.Graph(connected)
	groups=nx.connected_components(network)
	
	return groups
def detect_clusters_density(dist, box_dim):
	
	#run DBSCAN algorithm
	dbscan_output=DBSCAN(eps=args.dbscan_dist,metric='precomputed',min_samples=args.dbscan_nb).fit(dist)

	#build 'groups' structure i.e. a list whose element are all the clusters identified
	groups=[]
	for c_lab in numpy.unique(dbscan_output.labels_):
		tmp_pos=numpy.argwhere(dbscan_output.labels_==c_lab)
		if c_lab==-1:
			for p in tmp_pos:
				groups.append([p[0]])
		else:
			g=[]
			for p in tmp_pos:
				g.append(p[0])
			groups.append(g)

	return groups

################################################################################################################################################
# FUNCTIONS: calculate statistics
################################################################################################################################################

def rolling_avg(loc_list):
	
	loc_arr=numpy.asarray(loc_list)
	shape=(loc_arr.shape[-1]-args.nb_smoothing+1,args.nb_smoothing)
	strides=(loc_arr.strides[-1],loc_arr.strides[-1])   	
	return numpy.average(numpy.lib.stride_tricks.as_strided(loc_arr, shape=shape, strides=strides), -1)
def get_sizes_sampled():
	
	global proteins_cluster_size, proteins_cluster_group, proteins_sizes_sampled, proteins_groups_sampled
	
	#sizes sampled
	#=============
	#case: gro file
	#--------------
	if args.xtcfilename=="no":
		for p_index in range(0,proteins_nb):
			proteins_sizes_sampled.append(proteins_cluster_size[p_index][0])
		proteins_sizes_sampled=list(numpy.unique(proteins_sizes_sampled))
	#case: xtc file
	#--------------
	else:
		proteins_sizes_sampled=list(numpy.unique(proteins_cluster_size[0]))		#sizes sampled by 1st protein	
		for p_index in range(1,proteins_nb):									#sizes sampled by remaining proteins
			proteins_sizes_sampled=list(numpy.unique(proteins_sizes_sampled + list(numpy.unique(proteins_cluster_size[p_index]))))

	#groups sampled
	#==============
	if args.cluster_groups_file!="no":
		#case: gro file
		#--------------
		if args.xtcfilename=="no":
			for p_index in range(0,proteins_nb):
				proteins_groups_sampled.append(proteins_cluster_group[p_index][0])
			proteins_groups_sampled=list(numpy.unique(proteins_groups_sampled))
		#case: xtc file
		#--------------
		else:
			proteins_groups_sampled=list(numpy.unique(proteins_cluster_group[0]))
			for p_index in range(1,proteins_nb):
				proteins_groups_sampled=list(numpy.unique(proteins_groups_sampled + list(numpy.unique(proteins_cluster_group[p_index]))))

	#update containers
	for c_size in proteins_sizes_sampled:
		sizes_nb[c_size]=[]
		sizes_pc[c_size]=[]
	for g_index in proteins_groups_sampled:
		groups_nb[g_index]=[]
		groups_pc[g_index]=[]

	return
def update_color_dict():
	
	global sizes_colors_nb, sizes_colors_list
	global groups_colors_nb, groups_colors_list
		
	#colormap for sizes: extract colors from jet colormap
	#-------------------
	sizes_colors_value=plt.cm.jet(numpy.linspace(0, 1, len(proteins_sizes_sampled)))
	for c_size in proteins_sizes_sampled:
		c_index=proteins_sizes_sampled.index(c_size)
		sizes_colors_dict[c_size]=sizes_colors_value[c_index]
	for k in sorted(sizes_colors_dict.iterkeys()):
		sizes_colors_list.append(sizes_colors_dict[k])
	sizes_colors_nb=numpy.size(sizes_colors_dict.keys())

	#colormap for groups
	#-------------------
	if args.cluster_groups_file!="no":
		#case: user specified color map instead of colors
		if groups_colors_map!="custom":
			tmp_cmap=cm.get_cmap(groups_colors_map)
			groups_colors_value=tmp_cmap(numpy.linspace(0, 1, groups_number))
			for g_index in range(0, groups_number):
				groups_colors_dict[g_index]=groups_colors_value[g_index]

		#create list of colours ordered by group size
		groups_colors_nb=groups_number
		for g in sorted(groups_colors_dict.iterkeys()):
			groups_colors_list.append(groups_colors_dict[g])

		#case: add the group 'other' in grey
		if groups_number in proteins_groups_sampled:
			groups_colors_nb+=1
			groups_colors_dict[groups_number]="#C0C0C0"
			groups_colors_list.append("#C0C0C0")						#choice of appending or prepending.. (top or bottom of colour bar...)
							
	return
def calc_stat():
		
	global proteins_cluster_size_mat, proteins_cluster_group_mat

	#preprocess: create dictionary of matrix of sizes sampled by each residue at each frame
	proteins_cluster_size_mat=numpy.asarray(proteins_cluster_size.values())
	proteins_cluster_group_mat=numpy.asarray(proteins_cluster_group.values())

	#case: gro file
	#==============
	if args.xtcfilename=="no":
		#preprocess: create data structure
		tmp_groups_nb={}
		tmp_groups_pc={}

		#initialise group data
		if args.cluster_groups_file!="no":
			for g_index in proteins_groups_sampled:
				tmp_groups_nb[g_index]=0
				tmp_groups_pc[g_index]=0
				
		#calculate % represented by each size / size groups
		for c_size in proteins_sizes_sampled:
			tmp_sizes=list(proteins_cluster_size_mat[:,0])
			tmp_nb=int(tmp_sizes.count(c_size)/float(c_size))
			tmp_pc=tmp_sizes.count(c_size)/float(proteins_nb)*100
			sizes_nb[c_size].append(tmp_nb)
			sizes_pc[c_size].append(tmp_pc)
			if args.cluster_groups_file!="no":
				tmp_groups_nb[groups_sizes_dict[c_size]]+=tmp_nb
				tmp_groups_pc[groups_sizes_dict[c_size]]+=tmp_pc		
	
		#size groups
		if args.cluster_groups_file!="no":
			for g_index in proteins_groups_sampled:
				groups_nb[g_index].append(tmp_groups_nb[g_index])
				groups_pc[g_index].append(tmp_groups_pc[g_index])

	#case: xtc file
	#==============
	else:
		#evolution of % represented by each size / size group
		#----------------------------------------------------
		#preprocess: create data structure
		tmp_groups_nb={}
		tmp_groups_pc={}
	
		#calculate % represented by each size / size groups
		for frame in sorted(time_stamp.iterkeys()):
			frame_index=sorted(time_stamp.keys()).index(frame)
			#initialise size groups pc				
			for g_index in proteins_groups_sampled:
				tmp_groups_nb[g_index]=0
				tmp_groups_pc[g_index]=0
			#sizes
			tmp_sizes=list(proteins_cluster_size_mat[:,frame_index])
			tmp_max_size=0
			tmp_max_nb=0
			tmp_max_pc=0
			for c_size in proteins_sizes_sampled:
				tmp_nb=int(tmp_sizes.count(c_size)/float(c_size))
				tmp_pc=tmp_sizes.count(c_size)/float(proteins_nb)*100
				sizes_nb[c_size].append(tmp_nb)
				sizes_pc[c_size].append(tmp_pc)
				if tmp_nb>0 and c_size>tmp_max_size:
					tmp_max_size=c_size
					tmp_max_nb=tmp_nb
					tmp_max_pc=tmp_pc
				if args.cluster_groups_file!="no":
					tmp_groups_nb[groups_sizes_dict[c_size]]+=tmp_nb
					tmp_groups_pc[groups_sizes_dict[c_size]]+=tmp_pc
			
			#biggest cluster size
			biggest_size[frame_index]=tmp_max_size
			biggest_nb[frame_index]=tmp_max_nb
			biggest_pc[frame_index]=tmp_max_pc		
			
			#size groups
			if args.cluster_groups_file!="no":
				for g_index in proteins_groups_sampled:
					groups_nb[g_index].append(tmp_groups_nb[g_index])
					groups_pc[g_index].append(tmp_groups_pc[g_index])
							
		#longest stability of each size group
		#------------------------------------
		if args.cluster_groups_file!="no":
			for g_index in proteins_groups_sampled:
				#find max stability of current group index for each lipid
				tmp_proteins_stability={}
				for p_index in range(0,proteins_nb):
					if g_index in proteins_cluster_group_mat[p_index,:]:
						tmp_proteins_stability[p_index]=max(len(list(v)) for g,v in itertools.groupby(proteins_cluster_group_mat[p_index,:], lambda x: x == g_index) if g)
					else:
						tmp_proteins_stability[p_index]=0
				#store maximum stability
				groups_stability[g_index]=max(tmp_proteins_stability.values())

	return
def smooth_data():

	global time_sorted, time_smoothed, biggest_size_sorted, biggest_nb_sorted, biggest_pc_sorted, biggest_size_smoothed, biggest_nb_smoothed, biggest_pc_smoothed
		
	#sort data into ordered lists
	#-----------------------------
	for frame in sorted(time_stamp.keys()):
		time_sorted.append(time_stamp[frame])
		frame_index=sorted(time_stamp.keys()).index(frame)
		biggest_size_sorted.append(biggest_size[frame_index])
		biggest_nb_sorted.append(biggest_nb[frame_index])
		biggest_pc_sorted.append(biggest_pc[frame_index])

	#calculate running average on sorted lists
	#-----------------------------------------
	if args.nb_smoothing>1:
		#time
		time_smoothed=rolling_avg(time_sorted)		
		
		#biggest cluster size
		biggest_size_smoothed=rolling_avg(biggest_size_sorted)
		biggest_nb_smoothed=rolling_avg(biggest_nb_sorted)
		biggest_pc_smoothed=rolling_avg(biggest_pc_sorted)
		
		#size
		for c_size in proteins_sizes_sampled:
			sizes_nb_smoothed[c_size]=list(rolling_avg(sizes_nb[c_size]))
			sizes_pc_smoothed[c_size]=list(rolling_avg(sizes_pc[c_size]))
		
		#groups
		if args.cluster_groups_file!="no":
			for g_index in proteins_groups_sampled:
				groups_nb_smoothed[g_index]=list(rolling_avg(groups_nb[g_index]))
				groups_pc_smoothed[g_index]=list(rolling_avg(groups_pc[g_index]))
		
	return

################################################################################################################################################
# FUNCTIONS: write outputs
################################################################################################################################################

def write_warning():
	filename_details=os.getcwd() + '/' + str(args.output_folder) + '/warning.stat'
	output_stat = open(filename_details, 'w')		
	output_stat.write("[protein clustering statistics - written by cluster_prot v" + str(version_nb) + "]\n")
	output_stat.write("\n")	
	#general info
	output_stat.write("1. Nb of proteins: " + str(proteins_nb) + "\n")
	output_stat.write("2. Cluster detection Method:\n")
	if args.m_algorithm=="min":
		output_stat.write(" - connectivity based (min distances)\n")
		output_stat.write(" - contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
	elif args.m_algorithm=="cog":
		output_stat.write(" - connectivity based (cog distances)\n")
		output_stat.write(" - contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
	else:
		output_stat.write(" - density based (DBSCAN)\n")
		output_stat.write(" - search radius = " + str(args.dbscan_dist) + " Angstrom, nb of neighbours = " + str(args.dbscan_nb) + "\n")
	output_stat.write("\n")
	#warning message
	output_stat.write("Warning: a single cluster size (" + str(proteins_sizes_sampled[0]) + ") was detected throughout the trajectory. Check the -m, -c, -r or -n options (see cluster_prot -h).")
	output_stat.close()
	
	return

#case: xtc file
#==============
#sizes
#-----
#biggest
def write_xvg_biggest():
	filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_3_biggest/xvg/1_2_clusterprot_biggest.txt'
	filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_3_biggest/xvg/1_2_clusterprot_biggest.xvg'
	output_txt = open(filename_txt, 'w')
	output_txt.write("@[lipid tail order parameters statistics - written by order_param v" + str(version_nb) + "]\n")
	output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 1_2_clusterprot_biggest.xvg.\n")
	output_xvg = open(filename_xvg, 'w')
	output_xvg.write("@ title \"Number of protein clusters\"\n")
	output_xvg.write("@ xaxis  label \"time (ns)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length 3\n")
	output_xvg.write("@ s0 legend \"size\"\n")
	output_xvg.write("@ s1 legend \"%\"\n")
	output_xvg.write("@ s2 legend \"nb\"\n")
	output_txt.write("1_2_clusterprot_biggest.xvg,1,size,k\n")
	output_txt.write("1_2_clusterprot_biggest.xvg,2,%,c\n")
	output_txt.write("1_2_clusterprot_biggest.xvg,3,nb,r\n")
	output_txt.close()
	for frame in sorted(time_stamp.iterkeys()):
		frame_index=sorted(time_stamp.keys()).index(frame)
		results=str(time_stamp[frame]) + "	" + str(biggest_size[frame_index]) + "	" + str(round(biggest_pc[frame_index],2)) + "	" + str(biggest_nb[frame_index])
		output_xvg.write(results + "\n")
	output_xvg.close()
	return
def graph_xvg_biggest():
	#create filenames
	#----------------
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_3_biggest/png/1_2_clusterprot_biggest.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_3_biggest/1_2_clusterprot_biggest.svg'

	#create figure
	#-------------
	fig=plt.figure(figsize=(8, 6.2))
	fig.suptitle("Evolution of the size of the biggest protein cluster")

			
	#plot data: %
	#------------
	ax1 = fig.add_subplot(211)
	p_upper={}
	p_upper["size"]=plt.plot(time_sorted, biggest_size_sorted, color='k', linewidth=2.0, label="size")
	fontP.set_size("small")
	ax1.legend(prop=fontP)
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('cluster size', fontsize="small")

	#plot data: nb #TO DO: make 2 y axis for the bottom bit, to include nb of clusters
	#-------------
	ax2 = fig.add_subplot(212)
	p_lower={}
	p_lower["pc"]=plt.plot(time_sorted, biggest_pc_sorted, color='c', linewidth=2.0, label="% of proteins")
	#p_lower["nb"]=plt.plot(time_sorted, biggest_nb_sorted, color='r', linewidth=2.0, label="nb of clusters")
	fontP.set_size("small")
	ax2.legend(prop=fontP)
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('% of proteins', fontsize="small")

	#save figure
	#-----------
	ax1.set_ylim(0, max(biggest_size_sorted)+2)
	ax2.set_ylim(0, 100)
	ax1.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax1.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
	ax2.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax2.yaxis.set_major_locator(MaxNLocator(nbins=5))
	plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax2.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="small" )	
	plt.subplots_adjust(top=0.9, bottom=0.07, hspace=0.37, left=0.09, right=0.96)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
	return
def write_xvg_biggest_smoothed():
	filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_3_biggest/xvg/1_3_clusterprot_biggest_smooth.txt'
	filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_3_biggest/xvg/1_3_clusterprot_biggest_smooth.xvg'
	output_txt = open(filename_txt, 'w')
	output_txt.write("@[lipid tail order parameters statistics - written by order_param v" + str(version_nb) + "]\n")
	output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 1_3_clusterprot_biggest_smooth.xvg.\n")
	output_xvg = open(filename_xvg, 'w')
	output_xvg.write("@ title \"Number of protein clusters\"\n")
	output_xvg.write("@ xaxis  label \"time (ns)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length 3\n")
	output_xvg.write("@ s0 legend \"size\"\n")
	output_xvg.write("@ s1 legend \"%\"\n")
	output_xvg.write("@ s2 legend \"nb\"\n")
	output_txt.write("1_3_clusterprot_biggest_smooth.xvg,1,size,k\n")
	output_txt.write("1_3_clusterprot_biggest_smooth.xvg,2,%,c\n")
	output_txt.write("1_3_clusterprot_biggest_smooth.xvg,3,nb,r\n")
	for frame_index in range(0, len(time_smoothed)):
		results=str(time_smoothed[frame_index]) + "	" + str(biggest_size_smoothed[frame_index]) + "	" + str(round(biggest_pc_smoothed[frame_index],2)) + "	" + str(biggest_nb_smoothed[frame_index])
		output_xvg.write(results + "\n")
	output_xvg.close()
	return
def graph_xvg_biggest_smoothed():
	#create filenames
	#----------------
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_3_biggest/png/1_3_clusterprot_biggest_smoothed.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_3_biggest/1_3_clusterprot_biggest_smoothed.svg'

	#create figure
	#-------------
	fig=plt.figure(figsize=(8, 6.2))
	fig.suptitle("Evolution of the size of the biggest protein cluster")
		
	#plot data: %
	#------------
	ax1 = fig.add_subplot(211)
	p_upper={}
	p_upper["size"]=plt.plot(time_smoothed, biggest_size_smoothed, color='k', linewidth=2.0, label="size")
	fontP.set_size("small")
	ax1.legend(prop=fontP)
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('cluster size', fontsize="small")

	#plot data: nb #TO DO: make 2 y axis for the bottom bit, to include nb of clusters
	#-------------
	ax2 = fig.add_subplot(212)
	p_lower={}
	p_lower["pc"]=plt.plot(time_smoothed, biggest_pc_smoothed, color='c', linewidth=2.0, label="% of proteins")
	#p_lower["nb"]=plt.plot(time_smoothed, biggest_nb_smoothed, color='r', linewidth=2.0, label="nb of clusters")
	fontP.set_size("small")
	ax2.legend(prop=fontP)
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('% of proteins', fontsize="small")

	#save figure
	#-----------
	ax1.set_ylim(0, numpy.floor(max(biggest_size_sorted)+2))
	ax2.set_ylim(0, 100)
	ax1.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax1.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
	ax2.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax2.yaxis.set_major_locator(MaxNLocator(nbins=5))
	plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax2.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="small" )	
	plt.subplots_adjust(top=0.9, bottom=0.07, hspace=0.37, left=0.09, right=0.96)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
	return
#each size
def write_xvg_sizes():

	filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_2_plots_1D/xvg/1_2_clusterprot_1D.txt'
	filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_2_plots_1D/xvg/1_2_clusterprot_1D.xvg'
	output_txt = open(filename_txt, 'w')
	output_txt.write("@[lipid tail order parameters statistics - written by order_param v" + str(version_nb) + "]\n")
	output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 1_2_clusterprot_1D.xvg.\n")
	output_xvg = open(filename_xvg, 'w')
	output_xvg.write("@ title \"Number of protein clusters\"\n")
	output_xvg.write("@ xaxis  label \"time (ns)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length " + str(len(proteins_groups_sampled)*2) + "\n")
	#write caption: %
	for c_index in range(0,len(proteins_sizes_sampled)):
		c_size=proteins_sizes_sampled[c_index]
		output_xvg.write("@ s" + str(c_index) + " legend \"% " + str(c_size) + "\"\n")
		output_txt.write("1_2_clusterprot_1D.xvg," + str(c_index+1) + ",% " + str(c_size) + "," + mcolors.rgb2hex(sizes_colors_dict[c_size]) + "\n")
	#write caption: nb
	for c in range(len(proteins_sizes_sampled),len(proteins_sizes_sampled)*2):
		c_size=proteins_sizes_sampled[c-len(proteins_sizes_sampled)]
		output_xvg.write("@ s" + str(c) + " legend \"nb " + str(c_size) + "\"\n")
		output_txt.write("1_2_clusterprot_1D.xvg," + str(c+1) + ",nb " + str(c_size) + "," + mcolors.rgb2hex(sizes_colors_dict[c_size]) + "\n")
	output_txt.close()
	#write results
	for frame in sorted(time_stamp.iterkeys()):
		results=str(time_stamp[frame])
		frame_index=sorted(time_stamp.keys()).index(frame)
		for c_size in proteins_sizes_sampled:
			results+="	" + str(round(sizes_pc[c_size][frame_index],2))
		for c_size in proteins_sizes_sampled:
			results+="	" + str(round(sizes_nb[c_size][frame_index],2))
		output_xvg.write(results + "\n")
	output_xvg.close()
	
	return
def graph_xvg_sizes():
	
	#create filenames
	#----------------
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_2_plots_1D/png/1_2_clusterprot_1D.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_2_plots_1D/1_2_clusterprot_1D.svg'

	#create figure
	#-------------
	fig=plt.figure(figsize=(8, 6.2))
	fig.suptitle("Evolution of proteins distribution")
		
	#plot data: %
	#------------
	ax1 = fig.add_subplot(211)
	p_upper={}
	for c_size in proteins_sizes_sampled:
		tmp_label=str(c_size)
		p_upper[c_size]=plt.plot(time_sorted, sizes_pc[c_size], color=sizes_colors_dict[c_size], linewidth=2.0, label=tmp_label)
	fontP.set_size("small")
	ax1.legend(prop=fontP)
	plt.title("%", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('% of proteins', fontsize="small")

	#plot data: nb
	#-------------
	ax2 = fig.add_subplot(212)
	p_lower={}
	for c_size in proteins_sizes_sampled:
		tmp_label=str(c_size)
		p_lower[c_size]=plt.plot(time_sorted, sizes_nb[c_size], color=sizes_colors_dict[c_size], linewidth=2.0, label=tmp_label)
	fontP.set_size("small")
	ax2.legend(prop=fontP)
	plt.title("nb", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('nb of clusters', fontsize="small")

	#save figure
	#-----------
	ax1.set_ylim(0, 100)
	ax2.set_ylim(0, max(max(sizes_nb.values()))+1)
	ax1.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax1.yaxis.set_major_locator(MaxNLocator(nbins=5))
	ax2.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax2.yaxis.set_major_locator(MaxNLocator(nbins=5))
	plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax2.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="small" )	
	plt.subplots_adjust(top=0.9, bottom=0.07, hspace=0.37, left=0.09, right=0.96)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()

	return
def write_xvg_sizes_smoothed():
	filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_4_plots_1D_smoothed/xvg/1_4_clusterprot_1D_smoothed.txt'
	output_txt = open(filename_txt, 'w')
	output_txt.write("@[lipid tail order parameters statistics - written by order_param v" + str(version_nb) + "]\n")
	output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 1_4_clusterprot_1D_smoothed.xvg.\n")
	filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_4_plots_1D_smoothed/xvg/1_4_clusterprot_1D_smoothed.xvg'
	output_xvg = open(filename_xvg, 'w')
	output_xvg.write("@ title \"Number of protein clusters\"\n")
	output_xvg.write("@ xaxis  label \"time (ns)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length " + str(len(proteins_sizes_sampled)*2) + "\n")
	#write caption: %
	for c_index in range(0,len(proteins_sizes_sampled)):
		c_size=proteins_sizes_sampled[c_index]
		output_xvg.write("@ s" + str(c_index) + " legend \"% " + str(c_size) + "\"\n")
		output_txt.write("1_4_clusterprot_1D_smoothed.xvg," + str(c_index+1) + ",% " + str(c_size) + "," + mcolors.rgb2hex(sizes_colors_dict[c_size]) + "\n")
	#write caption: nb
	for c in range(len(proteins_sizes_sampled),len(proteins_sizes_sampled)*2):
		c_size=proteins_sizes_sampled[c-len(proteins_sizes_sampled)]
		output_xvg.write("@ s" + str(c) + " legend \"nb " + str(c_size) + "\"\n")
		output_txt.write("1_4_clusterprot_1D_smoothed.xvg," + str(c+1) + ",nb " + str(c_size) + "," + mcolors.rgb2hex(sizes_colors_dict[c_size]) + "\n")
	output_txt.close()
	#write results
	for frame_index in range(0, len(time_smoothed)):
		results=str(time_smoothed[frame_index])
		for c_size in proteins_sizes_sampled:
			results+="	" + str(round(sizes_pc_smoothed[c_size][frame_index],2))
		for c_size in proteins_sizes_sampled:
			results+="	" + str(round(sizes_nb_smoothed[c_size][frame_index],2))
		output_xvg.write(results + "\n")
	output_xvg.close()

	return
def graph_xvg_sizes_smoothed():
	
	#create filenames
	#----------------
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_4_plots_1D_smoothed/png/1_4_clusterprot_1D_smoothed.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_4_plots_1D_smoothed/1_4_clusterprot_1D_smoothed.svg'

	#create figure
	#-------------
	fig=plt.figure(figsize=(8, 6.2))
	fig.suptitle("Evolution of proteins distribution")
		
	#plot data: %
	#------------
	ax1 = fig.add_subplot(211)
	p_upper={}
	for c_size in proteins_sizes_sampled:
		tmp_label=str(c_size)
		p_upper[c_size]=plt.plot(time_smoothed, sizes_pc_smoothed[c_size], color=sizes_colors_dict[c_size], linewidth=2.0, label=tmp_label)
	fontP.set_size("small")
	ax1.legend(prop=fontP)
	plt.title("%", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('% of proteins', fontsize="small")

	#plot data: nb
	#-------------
	ax2 = fig.add_subplot(212)
	p_lower={}
	for c_size in proteins_sizes_sampled:
		tmp_label=str(c_size)
		p_lower[c_size]=plt.plot(time_smoothed, sizes_nb_smoothed[c_size], color=sizes_colors_dict[c_size], linewidth=2.0, label=tmp_label)
	fontP.set_size("small")
	ax2.legend(prop=fontP)
	plt.title("nb", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('nb of clusters', fontsize="small")

	#save figure
	#-----------
	ax1.set_ylim(0, 100)
	ax2.set_ylim(0, max(max(sizes_nb_smoothed.values()))+1)
	ax1.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax1.yaxis.set_major_locator(MaxNLocator(nbins=5))
	ax2.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax2.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
	plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax2.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="small" )	
	plt.subplots_adjust(top=0.9, bottom=0.07, hspace=0.37, left=0.09, right=0.96)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()

	return
#2D summary
def graph_aggregation_2D_sizes():
		
	#create filenames
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_1_plots_2D//png/1_1_clusterprot_2D.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_1_plots_2D//1_1_clusterprot_2D.svg'

	#create data
	lip_2D_evolution=numpy.zeros((proteins_nb,len(time_stamp.keys())))
	for p_index in range(0,proteins_nb):
		lip_2D_evolution[p_index,:]=numpy.asarray(proteins_cluster_size[p_index])

	#build color map
	color_map=mcolors.LinearSegmentedColormap.from_list('custom', sizes_colors_list, sizes_colors_nb)

	#determine nb of colours and their boundaries
	bounds=[]
	cb_ticks_lab=[]
	for c in sorted(sizes_colors_dict.iterkeys()):
		bounds.append(c-0.5)
		cb_ticks_lab.append(str(c))
	bounds.append(sorted(sizes_colors_dict.iterkeys())[-1]+0.5)
	norm=mpl.colors.BoundaryNorm(bounds, color_map.N)
				
	#create figure ('norm' requires at least 2 elements to work)
	fig=plt.figure(figsize=(9, 8))
	ax_plot=fig.add_axes([0.10, 0.1, 0.75, 0.77])	
	ax_plot.matshow(lip_2D_evolution, origin='lower', interpolation='nearest', cmap=color_map, aspect='auto', norm=norm)

	#create color bar
	ax_cbar=fig.add_axes([0.88, 0.1, 0.025, 0.77])
	cb=mpl.colorbar.ColorbarBase(ax_cbar, cmap=color_map, norm=norm, boundaries=bounds)

	#position and label color bar ticks
	cb_ticks_pos=[]
	for b in range(1,len(bounds)):
		cb_ticks_pos.append(bounds[b-1]+(bounds[b]-bounds[b-1])/2)
	cb_ticks_pos.append(bounds[-1])
	cb.set_ticks(cb_ticks_pos)
	cb.set_ticklabels(cb_ticks_lab)
	for t in cb.ax.get_yticklabels():
		t.set_fontsize('xx-small')
	
	#x axis ticks
	ax_plot.xaxis.set_label_position('bottom') 
	ax_plot.xaxis.set_ticks_position('bottom')
	xticks_pos=ax_plot.xaxis.get_ticklocs()[1:-1]
	tmp_xticks_lab=[""]
	for t in sorted(time_stamp.values()):
		tmp_xticks_lab.append('{0:0g}'.format(numpy.floor(t)))
	xticks_lab=[""]
	for t in xticks_pos:
		xticks_lab.append(tmp_xticks_lab[int(t)+1])
	ax_plot.xaxis.set_ticklabels(xticks_lab)

	#y axis ticks (increase the index by 1 to get 1-based numbers)
	ax_plot.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x,p: '{0:0g}'.format(x+1)))
	
	#set title and limits
	ax_plot.set_xlabel("time (ns)", fontsize="medium")
	ax_plot.set_ylabel("protein #", fontsize="medium")
	plt.setp(ax_plot.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax_plot.yaxis.get_majorticklabels(), fontsize="small" )
	ax_plot.yaxis.set_major_locator(MaxNLocator(prune='lower'))	
	ax_plot.set_title("Evolution of the cluster size in which proteins are involved", fontsize="medium")	
	ax_cbar.set_ylabel('cluster size',fontsize='small')
	
	#save figure
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
			
	return

#groups
#------
#each group
def write_xvg_groups():
	filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_2_plots_1D/xvg/2_2_clusterprot_1D.txt'
	filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_2_plots_1D/xvg/2_2_clusterprot_1D.xvg'
	output_txt = open(filename_txt, 'w')
	output_txt.write("@[lipid tail order parameters statistics - written by order_param v" + str(version_nb) + "]\n")
	output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 2_2_clusterprot_1D.xvg.\n")
	output_xvg = open(filename_xvg, 'w')
	output_xvg.write("@ title \"Number of protein clusters\"\n")
	output_xvg.write("@ xaxis  label \"time (ns)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length " + str(len(proteins_groups_sampled)*2) + "\n")
	#write caption: %
	for g in range(0,len(proteins_groups_sampled)):
		g_index=proteins_groups_sampled[g]
		if g_index==groups_number:
			output_xvg.write("@ s" + str(g) + " legend \"% other\"\n")
			output_txt.write("2_2_clusterprot_1D.xvg," + str(g+1) + ",% other," + mcolors.rgb2hex(groups_colors_dict[g_index]) + "\n")
		elif groups_boundaries[g_index][1]==100000:
			output_xvg.write("@ s" + str(g) + " legend \"% >=" + str(groups_boundaries[g_index][0]) + "\"\n")
			output_txt.write("2_2_clusterprot_1D.xvg," + str(g+1) + ",% >=" + str(groups_boundaries[g_index][0]) + "," + mcolors.rgb2hex(groups_colors_dict[g_index]) + "\n")
		else:
			output_xvg.write("@ s" + str(g) + " legend \"%" + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + "\"\n")
			output_txt.write("2_2_clusterprot_1D.xvg," + str(g+1) + ",% " + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + "," + mcolors.rgb2hex(groups_colors_dict[g_index]) + "\n")
	#write caption: nb
	for g in range(len(proteins_groups_sampled),len(proteins_groups_sampled)*2):
		g_index=proteins_groups_sampled[g-len(proteins_groups_sampled)]
		if g_index==groups_number:
			output_xvg.write("@ s" + str(g) + " legend \"nb other\"\n")
			output_txt.write("2_2_clusterprot_1D.xvg," + str(g+1) + ",nb other," + mcolors.rgb2hex(groups_colors_dict[g_index]) + "\n")
		elif groups_boundaries[g_index][1]==100000:
			output_xvg.write("@ s" + str(g) + " legend \"nb >=" + str(groups_boundaries[g_index][0]) + "\"\n")
			output_txt.write("2_2_clusterprot_1D.xvg," + str(g+1) + ",nb >=" + str(groups_boundaries[g_index][0]) + "," + mcolors.rgb2hex(groups_colors_dict[g_index]) + "\n")
		else:
			output_xvg.write("@ s" + str(g) + " legend \"nb" + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + "\"\n")
			output_txt.write("2_2_clusterprot_1D.xvg," + str(g+1) + ",nb " + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + "," + mcolors.rgb2hex(groups_colors_dict[g_index]) + "\n")
	output_txt.close()
	#write results
	for frame in sorted(time_stamp.iterkeys()):
		results=str(time_stamp[frame])
		frame_index=sorted(time_stamp.keys()).index(frame)
		for g_index in proteins_groups_sampled:
			results+="	" + str(round(groups_pc[g_index][frame_index],2))
		for g_index in proteins_groups_sampled:
			results+="	" + str(round(groups_nb[g_index][frame_index],2))
		output_xvg.write(results + "\n")
	output_xvg.close()
	
	return
def graph_xvg_groups():
	
	#create filenames
	#----------------
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_2_plots_1D/png/2_2_clusterprot_1D.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_2_plots_1D/2_2_clusterprot_1D.svg'

	#create figure
	#-------------
	fig=plt.figure(figsize=(8, 6.2))
	fig.suptitle("Evolution of proteins distribution")
		
	#plot data: %
	#------------
	ax1 = fig.add_subplot(211)
	p_upper={}
	for g_index in proteins_groups_sampled:
		if g_index==groups_number:
			tmp_label="other"
		elif groups_boundaries[g_index][1]==100000:
			tmp_label=">=" + str(groups_boundaries[g_index][0])
		else:
			tmp_label=str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1])
		p_upper[g_index]=plt.plot(time_sorted, groups_pc[g_index], color=groups_colors_dict[g_index], linewidth=2.0, label=tmp_label)
	fontP.set_size("small")
	ax1.legend(prop=fontP)
	plt.title("%", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('% of proteins', fontsize="small")

	#plot data: nb
	#-------------
	ax2 = fig.add_subplot(212)
	p_lower={}
	for g_index in proteins_groups_sampled:
		if g_index==groups_number:
			tmp_label="other"
		elif groups_boundaries[g_index][1]==100000:
			tmp_label=">=" + str(groups_boundaries[g_index][0])
		else:
			tmp_label=str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1])
		p_lower[g_index]=plt.plot(time_sorted, groups_nb[g_index], color=groups_colors_dict[g_index], linewidth=2.0, label=tmp_label)
	fontP.set_size("small")
	ax2.legend(prop=fontP)
	plt.title("nb", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('nb of clusters', fontsize="small")

	#save figure
	#-----------
	ax1.set_ylim(0, 100)
	ax2.set_ylim(0, max(max(groups_nb.values()))+1)
	ax1.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax1.yaxis.set_major_locator(MaxNLocator(nbins=5))
	ax2.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax2.yaxis.set_major_locator(MaxNLocator(nbins=5))
	plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax2.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="small" )	
	plt.subplots_adjust(top=0.9, bottom=0.07, hspace=0.37, left=0.09, right=0.96)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()

	return
def write_xvg_groups_smoothed():
	filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_3_plots_1D_smoothed/xvg/2_3_clusterprot_1D_smoothed.txt'
	output_txt = open(filename_txt, 'w')
	output_txt.write("@[lipid tail order parameters statistics - written by order_param v" + str(version_nb) + "]\n")
	output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 2_3_clusterprot_1D_smoothed.xvg.\n")
	filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_3_plots_1D_smoothed/xvg/2_3_clusterprot_1D_smoothed.xvg'
	output_xvg = open(filename_xvg, 'w')
	output_xvg.write("@ title \"Number of protein clusters\"\n")
	output_xvg.write("@ xaxis  label \"time (ns)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length " + str(len(proteins_groups_sampled)*2) + "\n")
	#write caption: %
	for g in range(0,len(proteins_groups_sampled)):
		g_index=proteins_groups_sampled[g]
		if g_index==groups_number:
			output_xvg.write("@ s" + str(g) + " legend \"% other\"\n")
			output_txt.write("2_3_clusterprot_1D_smoothed.xvg," + str(g+1) + ",% other," + mcolors.rgb2hex(groups_colors_dict[g_index]) + "\n")
		elif groups_boundaries[g_index][1]==100000:
			output_xvg.write("@ s" + str(g) + " legend \"% >=" + str(groups_boundaries[g_index][0]) + "\"\n")
			output_txt.write("2_3_clusterprot_1D_smoothed.xvg," + str(g+1) + ",% >=" + str(groups_boundaries[g_index][0]) + "," + mcolors.rgb2hex(groups_colors_dict[g_index]) + "\n")
		else:
			output_xvg.write("@ s" + str(g) + " legend \"%" + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + "\"\n")
			output_txt.write("2_3_clusterprot_1D_smoothed.xvg," + str(g+1) + ",% " + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + "," + mcolors.rgb2hex(groups_colors_dict[g_index]) + "\n")
	#write caption: nb
	for g in range(len(proteins_groups_sampled),len(proteins_groups_sampled)*2):
		g_index=proteins_groups_sampled[g-len(proteins_groups_sampled)]
		if g_index==groups_number:
			output_xvg.write("@ s" + str(g) + " legend \"nb other\"\n")
			output_txt.write("2_3_clusterprot_1D_smoothed.xvg," + str(g+1) + ",nb other," + mcolors.rgb2hex(groups_colors_dict[g_index]) + "\n")
		elif groups_boundaries[g_index][1]==100000:
			output_xvg.write("@ s" + str(g) + " legend \"nb >=" + str(groups_boundaries[g_index][0]) + "\"\n")
			output_txt.write("2_3_clusterprot_1D_smoothed.xvg," + str(g+1) + ",nb >=" + str(groups_boundaries[g_index][0]) + "," + mcolors.rgb2hex(groups_colors_dict[g_index]) + "\n")
		else:
			output_xvg.write("@ s" + str(g) + " legend \"nb" + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + "\"\n")
			output_txt.write("2_3_clusterprot_1D_smoothed.xvg," + str(g+1) + ",nb " + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + "," + mcolors.rgb2hex(groups_colors_dict[g_index]) + "\n")
	output_txt.close()
	#write results
	for frame_index in range(0, len(time_smoothed)):
		results=str(time_smoothed[frame_index])
		for g_index in proteins_groups_sampled:
			results+="	" + str(round(groups_pc_smoothed[g_index][frame_index],2))
		for g_index in proteins_groups_sampled:
			results+="	" + str(round(groups_nb_smoothed[g_index][frame_index],2))
		output_xvg.write(results + "\n")
	output_xvg.close()

	return
def graph_xvg_groups_smoothed():
	
	#create filenames
	#----------------
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_3_plots_1D_smoothed/png/2_3_clusterprot_1D_smoothed.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_3_plots_1D_smoothed/2_3_clusterprot_1D_smoothed.svg'

	#create figure
	#-------------
	fig=plt.figure(figsize=(8, 6.2))
	fig.suptitle("Evolution of proteins distribution")
		
	#plot data: %
	#------------
	ax1 = fig.add_subplot(211)
	p_upper={}
	for g_index in proteins_groups_sampled:
		if g_index==groups_number:
			tmp_label="other"
		elif groups_boundaries[g_index][1]==100000:
			tmp_label=">=" + str(groups_boundaries[g_index][0])
		else:
			tmp_label=str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1])
		p_upper[g_index]=plt.plot(time_smoothed, groups_pc_smoothed[g_index], color=groups_colors_dict[g_index], linewidth=2.0, label=tmp_label)
	fontP.set_size("small")
	ax1.legend(prop=fontP)
	plt.title("%", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('% of proteins', fontsize="small")

	#plot data: nb
	#-------------
	ax2 = fig.add_subplot(212)
	p_lower={}
	for g_index in proteins_groups_sampled:	
		if g_index==groups_number:
			tmp_label="other"
		elif groups_boundaries[g_index][1]==100000:
			tmp_label=">=" + str(groups_boundaries[g_index][0])
		else:
			tmp_label=str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1])
		p_lower[g_index]=plt.plot(time_smoothed, groups_nb_smoothed[g_index], color=groups_colors_dict[g_index], linewidth=2.0, label=tmp_label)
	fontP.set_size("small")
	ax2.legend(prop=fontP)
	plt.title("nb", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('nb of clusters', fontsize="small")

	#save figure
	#-----------
	ax1.set_ylim(0, 100)
	ax2.set_ylim(0, max(max(groups_nb_smoothed.values()))+1)
	ax1.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax1.yaxis.set_major_locator(MaxNLocator(nbins=5))
	ax2.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax2.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
	plt.setp(ax1.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax1.yaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax2.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax2.yaxis.get_majorticklabels(), fontsize="small" )	
	plt.subplots_adjust(top=0.9, bottom=0.07, hspace=0.37, left=0.09, right=0.96)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()

	return
#2D summary
def graph_aggregation_2D_groups():
	
	#create filenames
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_1_plots_2D//png/2_1_clusterprot_2D.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_1_plots_2D//2_1_clusterprot_2D.svg'

	#create data
	lip_2D_evolution=numpy.zeros((proteins_nb,len(time_stamp.keys())))
	for p_index in range(0,proteins_nb):
		lip_2D_evolution[p_index,:]=numpy.asarray(proteins_cluster_group[p_index])

	#build color map
	color_map=mcolors.LinearSegmentedColormap.from_list('custom', groups_colors_list, groups_colors_nb)
	
	#determine nb of colours and their boundaries
	bounds=[]
	cb_ticks_lab=[]
	for g_index in sorted(groups_colors_dict.iterkeys()):
		bounds.append(g_index-0.5)
		if g_index==groups_number:
			cb_ticks_lab.append("other")
		elif groups_boundaries[g_index][1]==100000:
			cb_ticks_lab.append(">=" + str(groups_boundaries[g_index][0]))
		else:
			cb_ticks_lab.append(str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]))
	bounds.append(sorted(groups_colors_dict.iterkeys())[-1]+0.5)
	norm=mpl.colors.BoundaryNorm(bounds, color_map.N)
				
	#create figure ('norm' requires at least 2 elements to work)
	fig=plt.figure(figsize=(9, 8))
	ax_plot=fig.add_axes([0.09, 0.1, 0.75, 0.77])	
	ax_plot.matshow(lip_2D_evolution, origin='lower', interpolation='nearest', cmap=color_map, aspect='auto', norm=norm)

	#create color bar
	ax_cbar=fig.add_axes([0.87, 0.1, 0.025, 0.77])
	cb=mpl.colorbar.ColorbarBase(ax_cbar, cmap=color_map, norm=norm, boundaries=bounds)

	#position and label color bar ticks
	cb_ticks_pos=[]
	for b in range(1,len(bounds)):
		cb_ticks_pos.append(bounds[b-1]+(bounds[b]-bounds[b-1])/2)
	cb_ticks_pos.append(bounds[-1])
	cb.set_ticks(cb_ticks_pos)
	cb.set_ticklabels(cb_ticks_lab)
	for t in cb.ax.get_yticklabels():
		t.set_fontsize('small')
				
	#x axis ticks
	ax_plot.xaxis.set_label_position('bottom') 
	ax_plot.xaxis.set_ticks_position('bottom')
	xticks_pos=ax_plot.xaxis.get_ticklocs()[1:-1]
	tmp_xticks_lab=[""]
	for t in sorted(time_stamp.values()):
		tmp_xticks_lab.append('{0:0g}'.format(numpy.floor(t)))
	xticks_lab=[""]
	for t in xticks_pos:
		xticks_lab.append(tmp_xticks_lab[int(t)+1])
	ax_plot.xaxis.set_ticklabels(xticks_lab)

	#y axis ticks (increase the index by 1 to get 1-based numbers)
	ax_plot.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x,p: '{0:0g}'.format(x+1)))
	
	#set title and limits
	ax_plot.set_xlabel("time (ns)", fontsize="medium")
	ax_plot.set_ylabel("protein #", fontsize="medium")
	plt.setp(ax_plot.xaxis.get_majorticklabels(), fontsize="small" )
	plt.setp(ax_plot.yaxis.get_majorticklabels(), fontsize="small" )
	ax_plot.yaxis.set_major_locator(MaxNLocator(prune='lower'))	
	ax_plot.set_title("Evolution of the cluster size group in which proteins are involved", fontsize="medium")	
	ax_cbar.set_ylabel('size groups',fontsize='medium')
	
	#save figure
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
			
	return
#stability
def write_stability_groups():
	
	filename_details=os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_0_clusterprot_stability.stat'
	output_stat = open(filename_details, 'w')		
	output_stat.write("[protein clustering statistics - written by cluster_prot v" + str(version_nb) + "]\n")
	output_stat.write("\n")

	#general info
	output_stat.write("1. Nb of proteins: " + str(proteins_nb) + "\n")
	output_stat.write("2. Cluster detection Method:\n")
	if args.m_algorithm=="min":
		output_stat.write(" - connectivity based (min distances)\n")
		output_stat.write(" - contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
	elif args.m_algorithm=="cog":
		output_stat.write(" - connectivity based (cog distances)\n")
		output_stat.write(" - contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
	else:
		output_stat.write(" - density based (DBSCAN)\n")
		output_stat.write(" - search radius = " + str(args.dbscan_dist) + " Angstrom, nb of neighbours = " + str(args.dbscan_nb) + "\n")
	output_stat.write("\n")
	output_stat.write("Maximum stability (in number of consecutive frames) for each cluster groups\n")
	output_stat.write("Note: frames skipped is not taken into account (the nb below correspnd to consecutive frames processed)\n")
	
	tmp_cap1=""
	tmp_cap2="-----"
	for g_index in proteins_groups_sampled:
		if g_index==groups_number:
			tmp_cap1+="	other"
		elif groups_boundaries[g_index][1]==100000:
			tmp_cap1+="	>=" + str(groups_boundaries[g_index][0])
		else:
			tmp_cap1+="	" + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1])
		tmp_cap2+="--------"

	output_stat.write("\n")
	output_stat.write(tmp_cap1 + "\n")
	output_stat.write(tmp_cap2 + "\n")
	results=str("")
	for g_index in proteins_groups_sampled:
		results+= "	" + str(groups_stability[g_index])
	output_stat.write(results + "\n")

	output_stat.close()
	
	return

#annotations
#===========
def write_frame_stat(f_nb, f_index, t):

	#case: gro file or xtc summary
	#=============================
	if f_index=="all" and t=="all":
		#sizes
		#-----
		#create file
		if args.xtcfilename=="no":
			filename_details=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/' + args.grofilename[:-4] + '_annotated_clustprot_sizes_sampled.stat'		
		else:
			filename_details=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/' + args.xtcfilename[:-4] + '_annotated_clustprot_sizes_sampled.stat'		
		output_stat = open(filename_details, 'w')		
	
		#general info
		output_stat.write("[protein clustering statistics - written by cluster_prot v" + str(version_nb) + "]\n")
		output_stat.write("\n")
		output_stat.write("1. nb of proteins: " + str(proteins_nb) + "\n")
		output_stat.write("\n")
		output_stat.write("2. cluster detection Method:\n")
		if args.m_algorithm=="min":
			output_stat.write(" - connectivity based (min distances)\n")
			output_stat.write(" - contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
		elif args.m_algorithm=="cog":
			output_stat.write(" - connectivity based (cog distances)\n")
			output_stat.write(" - contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
		else:
			output_stat.write(" - density based (DBSCAN)\n")
			output_stat.write(" - search radius = " + str(args.dbscan_dist) + " Angstrom, nb of neighbours = " + str(args.dbscan_nb) + "\n")
		if args.xtcfilename!="no":
			output_stat.write("\n")
			output_stat.write("3. nb frames processed:	" + str(nb_frames_processed) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
		
		#data results
		output_stat.write("\n")
		output_stat.write("Sizes sampled: ")
		tmp_res=str(proteins_sizes_sampled[0])
		for c_size in proteins_sizes_sampled[1:]:
			tmp_res+="," + str(c_size)
		output_stat.write(tmp_res + "\n")
		output_stat.write("\n")
		output_stat.close()
		
		#groups
		#------
		if args.cluster_groups_file!="no":
			#create file
			if args.xtcfilename=="no":
				filename_details=os.getcwd() + '/' + str(args.output_folder) + '/2_groups/' + args.grofilename[:-4] + '_annotated_clustprot_groups_sampled.stat'		
			else:
				filename_details=os.getcwd() + '/' + str(args.output_folder) + '/2_groups/' + args.xtcfilename[:-4] + '_annotated_clustprot_groups_sampled.stat'		
			output_stat = open(filename_details, 'w')		

			#general info
			output_stat.write("[protein clustering statistics - written by cluster_prot v" + str(version_nb) + "]\n")
			output_stat.write("\n")
			output_stat.write("1. nb of proteins: " + str(proteins_nb) + "\n")
			output_stat.write("\n")
			output_stat.write("2. cluster detection Method:\n")
			if args.m_algorithm=="min":
				output_stat.write(" - connectivity based (min distances)\n")
				output_stat.write(" - contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
			elif args.m_algorithm=="cog":
				output_stat.write(" - connectivity based (cog distances)\n")
				output_stat.write(" - contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
			else:
				output_stat.write(" - density based (DBSCAN)\n")
				output_stat.write(" - search radius = " + str(args.dbscan_dist) + " Angstrom, nb of neighbours = " + str(args.dbscan_nb) + "\n")
			if args.xtcfilename!="no":
				output_stat.write("\n")
				output_stat.write("3. nb frames processed:	" + str(nb_frames_processed) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")

			#group index definition
			output_stat.write("\n")
			output_stat.write("Group size ranges:\n")
			for g_index in range(0,groups_number):
				if groups_boundaries[g_index][1]==100000:
					output_stat.write(str(g_index) + "=" + str(groups_boundaries[g_index][0]) + "+\n")
				else:
					output_stat.write(str(g_index) + "=" + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + "\n")

			#data results
			output_stat.write("\n")
			output_stat.write("Groups sampled: ")
			tmp_res=str(proteins_groups_sampled[0])
			for g_index in proteins_groups_sampled[1:]:
				tmp_res+="," + str(g_index)
			output_stat.write(tmp_res + "\n")		
			output_stat.write("\n")
		output_stat.close()
		
	#case: xtc snapshot
	#==================
	else:
		#sizes
		#-----
		#create file
		if args.xtcfilename=="no":
			filename_details=os.getcwd() + '/' + str(args.output_folder) + '/3_snapshots/sizes/' + args.grofilename[:-4] + '_annotated_clusterprot_sizes.stat'
		else:
			filename_details=os.getcwd() + "/" + str(args.output_folder) + '/3_snapshots/sizes/' + args.xtcfilename[:-4] + '_annotated_clusterprot_sizes_' + str(int(t)).zfill(5) + 'ns.stat'
		output_stat = open(filename_details, 'w')		
	
		#general info
		output_stat.write("[protein clustering statistics - written by cluster_prot v" + str(version_nb) + "]\n")
		output_stat.write("\n")
		output_stat.write("1. bb of proteins: " + str(proteins_nb) + "\n")
		output_stat.write("2. cluster detection Method:\n")
		if args.m_algorithm=="min":
			output_stat.write(" - connectivity based (min distances)\n")
			output_stat.write(" - contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
		elif args.m_algorithm=="cog":
			output_stat.write(" - connectivity based (cog distances)\n")
			output_stat.write(" - contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
		else:
			output_stat.write(" - density based (DBSCAN)\n")
			output_stat.write(" - search radius = " + str(args.dbscan_dist) + " Angstrom, nb of neighbours = " + str(args.dbscan_nb) + "\n")
		if args.xtcfilename!="no":
			output_stat.write("\n")
			output_stat.write("3. time: " + str(t) + "ns (frame " + str(f_nb) + "/" + str(nb_frames_xtc) + ")\n")

		#what's in this file
		output_stat.write("\n")
		output_stat.write("Distribution of proteins by cluster size:\n")
		output_stat.write("\n")
		tmp_cap1="cluster_size"
		tmp_cap2="-----------"
		tmp_res="% proteins"
		for c_size in proteins_sizes_sampled:
			tmp_cap1+="	" + str(c_size)
			tmp_cap2+="--------"
		output_stat.write(tmp_cap1+"\n")
		output_stat.write(tmp_cap2+"\n")
		for c_size in proteins_sizes_sampled:		
			tmp_res+="	" + str(round(sizes_pc[c_size][f_index],1))
		output_stat.write(tmp_res + "\n")		
		output_stat.close()
		
		#groups
		#======
		if args.cluster_groups_file!="no":
			#create file
			if args.xtcfilename=="no":
				filename_details=os.getcwd() + '/' + str(args.output_folder) + '/3_snapshots/groups/' + args.grofilename[:-4] + '_annotated_clusterprot_groups.stat'
			else:
				filename_details=os.getcwd() + "/" + str(args.output_folder) + '/3_snapshots/groups/' + args.xtcfilename[:-4] + '_annotated_clusterprot_groups_' + str(int(t)).zfill(5) + 'ns.stat'
			output_stat = open(filename_details, 'w')		
		
			#general info
			output_stat.write("[protein clustering statistics - written by cluster_prot v" + str(version_nb) + "]\n")
			output_stat.write("\n")
			output_stat.write("1. nb of proteins: " + str(proteins_nb) + "\n")
			output_stat.write("2. cluster detection Method:\n")
			if args.m_algorithm=="min":
				output_stat.write(" - connectivity based (min distances)\n")
				output_stat.write(" - contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
			elif args.m_algorithm=="cog":
				output_stat.write(" - connectivity based (cog distances)\n")
				output_stat.write(" - contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
			else:
				output_stat.write(" - density based (DBSCAN)\n")
				output_stat.write(" - search radius = " + str(args.dbscan_dist) + " Angstrom, nb of neighbours = " + str(args.dbscan_nb) + "\n")
			if args.xtcfilename!="no":
				output_stat.write("\n")
				output_stat.write("3. time: " + str(t) + "ns (frame " + str(f_nb) + "/" + str(nb_frames_xtc) + ")\n")
		
			#group index definition
			output_stat.write("\n")
			output_stat.write("Group size ranges:\n")
			for g_index in range(0,groups_number):
				if groups_boundaries[g_index][1]==100000:
					output_stat.write(str(g_index) + "=" + str(groups_boundaries[g_index][0]) + "+\n")
				else:
					output_stat.write(str(g_index) + "=" + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + "\n")

			#what's in this file
			output_stat.write("\n")
			output_stat.write("Distribution of proteins by cluster size group:\n")
			output_stat.write("\n")
			tmp_cap1="cluster_group"
			tmp_cap2="-----------"
			tmp_res="% proteins"
			for c_group in proteins_groups_sampled:
				tmp_cap1+="	" + str(c_group)
				tmp_cap2+="--------"
			output_stat.write(tmp_cap1+"\n")
			output_stat.write(tmp_cap2+"\n")
			for c_group in proteins_groups_sampled:		
				tmp_res+="	" + str(round(groups_pc[c_group][f_index],1))
			output_stat.write(tmp_res + "\n")		
			output_stat.close()

	return
def write_frame_snapshot(f_index, t):
	
	#sizes
	#=====
	#store cluster info in beta factor field
	for p_index in range(0,proteins_nb):
			proteins_sele[p_index].set_bfactor(proteins_cluster_size[p_index][f_index])
	
	#write annotated file
	if args.xtcfilename=="no":
		all_atoms.write(os.getcwd() + '/' + str(args.output_folder) + '/3_snapshots/sizes/' + args.grofilename[:-4] + '_annotated_clusterprot_sizes', format="PDB")
	else:
		tmp_name=os.getcwd() + "/" + str(args.output_folder) + '/3_snapshots/sizes/' + args.xtcfilename[:-4] + '_annotated_clusterprot_sizes_' + str(int(t)).zfill(5) + 'ns.pdb'
		W=Writer(tmp_name, nb_atoms)
		W.write(all_atoms)
	
	#groups
	#======
	if args.cluster_groups_file!="no":
		#store cluster info in beta factor field
		for p_index in range(0,proteins_nb):
				proteins_sele[p_index].set_bfactor(proteins_cluster_group[p_index][f_index])
		
		#write annotated file
		if args.xtcfilename=="no":
			all_atoms.write(os.getcwd() + '/' + str(args.output_folder) + '/3_snapshots/groups/' + args.grofilename[:-4] + '_annotated_clusterprot_groups', format="PDB")
		else:
			tmp_name=os.getcwd() + "/" + str(args.output_folder) + '/3_snapshots/groups/' + args.xtcfilename[:-4] + '_annotated_clusterprot_groups_' + str(int(t)).zfill(5) + 'ns.pdb'
			W=Writer(tmp_name, nb_atoms)
			W.write(all_atoms)
		
	return
def write_frame_annotation(f_index,t):

	#sizes
	#=====
	#create file
	if args.xtcfilename=="no":
		filename_details=os.getcwd() + '/' + str(args.output_folder) + '/3_snapshots/sizes/' + args.grofilename[:-4] + '_annotated_clusterprot_sizes.txt'
	else:
		filename_details=os.getcwd() + "/" + str(args.output_folder) + '/3_snapshots/sizes/' + args.xtcfilename[:-4] + '_annotated_clusterprot_sizes_' + str(int(t)).zfill(5) + 'ns.txt'
	output_stat = open(filename_details, 'w')		

	#output VMD protein selection line
	tmp_prot_sele=proteins_sele_string_VMD[0]
	for p_index in range(1,proteins_nb):
		tmp_prot_sele+="." + proteins_sele_string_VMD[p_index]
	output_stat.write(tmp_prot_sele + "\n")
				
	#ouput min and max size
	output_stat.write(str(numpy.min(proteins_cluster_size_mat[:,f_index])) + "." + str(numpy.max(proteins_cluster_size_mat[:,f_index])) + "\n")
	
	#ouptut cluster size for each protein
	tmp_sizes="1"
	for p_index in range(0,proteins_nb):
		tmp_sizes+="." + str(proteins_cluster_size_mat[p_index,f_index])
	output_stat.write(tmp_sizes + "\n")
	output_stat.close()
	
	#groups
	#======
	if args.cluster_groups_file!="no":
		#create file
		if args.xtcfilename=="no":
			filename_details=os.getcwd() + '/' + str(args.output_folder) + '/3_snapshots/groups/' + args.grofilename[:-4] + '_annotated_clusterprot_groups.txt'
		else:
			filename_details=os.getcwd() + "/" + str(args.output_folder) + '/3_snapshots/groups/' + args.xtcfilename[:-4] + '_annotated_clusterprot_groups_' + str(int(t)).zfill(5) + 'ns.txt'
		output_stat = open(filename_details, 'w')		
	
		#output VMD protein selection line
		tmp_prot_sele=proteins_sele_string_VMD[0]
		for p_index in range(1,proteins_nb):
			tmp_prot_sele+="." + proteins_sele_string_VMD[p_index]
		output_stat.write(tmp_prot_sele + "\n")
		
		#ouput min and max size
		output_stat.write(str(numpy.min(proteins_cluster_group_mat[:,f_index])) + "." + str(numpy.max(proteins_cluster_group_mat[:,f_index])) + "\n")
		
		#ouptut cluster size for each protein
		tmp_groups="1"
		for p_index in range(0,proteins_nb):
			tmp_groups+="." + str(proteins_cluster_group_mat[p_index,f_index])
		output_stat.write(tmp_groups + "\n")
		output_stat.close()
	
	return
def write_xtc_snapshots():
	
	#NB: - this will always output the first and final frame snapshots
	#    - it will also intermediate frames according to the -w option	

	loc_nb_frames_processed=0
	for ts in U.trajectory:

		#case: frames before specified time boundaries
		#---------------------------------------------
		if ts.time/float(1000)<args.t_start:
			progress='\r -skipping frame ' + str(ts.frame) + '/' + str(nb_frames_xtc) + '        '
			sys.stdout.flush()
			sys.stdout.write(progress)

		#case: frames within specified time boundaries
		#---------------------------------------------
		elif ts.time/float(1000)>args.t_start and ts.time/float(1000)<args.t_end:
			progress='\r -writing snapshots...   frame ' + str(ts.frame) + '/' + str(nb_frames_xtc) + '        '
			sys.stdout.flush()
			sys.stdout.write(progress)
			if ((ts.frame-1) % args.frames_dt)==0:
				if ((loc_nb_frames_processed) % args.frames_write_dt)==0 or loc_nb_frames_processed==nb_frames_processed-1:
					write_frame_stat(ts.frame, loc_nb_frames_processed, ts.time/float(1000))
					write_frame_snapshot(loc_nb_frames_processed, ts.time/float(1000))
					write_frame_annotation(loc_nb_frames_processed, ts.time/float(1000))
				loc_nb_frames_processed+=1
		
		#case: frames after specified time boundaries
		#---------------------------------------------
		elif ts.time/float(1000)>args.t_end:
			break

	print ''

	return
def write_xtc_annotation():
	
	#sizes
	#=====
	#create file
	filename_details=os.getcwd() + '/' + str(args.output_folder) + '/4_VMD/' + args.xtcfilename[:-4] + '_annotated_clustprot_sizes.txt'
	output_stat = open(filename_details, 'w')		

	#output VMD protein selection line
	tmp_prot_sele=proteins_sele_string_VMD[0]
	for p_index in range(1,proteins_nb):
		tmp_prot_sele+="." + proteins_sele_string_VMD[p_index]
	output_stat.write(tmp_prot_sele + "\n")
	
	#ouput min and max size
	output_stat.write(str(min(proteins_sizes_sampled)) + "." + str(max(proteins_sizes_sampled)) + "\n")
	
	#ouptut cluster size for each protein
	proteins_cluster_size_mat=numpy.asarray(proteins_cluster_size.values())
	for frame in sorted(time_stamp.iterkeys()):
		tmp_sizes=str(frame)
		frame_index=sorted(time_stamp.keys()).index(frame)
		for p_index in range(0,proteins_nb):
			tmp_sizes+="." + str(proteins_cluster_size_mat[p_index,frame_index])
		output_stat.write(tmp_sizes + "\n")
	output_stat.close()
	
	#groups
	#======
	if args.cluster_groups_file!="no":
		#create file
		filename_details=os.getcwd() + '/' + str(args.output_folder) + '/4_VMD/' + args.xtcfilename[:-4] + '_annotated_clustprot_groups.txt'
		output_stat = open(filename_details, 'w')		
	
		#output VMD protein selection line
		tmp_prot_sele=proteins_sele_string_VMD[0]
		for p_index in range(1,proteins_nb):
			tmp_prot_sele+="." + proteins_sele_string_VMD[p_index]
		output_stat.write(tmp_prot_sele + "\n")
		
		#ouput min and max size
		output_stat.write(str(min(proteins_groups_sampled)) + "." + str(max(proteins_groups_sampled)) + "\n")
		
		#ouptut cluster size for each protein
		proteins_cluster_group_mat=numpy.asarray(proteins_cluster_group.values())
		for frame in sorted(time_stamp.iterkeys()):
			tmp_groups=str(frame)
			frame_index=sorted(time_stamp.keys()).index(frame)
			for p_index in range(0,proteins_nb):
				tmp_groups+="." + str(proteins_cluster_group_mat[p_index,frame_index])
			output_stat.write(tmp_groups + "\n")
		output_stat.close()

	return

################################################################################################################################################
# DATA STRUCTURES
################################################################################################################################################

#time
time_stamp={}

#cluster size/size group each lipid is involved in at each frame
sizes_nb={}
sizes_pc={}
groups_nb={}
groups_pc={}
proteins_sizes_sampled=[]
proteins_groups_sampled=[]
proteins_cluster_size={}
proteins_cluster_group={}
proteins_cluster_size_mat=numpy.zeros((proteins_nb,1))
proteins_cluster_group_mat=numpy.zeros((proteins_nb,1))
for p_index in range(0,proteins_nb):
	proteins_cluster_size[p_index]=[]
	proteins_cluster_group[p_index]=[]
	
sizes_colors_nb=0
sizes_colors_dict={}
sizes_colors_list=[]

#sizes sampled by each specie
if args.xtcfilename!="no":
	time_sorted=[]
	groups_stability={}
	biggest_size={}
	biggest_nb={}
	biggest_pc={}
	biggest_size_sorted=[]
	biggest_nb_sorted=[]
	biggest_pc_sorted=[]
	
#smooth data
if args.nb_smoothing>1:
	time_smoothed=[]
	biggest_size_smoothed={}
	biggest_nb_smoothed={}
	biggest_pc_smoothed={}	
	sizes_nb_smoothed={}
	sizes_pc_smoothed={}
	groups_nb_smoothed={}
	groups_pc_smoothed={}

################################################################################################################################################
# ALGORITHM : Browse trajectory and process relevant frames
################################################################################################################################################

print "\nDetecting lipid clusters..."

#case: gro file
#==============
if args.xtcfilename=="no":
	#store dummy time
	time_stamp[1]=0
							
	#detect clusters
	if args.m_algorithm!="density":
		tmp_groups=detect_clusters_connectivity(get_distances(U.trajectory.ts.dimensions), U.dimensions)
	else:
		tmp_groups=detect_clusters_density(get_distances(U.trajectory.ts.dimensions), U.dimensions)
	
	#case: store cluster size only for each protein
	if args.cluster_groups_file=="no":
		for g in tmp_groups:
			for p_index in g:
				proteins_cluster_size[p_index].append(numpy.size(g))
	#case: store cluster size and group size for each protein
	else:
		for g in tmp_groups:
			for p_index in g:
				proteins_cluster_size[p_index].append(numpy.size(g))
				proteins_cluster_group[p_index].append(groups_sizes_dict[numpy.size(g)])
			
#case: xtc file
#==============
else:
	for ts in U.trajectory:

		#case: frames before specified time boundaries
		#---------------------------------------------
		if ts.time/float(1000)<args.t_start:
			progress='\r -skipping frame ' + str(ts.frame) + '/' + str(nb_frames_xtc) + '        '
			sys.stdout.flush()
			sys.stdout.write(progress)

		#case: frames within specified time boundaries
		#---------------------------------------------
		elif ts.time/float(1000)>args.t_start and ts.time/float(1000)<args.t_end:
			progress='\r -processing frame ' + str(ts.frame) + '/' + str(nb_frames_xtc) + '        '
			sys.stdout.flush()
			sys.stdout.write(progress)
			if ((ts.frame-1) % args.frames_dt)==0:
				nb_frames_processed+=1
				#store time
				time_stamp[ts.frame]=ts.time/float(1000)
						
				#detect clusters
				if args.m_algorithm!="density":
					tmp_groups=detect_clusters_connectivity(get_distances(U.trajectory.ts.dimensions), U.dimensions)
				else:
					tmp_groups=detect_clusters_density(get_distances(U.trajectory.ts.dimensions), U.dimensions)
				
				#case: store cluster size only for each lipd
				if args.cluster_groups_file=="no":
					for g in tmp_groups:
						for p_index in g:
							proteins_cluster_size[p_index].append(numpy.size(g))
				#case: store cluster size and group size for each lipd
				else:
					for g in tmp_groups:
						for p_index in g:
							proteins_cluster_size[p_index].append(numpy.size(g))
							proteins_cluster_group[p_index].append(groups_sizes_dict[numpy.size(g)])
		
		#case: frames after specified time boundaries
		#--------------------------------------------
		elif ts.time/float(1000)>args.t_end:
			break
									
	print ""

################################################################################################################################################
# CALCULATE STATISTICS
################################################################################################################################################

print "\nCalculating statistics..."
get_sizes_sampled()
update_color_dict()
calc_stat()
if args.xtcfilename!="no":
	smooth_data()

################################################################################################################################################
# PRODUCE OUTPUTS
################################################################################################################################################

print "\nWriting outputs..."

#case: gro file
#==============
if args.xtcfilename=="no":
	if len(proteins_sizes_sampled)>1:
		print " -writing statistics..."
		write_frame_stat(1,"all","all")
		print " -writing annotated pdb..."
		write_frame_snapshot(0,0)
		write_frame_annotation(0,0)
	else:
		print "\n"
		print "Warning: a single cluster size (", str(proteins_sizes_sampled[0]), ") was detected throughout the trajectory. Check the --algorithm, --cutoff, --radius or --neighbours options (see cluster_prot -h)."
		write_warning()

#case: xtc file
#==============
else:
	if len(proteins_sizes_sampled)>1:
		#writing statistics
		print " -writing statistics..."
		write_frame_stat(1,"all","all")
		#output cluster snaphots
		write_xtc_snapshots()
		#write annotation files for VMD
		print " -writing VMD annotation files..."
		write_xtc_annotation()
		#write xvg and graphs
		print " -writing xvg and graphs..."
		graph_aggregation_2D_sizes()
		write_xvg_biggest()
		graph_xvg_biggest()
		write_xvg_sizes()
		graph_xvg_sizes()
		if args.nb_smoothing>1:
			write_xvg_sizes_smoothed()
			graph_xvg_sizes_smoothed()
			write_xvg_biggest_smoothed()
			graph_xvg_biggest_smoothed()
		if args.cluster_groups_file!="no":
			graph_aggregation_2D_groups()
			write_stability_groups()
			write_xvg_groups()
			graph_xvg_groups()
			if args.nb_smoothing>1:
				write_xvg_groups_smoothed()
				graph_xvg_groups_smoothed()
	else:
		print "\n"
		print "Warning: a single cluster size (", str(proteins_sizes_sampled[0]), ") was detected throughout the trajectory. Check the -m, -c, -r or -n options (see cluster_prot -h)."
		write_warning()
		
#exit
#====
print "\nFinished successfully! Check output in ./" + args.output_folder + "/"
print ""
sys.exit(0)
