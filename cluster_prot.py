#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

################################################################################################################################################
# RETRIEVE USER INPUTS
################################################################################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.1.2-dev5"
parser = argparse.ArgumentParser(prog='cluster_prot', usage='', add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, description=\
'''
**********************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/cluster_prot
**********************************************
	
[ DESCRITPION ]

This script calculates the evolution of proteins clustering status.
	
It produces 3 outputs ((a) and (b) require --groups to be specified, see note 5):
 (a) 2D plots: time evolution of the cluster size each protein is involved in
 (b) 1D plots: time evolution of the % of protein represented by each group
 (c) stability statistics: max nb of consecutive frames each group existed for

The metrics can be calculated for a single frame or for an entire trajectory - and in case
a trajectory is supplied the data for individual frame snapshots can also be produced at a
frequency specified by the option -w (if you don't specify an argument to -w snapshots
will only be written for the first and last frame). Note that writing files to disk will
considerably slow down the execution of the script.

Detection of transmembrane protein clusters
-------------------------------------------
Two clustering algorithms can be used to identify protein clusters.
->Connectivity based (relies on networkX module):
  A protein is considered in a cluster if it is within a distance less than --nx_cutoff
  from another protein. This means that a single protein can act as a connector between
  two otherwise disconnected protein clusters.
  This algorithm can be ran using either the minimum distante between proteins (default, 
  --algorithm 'min') or the distance between their center of geometry (--algorithm 'cog').
  The 'min' option scales as the square of the number of proteins and can thus be very
  slow for large systems.

->Density based (relies on the sklearn module and its implementation of DBSCAN):
  A protein is considered in a cluster if is surrounded by at least --db_neighbours other
  proteins within a radius of --db_radius.
  This density based approach is usually less suited to the detection of protein
  clusters but as a general rule the more compact the clusters, the smaller --db_radius
  the higher --db_neighbours can be - for details on this algorithm see its online
  documentation.
  This algorithm is selected by setting the --algorithm option to 'density'.

The identified protein clusters are considered to be transmembrane only if the closest
lipid headgroup neighbours to the cluster particles are all within the same leaflet.
In addition to the sizes identified, size groups can be defined - see note 7.

VMD visualisation
-----------------
The perturbation calculated can be visualised with VMD either for a single frame or an
entire trajectory. Note that in the case of a trajectory only the processed frame will be
annotated (every 10 by defaults) - you might want to pre-process your trajectory to remove
frames and then set the -t option to 1 when running the script.
 ->frame:
   The thickness or order parameter info is stored in the beta factor column. Just open
   the PDB file with VMD and choose Draw Style > Coloring Method > Beta.

 ->trajectory:
   The thickness or order parameter info is stored in a .txt file in folders '/3_VMD/'.
   To load it into your trajectory in VMD use the appropriate routine of the script
   script 'vmd_parser.tcl' (see https://github.com/jhelie/vmd_parser).


[ REQUIREMENTS ]

The following python modules are needed :
 - MDAnalysis
 - matplotlib
 - networkX (if --algorithm set to 'min' or 'cog')
 - sklearn (if --algorithm set to 'density')


[ NOTES ]

1. It's a good idea to trjconv the xtc first and only outputs the proteins,  as the 
   script will run MUCH faster. Also, use the -pbc mol option.	

2. If you're interested to find out whether proteins go from a TM orientation to an
   interfacial orientation, the script can identify bilayer leafletd and track the
   accummulation of proteins on each leaflet (but clustering information is not calculated
   for interfacial proteins).

   To identify leaflets and track interfacial proteins you must specify the --leaflets
   option.
   
   Identification of the bilayer leaflets can be further controlled via 3 ways:
   (a) beads
    By default, the particles taken into account to define leaflet are:e
    -> name PO4 or name PO3 or name B1A
   
    Note that only lipids which contain one of the beads mentioned in the selection string
    will be taken into account. If you wish to specify your own selection string (e.g. to
    choose different beads or add a bead not in the default list in order to take into
    account a particular lipid specie) you can do so by supplying a file via the --beads
    option. This file should contain a single line that can be passed as the argument
    to MDAnalysis selectAtoms() routine and should not contain any quotation marks, e.g.:
     -> name PO4 or name PO3 or name B1A or name AM1
        
   (b) leaflet finding method
    By default leaflets are identified using the MDAnalysis LeafletFinder routine and the
    the optimum cutoff to identify 2 lipids groups is determined using the optimize_cutoff
    routine.
    This optimisation process can take time in large systems and you can specify your own
    cutoff value to skip this step. For instance to use a 15 Angstrom cutoff value:
     -> '--leaflet 15'
   
    In very large systems (more then ~50,000 phospholipids) LeafletFinder (or rather the
    networkX module that it relies on) can fail. To  avoid this you can choose not to use
    this routine by specifying:
     -> '--leaflet large'
    In this case lipids whose headgroups z value is above the average lipids z value will
    be considered to make up the upper leaflet and those whose headgroups z value is below
    the average will be considered to be in the lower leaflet.
    This means that the bilayer should be as flat as possible in the gro file supplied in
    order to get a meaningful outcome.

   (c) flipflopping lipids
    In case lipids flipflop during the trajectory, a file listing them can be supplied
    with the --flipflops option. Each line of this file should follow the format:
     -> 'resname,resid,starting_leaflet,z_bead'
    where starting_leaflet is either 'upper' or 'lower' - e.g. 'POPC,145,lower,PO4'. The
    z_bead is used to track the position of the lipid.
    If flipflopping lipids are not specified they may add significant noise to results.   

3. Proteins are detected automatically but you can specify an input file to define your
   own selection with the -p option.
   In this case the supplied file should contain on each line a protein selection string
   that can be passed as the argument of the MDAnalysis selectAtoms() routine - for 
   instance 'bynum 1:344'.

4. Clusters statistics can be binned into size groups. The size groups are defined by
   supplying a file with the --groups option, whose lines all follow the format:
    -> 'lower_size,upper_size, colour'

   Size groups definition should follow the following rules:
    -to specify an open ended group use 'max', e.g. '3,max,colour'
    -groups should be ordered by increasing size and their boundaries should not overlap
    -boundaries are inclusive so you can specify one size groups with 'size,size,colour'
    -colours must be specified for each group (see (b) in note 8)
    -any cluster size not falling within the specified size groups will be labeled as
     'other' and coloured in grey (#E6E6E6)
	
5. The colours associated to each TM cluster size identified are based on the matplotlib
   'jet' colour map by default. You can specify your own colours as follows:
   (a) Colour of individual cluster sizes
    Colours of individual cluster sizes use the matplotlib 'jet' colour scheme and cannot 
    be modified. 

   (b) Cluster sizes colours
    Colours of individual cluster sizes use the matplotlib 'jet' colour scheme and cannot 
    be modified. You need to specify a size range via --colours_sizes on which to apply
    the colour map (this is because the script does eveything in one pass and cannot guess
    what the cluster sizes will be in advance).
     -> '--colours_sizes 2,6'
    Cluster sizes below the range will use the same colour as that of the lower size and
    those above the range will use the same colour as that of the upper size.

    Exact colour control can be achieved by using size groups (see note 4), so if you want
    to control individual cluster sizes colours just specify the relevant group file.

   (c) Colour definition
    Colours can be specified using single letter code (rgbcmykw), hex code  or the name of
    a colour map (see the matplotlib website for a list of the available colour maps).
    In case a colour map is used, its name must be specified as the colour for each lipid
    specie or size group.


[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: structure file [.gro] (required)
-x			: trajectory file [.xtc]
-o			: name of output folder
-b			: beginning time (ns)
-e			: ending time (ns)	
-t 		10	: process every t-frames
-w			: write snapshots every [w] processed frames (see 'DESCRIPTION')
--smooth		: nb of points to use for data smoothing

Lipids identification (see note 2)
-----------------------------------------------------
--beads			: leaflet identification technique, see note 2(a)
--flipflops		: input file with flipflopping lipids, see note 2(c)
--leaflets	no	: leaflet identification, see note 2(b)

Protein clusters identification
-----------------------------------------------------
--groups		: cluster groups definition file, see note 4
--proteins		: protein selection file, (optional, see note 3)
--colours_sizes	1,9	: range of cluster sizes to colour, see note 5
--algorithm	min	: 'cog','min' or 'density', see 'DESCRIPTION'
--nx_cutoff 	8	: networkX cutoff distance for protein-protein contact (Angstrom)
--db_radius 	20	: DBSCAN search radius (Angstrom)
--db_neighbours	3	: DBSCAN minimum number of neighbours within a circle of radius --db_radius	
 
Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
  
''')

#data options
parser.add_argument('-f', nargs=1, dest='grofilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-x', nargs=1, dest='xtcfilename', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-o', nargs=1, dest='output_folder', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-b', nargs=1, dest='t_start', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-e', nargs=1, dest='t_end', default=[10000000000000], type=int, help=argparse.SUPPRESS)
parser.add_argument('-t', nargs=1, dest='frames_dt', default=[10], type=int, help=argparse.SUPPRESS)
parser.add_argument('-w', nargs='?', dest='frames_write_dt', const=1000000000000000, default="no", help=argparse.SUPPRESS)
parser.add_argument('--smooth', nargs=1, dest='nb_smoothing', default=[0], type=int, help=argparse.SUPPRESS)

#lipids identification
parser.add_argument('--beads', nargs=1, dest='beadsfilename', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('--flipflops', nargs=1, dest='selection_file_ff', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('--leaflets', nargs=1, dest='cutoff_leaflet', default=['optimise'], help=argparse.SUPPRESS)

#protein clusters identification

parser.add_argument('--groups', nargs=1, dest='cluster_groups_file', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('--proteins', nargs=1, dest='selection_file_prot', default=['auto'], help=argparse.SUPPRESS)
parser.add_argument('--colours_sizes', nargs=1, dest='colours_sizes', default=['1,9'], help=argparse.SUPPRESS)
parser.add_argument('--algorithm', dest='m_algorithm', choices=['cog','min','density'], default='min', help=argparse.SUPPRESS)
parser.add_argument('--nx_cutoff', nargs=1, dest='cutoff_connect', default=[8], type=float, help=argparse.SUPPRESS)
parser.add_argument('--db_radius', nargs=1, dest='dbscan_dist', default=[20], type=float, help=argparse.SUPPRESS)
parser.add_argument('--db_neighbours', nargs=1, dest='dbscan_nb', default=[3], type=int, help=argparse.SUPPRESS)

#other options
parser.add_argument('--buffer_size', nargs=1, dest='buffer_size', default=[100], type=int, help=argparse.SUPPRESS)
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

#parse user inputs
#-----------------
args = parser.parse_args()
#data options
args.grofilename = args.grofilename[0]
args.xtcfilename = args.xtcfilename[0]
args.output_folder = args.output_folder[0]
args.t_start=args.t_start[0]
args.t_end=args.t_end[0]
args.frames_dt=args.frames_dt[0]
if args.frames_write_dt != "no":
	args.frames_write_dt = int(args.frames_write_dt)
args.nb_smoothing = args.nb_smoothing[0]
#lipids identification
args.beadsfilename = args.beadsfilename[0]
args.cutoff_leaflet = args.cutoff_leaflet[0]
args.selection_file_ff = args.selection_file_ff[0]
#protein clusters identification
args.cluster_groups_file = args.cluster_groups_file[0]
args.selection_file_prot = args.selection_file_prot[0]
args.colours_sizes = args.colours_sizes[0]
args.cutoff_connect = args.cutoff_connect[0]
args.dbscan_dist = args.dbscan_dist[0]
args.dbscan_nb = args.dbscan_nb[0]
#other options
args.buffer_size = args.buffer_size[0]

#process options
#---------------
global vmd_counter
global vmd_cluster_size
global vmd_cluster_group
global lipids_ff_nb
global colours_sizes_range
global colour_group_other
global colour_leaflet_lower
global colour_leaflet_upper

vmd_counter = 0
vmd_cluster_size = ""
vmd_cluster_group = ""
lipids_ff_nb = 0

colour_group_other = "#E6E6E6"											#very light grey
colour_leaflet_lower = "#808080"										#dark grey
colour_leaflet_upper = "#C0C0C0"										#light grey

#colour of cluster sizes
tmp_col_size = args.colours_sizes.split(',')
if len(tmp_col_size) != 2:
	print "Error: wrong format for the option --colours_sizes, it should be 'min,max' (see bilayer_perturbations --help, note 8)."
	sys.exit(1)
elif int(tmp_col_size[0]) > int(tmp_col_size[1]):
	print "Error: wrong format for the option --colours_sizes, it should be 'min,max' (see bilayer_perturbations --help, note 8)."
	sys.exit(1)
else:
	colours_sizes_range = [int(tmp_col_size[0]), int(tmp_col_size[1])]
	
#leaflet identification
if args.cutoff_leaflet != "no":
	if args.cutoff_leaflet != "large" and args.cutoff_leaflet != "no" and args.cutoff_leaflet != "optimise":
		try:
			args.cutoff_leaflet = float(args.cutoff_leaflet)
		except:
			print "Error: the argument of the --leaflets option should be a number or 'large', see note 2"
			sys.exit(1)


#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================
#generic science modules
try:
	import math
except:
	print "Error: you need to install the maths module."
	sys.exit(1)
try:
	import numpy as np
except:
	print "Error: you need to install the np module."
	sys.exit(1)
try:
	import scipy
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.colors as mcolors
	mcolorconv = mcolors.ColorConverter()
	import matplotlib.cm as cm				#colours library
	import matplotlib.ticker
	from matplotlib.ticker import MaxNLocator
	from matplotlib.font_manager import FontProperties
	fontP=FontProperties()
except:
	print "Error: you need to install the matplotlib module."
	sys.exit(1)
try:
	import pylab as plt
except:
	print "Error: you need to install the pylab module."
	sys.exit(1)

#MDAnalysis module
try:
	import MDAnalysis
	from MDAnalysis import *
	import MDAnalysis.analysis
	import MDAnalysis.analysis.leaflet
	import MDAnalysis.analysis.distances
	#set MDAnalysis to use periodic boundary conditions
	MDAnalysis.core.flags['use_periodic_selections'] = True
	MDAnalysis.core.flags['use_KDTree_routines'] = False
except:
	print "Error: you need to install the MDAnalysis module first. See http://mdanalysis.googlecode.com"
	sys.exit(1)

#=========================================================================================
# sanity check
#=========================================================================================
if not os.path.isfile(args.grofilename):
	print "Error: file " + str(args.grofilename) + " not found."
	sys.exit(1)
if args.cluster_groups_file != "no" and not os.path.isfile(args.cluster_groups_file):
	print "Error: file " + str(args.cluster_groups_file) + " not found."
	sys.exit(1)
if args.selection_file_ff != "no" and not os.path.isfile(args.selection_file_ff):
	print "Error: file " + str(args.selection_file_ff) + " not found."
	sys.exit(1)
if args.cluster_groups_file != "no" and not os.path.isfile(args.cluster_groups_file):
	print "Error: file " + str(args.cluster_groups_file) + " not found."
	sys.exit(1)
if args.beadsfilename != "no" and not os.path.isfile(args.beadsfilename):
	print "Error: file " + str(args.beadsfilename) + " not found."
	sys.exit(1)
if args.t_end < args.t_start:
	print "Error: the starting time (" + str(args.t_start) + "ns) for analysis is later than the ending time (" + str(args.t_end) + "ns)."
	sys.exit(1)
if args.buffer_size < -1:
	print "Error: the option --buffer_size should be greater than -1 (set to " + str(args.buffer_size) + ")."
	sys.exit(1)

if args.xtcfilename == "no":
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

if args.m_algorithm != "density":
	if '--db_radius' in sys.argv:
		print "Error: --db_radius option specified but --algorithm option set to '" + str(args.m_algorithm) + "'."
		sys.exit(1)
	elif '--db_neighbours' in sys.argv:
		print "Error: --db_neighbours option specified but --algorithm option set to '" + str(args.m_algorithm) + "'."
		sys.exit(1)
else:
	if '--nx_cutoff' in sys.argv:
		print "Error: --nx_cutoff option specified but --algorithm option set to 'density'."
		sys.exit(1)

if args.cutoff_leaflet != "no" and args.cutoff_leaflet != "optimise":
	try:
		args.cutoff_leaflet = float(args.cutoff_leaflet)
	except:
		print "Error: the argument of the --leaflets option should be a number or 'large', see note 2"
		sys.exit(1)

#=========================================================================================
# create folders and log file
#=========================================================================================
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
	if args.selection_file_ff != "no":
		shutil.copy2(args.selection_file_ff,args.output_folder + "/")	
	if args.selection_file_prot != "no" and args.selection_file_prot != "auto":
		shutil.copy2(args.selection_file_prot,args.output_folder + "/")
	if args.cluster_groups_file != "no":
		shutil.copy2(args.cluster_groups_file,args.output_folder + "/")
	if args.beadsfilename != "no":
		shutil.copy2(args.beadsfilename,args.output_folder + "/")

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================

def set_lipids_beads():													#DONE

	global leaflet_sele_string

	#set default beads
	leaflet_sele_string = "name PO4 or name PO3 or name B1A"

	#use users input
	if args.beadsfilename != "no":
		with open(args.beadsfilename) as f:
			lines = f.readlines()
		if len(lines) > 1:
			print "Error: the file " + str(args.beadsfilename) + " should conly ontain 1 line (" + str(len(lines)) + " found), see note 2(a)."
			sys.exit(1)
		else:
			if lines[0][-1] == "\n":
				lines[0] = lines[0][:-1]
			leaflet_sele_string = lines[0]

	return
def load_MDA_universe():												#DONE
	
	global U
	global all_atoms
	global nb_atoms
	global nb_frames_xtc
	global frames_to_process
	global frames_to_write
	global nb_frames_to_process
	global f_start
	f_start = 0
	if args.xtcfilename == "no":
		print "\nLoading file..."
		U = Universe(args.grofilename)
		all_atoms = U.selectAtoms("all")
		nb_atoms = all_atoms.numberOfAtoms()
		nb_frames_xtc = 1
		frames_to_process = [0]
		frames_to_write = [True]
		nb_frames_to_process = 1
	else:
		print "\nLoading trajectory..."
		U = Universe(args.grofilename, args.xtcfilename)
		all_atoms = U.selectAtoms("all")
		nb_atoms = all_atoms.numberOfAtoms()
		nb_frames_xtc = U.trajectory.numframes
		U.trajectory.rewind()
		#sanity check
		if U.trajectory.time/float(1000) < args.t_start:
			print "Error: the trajectory duration (" + str(U.trajectory.time/float(1000)) + "ns) is shorted than the starting stime specified (" + str(args.t_start) + "ns)."
			sys.exit(1)
		if U.trajectory.numframes < args.frames_dt:
			print "Warning: the trajectory contains fewer frames (" + str(nb_frames_xtc) + ") than the frame step specified (" + str(args.frames_dt) + ")."

		#create list of index of frames to process
		if args.t_start > 0:
			for ts in U.trajectory:
				#debug
				print ts.frame
				progress = '\r -skipping frame ' + str(ts.frame) + '/' + str(nb_frames_xtc) + '        '
				sys.stdout.flush()
				sys.stdout.write(progress)
				if ts.time/float(1000) < args.t_start:
					f_start = ts.frame-1
					break
		if (nb_frames_xtc - f_start)%args.frames_dt == 0:
			tmp_offset = 0
		else:
			tmp_offset = 1
		frames_to_process = map(lambda f:f_start + args.frames_dt*f, range(0,(nb_frames_xtc - f_start)//args.frames_dt+tmp_offset))
		nb_frames_to_process = len(frames_to_process)
			
		#create list of frames to write
		if args.frames_write_dt == "no":
			frames_to_write = [False for f_index in range(0, nb_frames_to_process)]
		else:
			frames_to_write = [True if (f_index % args.frames_write_dt == 0 or f_index == (nb_frames_xtc - f_start)//args.frames_dt) else False for f_index in range(0, nb_frames_to_process)]

	#check for the presence of proteins
	test_prot = U.selectAtoms("protein")
	if test_prot.numberOfAtoms() == 0:
			print "Error: no protein found in the system."
			sys.exit(1)
			
	#check for the presence of lipids
	if args.cutoff_leaflet != "no":
		test_beads = U.selectAtoms(leaflet_sele_string)
		if test_beads.numberOfAtoms() == 0:
			print "Error: invalid selection string '" + str(leaflet_sele_string) + "'"
			print "-> no lipid particles selected, check the --leaflets and/or --beads options."
			sys.exit(1)

	return
def identify_ff():														#DONE
	print "\nReading selection file for flipflopping lipids..."
	
	#declare variables
	global lipids_ff_nb
	global lipids_ff_info
	global lipids_ff_resnames
	global lipids_ff_leaflet
	global lipids_ff_u2l_index
	global lipids_ff_l2u_index
	global lipids_sele_ff
	global lipids_sele_ff_bead
	global lipids_sele_ff_bonds
	global lipids_sele_ff_VMD_string
	global leaflet_sele_string
	lipids_ff_nb = 0
	lipids_ff_info = {}
	lipids_ff_resnames = []
	lipids_ff_leaflet = []
	lipids_ff_u2l_index = []
	lipids_ff_l2u_index = []
	lipids_sele_ff = {}
	lipids_sele_ff_bead = {}
	lipids_sele_ff_bonds = {}
	lipids_sele_ff_VMD_string={}
		
	with open(args.selection_file_ff) as f:
		lines = f.readlines()
	lipids_ff_nb = len(lines)
	print " -found " + str(lipids_ff_nb) + " flipflopping lipids"
	leaflet_sele_string = leaflet_sele_string + " and not ("
	for l_index in range(0,lipids_ff_nb):
		line = lines[l_index]
		if line[-1] == "\n":
			line = line[:-1]
		try:
			line_content = line.split(',')
			if len(line_content) != 4:
				print "Error: wrong format for line " + str(l_index+1) + " in " + str(args.selection_file_ff) + ", see note 4 in bilayer_perturbations --help."
				print " ->", line
				sys.exit(1)
			#read current lipid details
			lip_resname = line_content[0]
			lip_resnum = int(line_content[1])
			lip_leaflet = line_content[2]
			lip_bead = line_content[3]
			lipids_ff_info[l_index] = [lip_resname,lip_resnum,lip_leaflet,lip_bead]
						
			#update: starting leaflets
			if lip_leaflet not in lipids_ff_leaflet:
				lipids_ff_leaflet.append(lip_leaflet)

			#update: index in directional lists
			if lip_leaflet == "upper":
				lipids_ff_u2l_index.append(l_index)
			elif lip_leaflet == "lower":
				lipids_ff_l2u_index.append(l_index)
			else:
				print "->unknown starting leaflet '" + str(lip_leaflet) + "'."
				sys.exit(1)
			
			#update: resnames
			if lip_resname not in lipids_ff_resnames:
				lipids_ff_resnames.append(lip_resname)
	
			#update: leaflet selection string
			if l_index==0:
				leaflet_sele_string+="(resname " + str(lip_resname) + " and resnum " + str(lip_resnum) + ")"
			else:
				leaflet_sele_string+=" or (resname " + str(lip_resname) + " and resnum " + str(lip_resnum) + ")"

			#create selections
			lipids_sele_ff[l_index] = U.selectAtoms("resname " + str(lip_resname) + " and resnum " + str(lip_resnum))
			lipids_sele_ff_bead[l_index] = lipids_sele_ff[l_index].selectAtoms("name " + str(lip_bead))
			lipids_sele_ff_VMD_string[l_index]="resname " + str(lipids_ff_info[l_index][0]) + " and resid " + str(lipids_ff_info[l_index][1])
			if lipids_sele_ff[l_index].numberOfAtoms() == 0:
				print "Error:"
				print line
				print "-> no such lipid found."
				sys.exit(1)	
		except:
			print "Error: invalid flipflopping lipid selection string on line " + str(l_index+1) + ": '" + line + "'"
			sys.exit(1)
	leaflet_sele_string+=")"		

	return
def identify_proteins():												#DONE
	print "\nIdentifying proteins..."
	
	#import modules
	if args.m_algorithm == "density":
		global DBSCAN
		from sklearn.cluster import DBSCAN
	else:
		global nx
		import networkx as nx

	#declare variables
	global proteins_nb
	global proteins_sele
	global proteins_sele_string
	global proteins_sele_string_VMD
	global proteins_boundaries
	global proteins_nb_atoms
	global nb_atom_per_protein
	proteins_nb = 0
	proteins_sele = {}
	proteins_sele_string = {}
	proteins_sele_string_VMD = {}
	proteins_boundaries = {}
	
	#check for protein presence
	if U.selectAtoms("protein").numberOfAtoms() == 0:
		print "Error: no protein detected."
		sys.exit(1)
	
	#case: selection file provided
	if args.selection_file_prot != "auto":
		print " -reading protein selection file..."
		with open(args.selection_file_prot) as f:
			lines = f.readlines()
		proteins_nb=len(lines)
		proteins_sele["all"] = MDAnalysis.core.AtomGroup.AtomGroup([])
		for p_index in range(0,proteins_nb):
			line = lines[p_index]
			if line[-1] == "\n":
				line = line[:-1]
			progress='\r -creating proteins selections: ' + str(p_index+1) + '/' + str(proteins_nb) + '        '
			sys.stdout.flush()
			sys.stdout.write(progress)
			try:
				print " p[" + str(p_index) + "]=U.selectAtoms(" + line + ")"
				proteins_sele[p_index] = U.selectAtoms(line[1:-2])
				proteins_sele["all"] += proteins_sele[p_index]
				proteins_boundaries[p_index] = [proteins_sele[p_index].indices()[0] + 1, proteins_sele[p_index].indices()[proteins_sele[p_index].numberOfAtoms()]+1]
				proteins_sele_string[p_index] = "bynum " + str(proteins_boundaries[p_index][0]) + ":" + str(proteins_boundaries[p_index][1])
				proteins_sele_string_VMD[p_index] = "serial " + str(proteins_boundaries[p_index][0]) + " to " + str(proteins_boundaries[p_index][1])
			except:
				print "Error:"
				print line
				print "->invalid selection string."
				sys.exit(1)
		proteins_nb_atoms = proteins_sele["all"].numberOfAtoms()
	
	#case: automatic detection
	else:
		#declare local variables
		proteins_ca_nb = {}
		proteins_ca_nmax = 0
		proteins_ca_group = {}
		proteins_boundaries = {}
	
		#retrieve 1st atom info
		proteins_sele["all"] = U.selectAtoms("protein")
		proteins_nb_atoms = proteins_sele["all"].numberOfAtoms()
		prec_resnum = proteins_sele["all"][0].resnum
		prec_segid = proteins_sele["all"][0].segid
		prec_atnum = proteins_sele["all"][0].number+1
		prev_atnum = proteins_sele["all"][0].number+1						#atom corresponding to the beginning of the current protein
		#browse following atoms
		for a in proteins_sele["all"][1:]:
			delta_res = a.resnum-prec_resnum
			delta_atm = a.number+1-prec_atnum
			if delta_res < 0 or a.segid != prec_segid or delta_atm > 1:
				proteins_boundaries[proteins_nb] = [prev_atnum,prec_atnum]
				proteins_nb += 1
				prev_atnum = a.number + 1
			prec_resnum = a.resnum
			prec_atnum = a.number + 1
			prec_segid = a.segid		
		#add last protein section
		if prev_atnum < proteins_sele["all"][proteins_nb_atoms-1].number:
			proteins_boundaries[proteins_nb] = [prev_atnum,proteins_sele["all"][proteins_nb_atoms-1].number+1]
			proteins_nb += 1
		
		#display results
		print " -protein found:", proteins_nb
		print " -protein boundaries (atom numbers): see protein.sele file"
		#create protein selections and save into a txt file
		filename_sele=os.getcwd() + '/' + str(args.output_folder) + '/proteins.sele'
		output_stat = open(filename_sele, 'w')	
		output_stat.write("#This file was generated by the script bilayer_perturbations v" + str(version_nb) +"\n")
		output_stat.write("#The lines below correspond to MDAnalysis section string, e.g. U.selectAtoms(LINE)\n")
		output_stat.write("\n")	
		for p_index in range(0, proteins_nb):
			progress='\r -creating proteins selections: ' + str(p_index+1) + '/' + str(proteins_nb) + '        '
			sys.stdout.flush()
			sys.stdout.write(progress)
			proteins_sele_string[p_index] = "bynum " + str(proteins_boundaries[p_index][0]) + ":" + str(proteins_boundaries[p_index][1])
			proteins_sele_string_VMD[p_index] = "serial " + str(proteins_boundaries[p_index][0]) + " to " + str(proteins_boundaries[p_index][1])
			proteins_sele[p_index] = U.selectAtoms(proteins_sele_string[p_index])
			output_stat.write(proteins_sele_string[p_index] + "\n")
		output_stat.close()

	nb_atom_per_protein = proteins_sele[0].numberOfAtoms()
	print ""

	return
def identify_leaflets():												#DONE
	print "\nIdentifying leaflets..."
	
	#declare variables
	global leaflet_sele
	global leaflet_sele_atoms
	leaflet_sele = {}
	leaflet_sele_atoms = {}
	for l in ["lower","upper","both"]:
		leaflet_sele[l] = {}
		leaflet_sele_atoms[l] = {}
	
	#check the leaflet selection string is valid
	test_beads = U.selectAtoms(leaflet_sele_string)
	if test_beads.numberOfAtoms() == 0:
		print "Error: invalid selection string '" + str(leaflet_sele_string) + "'"
		print "-> no particles selected."
		sys.exit(1)

	#use LeafletFinder:
	if args.cutoff_leaflet != 'large':
		if args.cutoff_leaflet == 'optimise':
			print " -optimising cutoff..."
			cutoff_value = MDAnalysis.analysis.leaflet.optimize_cutoff(U, leaflet_sele_string)
			L = MDAnalysis.analysis.leaflet.LeafletFinder(U, leaflet_sele_string, cutoff_value[0])
		else:
			L = MDAnalysis.analysis.leaflet.LeafletFinder(U, leaflet_sele_string, args.cutoff_leaflet)
	
		if np.shape(L.groups())[0]<2:
			print "Error: imposssible to identify 2 leaflets."
			sys.exit(1)
		if L.group(0).centerOfGeometry()[2] > L.group(1).centerOfGeometry()[2]:
			leaflet_sele["upper"] = L.group(0)
			leaflet_sele["lower"] = L.group(1)
		else:
			leaflet_sele["upper"] = L.group(1)
			leaflet_sele["lower"] = L.group(0)
		leaflet_sele["both"] = leaflet_sele["lower"] + leaflet_sele["upper"]
		if np.shape(L.groups())[0] == 2:
			print " -found 2 leaflets: ", leaflet_sele["upper"].numberOfResidues(), '(upper) and ', leaflet_sele["lower"].numberOfResidues(), '(lower) lipids'
		else:
			other_lipids=0
			for g in range(2, np.shape(L.groups())[0]):
				other_lipids += L.group(g).numberOfResidues()
			print " -found " + str(np.shape(L.groups())[0]) + " groups: " + str(leaflet_sele["upper"].numberOfResidues()) + "(upper), " + str(leaflet_sele["lower"].numberOfResidues()) + "(lower) and " + str(other_lipids) + " (others) lipids respectively"
	#use cog:
	else:
		leaflet_sele["both"] = U.selectAtoms(leaflet_sele_string)
		tmp_lipids_avg_z = leaflet_sele["both"].centerOfGeometry()[2]
		leaflet_sele["upper"] = leaflet_sele["both"].selectAtoms("prop z > " + str(tmp_lipids_avg_z))
		leaflet_sele["lower"] = leaflet_sele["both"].selectAtoms("prop z < " + str(tmp_lipids_avg_z))
		print " -found 2 leaflets: ", leaflet_sele["upper"].numberOfResidues(), '(upper) and ', leaflet_sele["lower"].numberOfResidues(), '(lower) lipids'
		
	return
def initialise_groups():												#DONE

	colormaps_possible = ['Spectral', 'summer', 'coolwarm', 'pink_r', 'Set1', 'Set2', 'Set3', 'brg_r', 'Dark2', 'hot', 'PuOr_r', 'afmhot_r', 'terrain_r', 'PuBuGn_r', 'RdPu', 'gist_ncar_r', 'gist_yarg_r', 'Dark2_r', 'YlGnBu', 'RdYlBu', 'hot_r', 'gist_rainbow_r', 'gist_stern', 'gnuplot_r', 'cool_r', 'cool', 'gray', 'copper_r', 'Greens_r', 'GnBu', 'gist_ncar', 'spring_r', 'gist_rainbow', 'RdYlBu_r', 'gist_heat_r', 'OrRd_r', 'CMRmap', 'bone', 'gist_stern_r', 'RdYlGn', 'Pastel2_r', 'spring', 'terrain', 'YlOrRd_r', 'Set2_r', 'winter_r', 'PuBu', 'RdGy_r', 'spectral', 'flag_r', 'jet_r', 'RdPu_r', 'Purples_r', 'gist_yarg', 'BuGn', 'Paired_r', 'hsv_r', 'bwr', 'cubehelix', 'YlOrRd', 'Greens', 'PRGn', 'gist_heat', 'spectral_r', 'Paired', 'hsv', 'Oranges_r', 'prism_r', 'Pastel2', 'Pastel1_r', 'Pastel1', 'gray_r', 'PuRd_r', 'Spectral_r', 'gnuplot2_r', 'BuPu', 'YlGnBu_r', 'copper', 'gist_earth_r', 'Set3_r', 'OrRd', 'PuBu_r', 'ocean_r', 'brg', 'gnuplot2', 'jet', 'bone_r', 'gist_earth', 'Oranges', 'RdYlGn_r', 'PiYG', 'CMRmap_r', 'YlGn', 'binary_r', 'gist_gray_r', 'Accent', 'BuPu_r', 'gist_gray', 'flag', 'seismic_r', 'RdBu_r', 'BrBG', 'Reds', 'BuGn_r', 'summer_r', 'GnBu_r', 'BrBG_r', 'Reds_r', 'RdGy', 'PuRd', 'Accent_r', 'Blues', 'Greys', 'autumn', 'cubehelix_r', 'nipy_spectral_r', 'PRGn_r', 'Greys_r', 'pink', 'binary', 'winter', 'gnuplot', 'RdBu', 'prism', 'YlOrBr', 'coolwarm_r', 'rainbow_r', 'rainbow', 'PiYG_r', 'YlGn_r', 'Blues_r', 'YlOrBr_r', 'seismic', 'Purples', 'bwr_r', 'autumn_r', 'ocean', 'Set1_r', 'PuOr', 'PuBuGn', 'nipy_spectral', 'afmhot']

	global colours_groups
	global colours_groups_list
	global groups_labels
	global groups_number
	global groups_boundaries
	global groups_sizes_dict

	groups_labels = {}
	groups_boundaries = {}
	groups_sizes_dict = {}
	colours_groups_dict = {}
	colours_groups_list = []
	colours_groups_map = "custom"
	
	#read file
	print "\nReading cluster groups definition file..."
	with open(args.cluster_groups_file) as f:
		lines = f.readlines()
	groups_number = len(lines)
	for g_index in range(0,groups_number):
		line = lines[g_index]
		if line[-1] == "\n":
			line = line[:-1]
		l_content = line.split(',')
		#check format
		if len(l_content) != 3:
			print "Error: the format of line " + str(g_index+1) + " should be 'min,max,colour' (see bilayer_perturbations --help, note 7)."
			print "->", line
			sys.exit(1)
		tmp_beg = int(l_content[0])
		tmp_end = l_content[1]
		colours_groups_dict[g_index] = l_content[2]						#attribute colour to group
		if l_content[2] not in colormaps_possible:	
			colours_groups_list.append(groups_colors_value[g_index])
		if tmp_end == "max":
			tmp_end = 100000										#put a stupidly big size to cap the open ended group
		else:
			tmp_end = int(tmp_end)
		groups_boundaries[g_index] = [tmp_beg,tmp_end]
		
	#display results
	print " -found " + str(groups_number) + " cluster groups:"
	for g_index in range(0,groups_number):
		if groups_boundaries[g_index][1] == 100000:
			print "   g" + str(g_index) + " = " + str(groups_boundaries[g_index][0]) + "+, " + str(colours_groups_dict[g_index])
		else:
			print "   g" + str(g_index) + " = " + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + ", " + str(colours_groups_dict[g_index])

	#check for boundaries overlapping
	prev_beg = groups_boundaries[0][0]
	prev_end = groups_boundaries[0][1]
	if prev_end < prev_beg:
		print "Error: the max size is smaller than the min size for specified cluster groups " + str(groups_boundaries[0]) + "."
		sys.exit(1)
	for g_index in range(1,groups_number):
		if groups_boundaries[g_index][1] < groups_boundaries[g_index][0]:
			print "Error: the max size is smaller than the min size for group " + str(g_index) + "(" + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + ")."
			sys.exit(1)
		if groups_boundaries[g_index][0] <= prev_end:
			print "Error: specified cluster groups " + str(prev_beg) + "-" + str(prev_end) + " and " + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + " overlap or are not in increasing order (boundaries are inclusive, see note 7 in bilayer_perturbations --help)."
			sys.exit(1)
		prev_beg = groups_boundaries[g_index][0]
		prev_end = groups_boundaries[g_index][1]
	
	#create equivalency table between groups and sizes
	for g_index in range(0,groups_number):
		bb = groups_boundaries[g_index]
		tmp_beg = bb[0]
		tmp_end = bb[1]
		for tmp_size in range(tmp_beg, tmp_end+1):
			groups_sizes_dict[tmp_size] = g_index
	for tmp_size in list(set(range(1,max(groups_sizes_dict.keys()))) - set(groups_sizes_dict.keys())): 		#this handles potentially unaccounted for sizes up to the maximum specified by the user
		groups_sizes_dict[tmp_size] = groups_number
	if max(groups_sizes_dict.keys())!= 100000:															    #this handles potentially unaccounted for sizes above the maximum specified by the user (in case it's not an open group)
		for tmp_size in range(max(groups_sizes_dict.keys())+1,100001):
			groups_sizes_dict[tmp_size] = groups_number

	#check whether a colour map was specified
	if groups_number > 1 and len(np.unique(colours_groups_dict.values())) == 1:
		if np.unique(colours_groups_dict.values())[0] in colormaps_possible:
			colours_groups_map = np.unique(colours_groups_dict.values())[0]
		else:
			print "Error: either the same color was specified for all groups or the color map '" + str(np.unique(colours_groups_dict.values())[0]) + "' is not valid."
			sys.exit(1)

	#generate colours from colour map if necessary
	if colours_groups_map != "custom":
		tmp_cmap = cm.get_cmap(colors_groups_map)
		groups_colors_value = tmp_cmap(np.linspace(0, 1, groups_number))
		for g_index in range(0, groups_number):
			colours_groups_dict[g_index] = groups_colors_value[g_index]
			colours_groups_list.append(groups_colors_value[g_index])
			
	#create label for each group
	for g_index in range(0, groups_number):
		if groups_boundaries[g_index][1] == 100000:
			groups_labels[g_index] = str(groups_boundaries[g_index][0]) + "+"
		elif groups_boundaries[g_index][0] == groups_boundaries[g_index][1]:
			groups_labels[g_index] = str(groups_boundaries[g_index][0])
		else:
			groups_labels[g_index] = str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1])

	#add colours for the group "other", "lower" and "upper"
	groups_labels[groups_number] = "other"
	groups_labels[groups_number+1] = "lower"
	groups_labels[groups_number+2] = "upper"
	colours_groups_dict[groups_number] = colour_group_other
		
	return

#=========================================================================================
# data structures
#=========================================================================================

def data_struct_time():													#DONE

	global frames_nb
	global frames_time
	frames_nb = np.zeros(nb_frames_to_process)
	frames_time = np.zeros(nb_frames_to_process)

	return
def data_struct_stats():												#DONE
	
	#NB:
	# 1. relative distribution of protein population at each frame
	# 2. for groups groups_number = "other", -1 = lower, groups_number + 1 = upper
	
	global cluster_sizes_nb
	global cluster_sizes_pc
	global cluster_biggest_nb
	global cluster_biggest_pc
	global cluster_biggest_size
	
	cluster_sizes_nb = {}
	cluster_sizes_pc = {}
	cluster_biggest_nb = np.zeros(nb_frames_to_process)
	cluster_biggest_pc = np.zeros(nb_frames_to_process)
	cluster_biggest_size = np.zeros(nb_frames_to_process)
	if args.cluster_groups_file != "no":
		global cluster_groups_nb
		global cluster_groups_pc
		cluster_groups_nb = {g_index: np.zeros(nb_frames_to_process) for g_index in range(-1,groups_number + 2)}
		cluster_groups_pc = {g_index: np.zeros(nb_frames_to_process) for g_index in range(-1,groups_number + 2)}

		if args.xtcfilename != "no":
			global proteins_groups_stability
			proteins_groups_stability = {}
	
	return
def data_struct_proteins():												#DONE
	
	#NB:
	# 1. size of the cluster size/index of the size group each protein is involved in at each frame
	# 2. for sizes -1 = lower, proteins_nb = upper
	# 3. for groups groups_number = "other", - 1 = lower, groups_number + 1 = upper
	
	global proteins_cluster_status_sizes

	proteins_cluster_status_sizes = np.zeros((proteins_nb, nb_frames_to_process))
	if args.cluster_groups_file != "no":
		global proteins_cluster_status_groups
		proteins_cluster_status_groups = np.zeros((proteins_nb, nb_frames_to_process))
	
	return
		
#=========================================================================================
# core functions
#=========================================================================================

def get_distances(box_dim):												#DONE
	
	#method: use minimum distance between proteins
	#---------------------------------------------
	if args.m_algorithm == "min":
		#pre-process: get protein coordinates
		tmp_proteins_coords = np.zeros((proteins_nb, nb_atom_per_protein, 3))
		for p_index in range(0, proteins_nb):
			tmp_proteins_coords[p_index,:] = proteins_sele[p_index].coordinates()

		#store min distance between each proteins
		dist_matrix = 100000 * np.ones((proteins_nb,proteins_nb))
		for n in range(proteins_nb,1,-1):
			dist_matrix[proteins_nb-n,proteins_nb-n+1:proteins_nb] = map(lambda pp: np.min(MDAnalysis.analysis.distances.distance_array(np.float32(tmp_proteins_coords[proteins_nb-n,:]), np.float32(tmp_proteins_coords[pp,:]), box_dim)), range(proteins_nb-n+1,proteins_nb))
			dist_matrix[proteins_nb-n+1:proteins_nb,proteins_nb-n] = dist_matrix[proteins_nb-n,proteins_nb-n+1:proteins_nb]
											
	#method: use distance between cog
	#--------------------------------
	else:
		tmp_proteins_cogs = np.asarray(map(lambda p_index: calculate_cog(proteins_sele[p_index], box_dim), range(0,proteins_nb)))
		dist_matrix = MDAnalysis.analysis.distances.distance_array(np.float32(tmp_proteins_cogs), np.float32(tmp_proteins_cogs), box_dim)

	return dist_matrix
def calculate_cog(tmp_coords, box_dim):									#DONE
	
	#this method allows to take pbc into account when calculcating the center of geometry 
	#see: http://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
	
	cog_coord = np.zeros(3)
	tmp_nb_atoms = np.shape(tmp_coords)[0]
	
	for n in range(0,3):
		tet = tmp_coords[:,n] * 2 * math.pi / float(box_dim[n])
		xsi = np.cos(tet)
		zet = np.sin(tet)
		tet_avg = math.atan2(-np.average(zet),-np.average(xsi)) + math.pi
		cog_coord[n] = tet_avg * box_dim[n] / float(2*math.pi)
	
	return cog_coord
def detect_clusters_connectivity(dist, box_dim):						#DONE
	
	#use networkx algorithm
	connected=(dist<args.cutoff_connect)
	network=nx.Graph(connected)
	groups=nx.connected_components(network)
	
	return groups
def detect_clusters_density(dist, box_dim):								#DONE
	
	#run DBSCAN algorithm
	dbscan_output = DBSCAN(eps = args.dbscan_dist, metric = 'precomputed', min_samples = args.dbscan_nb).fit(dist)

	#build 'groups' structure i.e. a list whose element are all the clusters identified
	groups = []
	for c_lab in np.unique(dbscan_output.labels_):
		tmp_pos = np.argwhere(dbscan_output.labels_ == c_lab)
		if c_lab == -1:
			groups += map(lambda p:p[0] , tmp_pos)
		else:
			groups.append(map(lambda p:p[0] , tmp_pos))

	return groups
def process_clusters(clusters, f_index):								#DONE

	global vmd_cluster_size
	
	#case: store cluster size only
	#-----------------------------
	if args.cluster_groups_file == "no":
		for cluster in clusters:
			c_size = np.size(cluster)
			proteins_cluster_status_sizes[cluster, f_index] = c_size

	#case: store cluster size and group size
	#---------------------------------------
	else:
		for cluster in clusters:
			c_size = np.size(cluster)
			g_index = groups_sizes_dict[c_size]
			proteins_cluster_status_sizes[cluster, f_index] = c_size
			proteins_cluster_status_groups[cluster, f_index] = g_index

	#create annotation line for current frame
	#----------------------------------------
	if args.buffer_size != -1:
		#size
		tmp_size = str(frames_nb[f_index])
		for p_index in range(0,proteins_nb):
			tmp_size += "." + str(proteins_cluster_status_sizes[p_index,f_index])
		vmd_cluster_size += tmp_size + "\n"
		if vmd_counter == args.buffer_size:
			with open(output_xtc_annotate_cluster_size, 'a') as f:
				f.write(vmd_cluster_size)
				vmd_cluster_size = ""

		#groups
		if args.cluster_groups_file!="no":
			global vmd_cluster_group
			tmp_group = str(frames_nb[f_index])
			for p_index in range(0,proteins_nb):
				tmp_group += "." + str(proteins_cluster_status_groups[p_index,f_index])
			vmd_cluster_group += tmp_group + "\n"
			if vmd_counter == args.buffer_size:
				with open(output_xtc_annotate_cluster_group, 'a') as f:
					f.write(vmd_cluster_group)
					vmd_cluster_group = ""

	return
def process_clusters_TM(clusters, f_index, box_dim):					#DONE

	global vmd_cluster_size
		
	tmp_lip_coords = {l: leaflet_sele[l].coordinates() for l in ["lower","upper"]}
	
	#case: store cluster size only
	#-----------------------------
	if args.cluster_groups_file == "no":			
		for cluster in clusters:
			c_size = np.size(cluster)

			#create selection for current cluster
			c_sele = MDAnalysis.core.AtomGroup.AtomGroup([])	
			for p_index in cluster:
				c_sele += proteins_sele[p_index]
			tmp_c_sele_coordinates = c_sele.coordinates()
			#find closest PO4 particles for each particles of clusters, if all are in the same leaflet then it's surfacic [NB: this is done at the CLUSTER level (the same criteria at the protein level would probably fail)]
			dist_min_lower = np.min(MDAnalysis.analysis.distances.distance_array(tmp_c_sele_coordinates, tmp_lip_coords["lower"], box_dim), axis = 1)
			dist_min_upper = np.min(MDAnalysis.analysis.distances.distance_array(tmp_c_sele_coordinates, tmp_lip_coords["upper"], box_dim), axis = 1)
			dist = dist_min_upper - dist_min_lower
			   
			#case: interfacial lower
			if np.size(dist[dist>0]) == np.size(dist):
				proteins_cluster_status_sizes[cluster, f_index] = -1
			#case: interfacial upper
			elif np.size(dist[dist>0]) == 0:
				proteins_cluster_status_sizes[cluster, f_index] = proteins_nb+1
			#case: TM
			else:
				proteins_cluster_status_sizes[cluster, f_index] = c_size				

	#case: store cluster size and group size
	#---------------------------------------
	else:
		for cluster in clusters:
			c_size = np.size(cluster)
			g_index = groups_sizes_dict[c_size]

			#create selection for current cluster
			c_sele = MDAnalysis.core.AtomGroup.AtomGroup([])	
			for p_index in cluster:
				c_sele += proteins_sele[p_index]
			tmp_c_sele_coordinates = c_sele.coordinates()
			#find closest PO4 particles for each particles of clusters, if all are in the same leaflet then it's surfacic [NB: this is done at the CLUSTER level (the same criteria at the protein level would probably fail)]
			dist_min_lower = np.min(MDAnalysis.analysis.distances.distance_array(tmp_c_sele_coordinates, tmp_lip_coords["lower"], box_dim), axis = 1)
			dist_min_upper = np.min(MDAnalysis.analysis.distances.distance_array(tmp_c_sele_coordinates, tmp_lip_coords["upper"], box_dim), axis = 1)
			dist = dist_min_upper - dist_min_lower
			   
			#case: interfacial lower
			if np.size(dist[dist>0]) == np.size(dist):
				proteins_cluster_status_sizes[cluster, f_index] = -1
				proteins_cluster_status_groups[cluster, f_index] = -1
			#case: interfacial upper
			elif np.size(dist[dist>0]) == 0:
				proteins_cluster_status_sizes[cluster, f_index] = proteins_nb+1
				proteins_cluster_status_groups[cluster, f_index] = groups_number + 1
			#case: TM
			else:
				proteins_cluster_status_sizes[cluster, f_index] = c_size				
				proteins_cluster_status_groups[cluster, f_index] = g_index
	
	#create annotation line for current frame
	#----------------------------------------
	if args.buffer_size != -1:
		#size
		tmp_size = str(frames_nb[f_index])
		for p_index in range(0,proteins_nb):
			tmp_size += "." + str(proteins_cluster_status_sizes[p_index,f_index])
		vmd_cluster_size += tmp_size + "\n"
		if vmd_counter == args.buffer_size:
			with open(output_xtc_annotate_cluster_size, 'a') as f:
				f.write(vmd_cluster_size)
				vmd_cluster_size = ""

		#groups
		if args.cluster_groups_file!="no":
			global vmd_cluster_group
			tmp_group = str(frames_nb[f_index])
			for p_index in range(0,proteins_nb):
				tmp_group += "." + str(proteins_cluster_status_groups[p_index,f_index])
			vmd_cluster_group += tmp_group + "\n"
			if vmd_counter == args.buffer_size:
				with open(output_xtc_annotate_cluster_group, 'a') as f:
					f.write(vmd_cluster_group)
					vmd_cluster_group = ""
	
	return
def get_sizes_sampled():												#DONE
		
	global cluster_sizes_sampled
	global cluster_TM_sizes_sampled
	global cluster_sizes_nb
	global cluster_sizes_pc
	cluster_sizes_sampled = sorted(np.unique(proteins_cluster_status_sizes))
	cluster_TM_sizes_sampled = [c_size for c_size in cluster_sizes_sampled if (c_size != -1 and c_size != proteins_nb+1)]
	cluster_sizes_nb = {c_size: np.zeros(nb_frames_to_process) for c_size in cluster_sizes_sampled + [-1,proteins_nb+1]}
	cluster_sizes_pc = {c_size: np.zeros(nb_frames_to_process) for c_size in cluster_sizes_sampled + [-1,proteins_nb+1]}

	if args.cluster_groups_file!="no":
		global cluster_groups_sampled
		cluster_groups_sampled = sorted(np.unique(proteins_cluster_status_groups))
		
	return
def update_color_dict():												#DONE
	
	global colours_sizes_nb
	global colours_sizes_dict
	global colours_sizes_list
		
	#colours: sizes (colours from jet colormap)
	#-------
	colours_sizes_nb = colours_sizes_range[1]-colours_sizes_range[0]+1
	colours_sizes_dict = {}
	colours_sizes_list = []

	tmp_cmap = cm.get_cmap('jet')
	colours_sizes_value = tmp_cmap(np.linspace(0, 1, colours_sizes_range[1]-colours_sizes_range[0]+1))
	for c_index in range(0, colours_sizes_range[1]-colours_sizes_range[0]+1):
		c_size = colours_sizes_range[0] + c_index
		colours_sizes_dict[c_size] = colours_sizes_value[c_index]
		colours_sizes_list.append(colours_sizes_value[c_index])
	#interfacial proteins on the lower leaflet (prepend -> bottom of colour bar)
	colours_sizes_dict[-1] = colour_leaflet_lower
	colours_sizes_list.insert(0, colour_leaflet_lower)
	#interfacial proteins on the upper leaflet (append -> top of colour bar)
	colours_sizes_dict[proteins_nb+1] = colour_leaflet_upper
	colours_sizes_list.append(colour_leaflet_upper)

	#colours: groups
	#-------
	if args.cluster_groups_file != "no":
		global colours_groups_nb
		colours_groups_nb = len(cluster_groups_sampled)
		#proteins in cluster of a different size		
		if groups_number in cluster_groups_sampled:
			colours_groups_list.append(colour_group_other)
		#interfacial proteins on the lower leaflet (prepend -> bottom of colour bar)
		colours_groups_list.insert(0, colour_leaflet_lower)
		#interfacial proteins on the upper leaflet (append -> top of colour bar)
		colours_groups_list.append(colour_leaflet_upper)

	return
def calculate_statistics():												#DONE

	#nb and corresponding % of each size / group
	#-------------------------------------------
	for f_index in range(0,nb_frames_to_process):
		
		#retrieve cluster status of all proteins for current frame
		tmp_sizes = list(proteins_cluster_status_sizes[:,f_index])
		
		#initialise values for biggest cluster stat (nb,% and size)
		tmp_max_nb = 0
		tmp_max_pc = 0
		tmp_max_size = float("-inf")
		
		#store current frame statistics for each size ever sampled
		for c_size in cluster_sizes_sampled:
			tmp_nb = tmp_sizes.count(c_size)
			tmp_pc = tmp_nb / float(proteins_nb) * 100
			if c_size != -1 and c_size != proteins_nb+1:					#take into account the cluster size except for interfacial peptides
				tmp_nb = int(tmp_nb / float(c_size))
				
			cluster_sizes_nb[c_size][f_index] = tmp_nb
			cluster_sizes_pc[c_size][f_index] = tmp_pc
			if tmp_nb > 0 and c_size > tmp_max_size and c_size != proteins_nb+1:
				tmp_max_nb = tmp_nb
				tmp_max_pc = tmp_pc
				tmp_max_size = c_size
			if args.cluster_groups_file != "no":
				g_index = groups_sizes_dict[c_size]
				cluster_groups_nb[g_index][f_index][g_index] += tmp_nb
				cluster_groups_pc[g_index][f_index][g_index] += tmp_pc
		
		#store biggest cluster stats
		cluster_biggest_nb[f_index] = tmp_max_nb
		cluster_biggest_pc[f_index] = tmp_max_pc		
		cluster_biggest_size[f_index] = tmp_max_size			

	#longest stability of each size group
	#------------------------------------
	if args.xtcfilename != "no" and args.cluster_groups_file != "no":
		for g_index in cluster_TM_groups_sampled:
			#find max stability of current group index for each protein
			tmp_proteins_stability = {}
			for p_index in range(0,proteins_nb):
				if g_index in proteins_cluster_group_mat[p_index,:]:
					tmp_proteins_stability[p_index] = max(len(list(v)) for g,v in itertools.groupby(proteins_cluster_status_groups[p_index,:], lambda x: x == g_index) if g)
				else:
					tmp_proteins_stability[p_index] = 0

			#store maximum stability
			proteins_groups_stability[g_index] = max(tmp_proteins_stability.values())

	return
def rolling_avg(loc_list):												#DONE
	
	loc_arr=np.asarray(loc_list)
	shape=(loc_arr.shape[-1]-args.nb_smoothing+1,args.nb_smoothing)
	strides=(loc_arr.strides[-1],loc_arr.strides[-1])   	
	return np.average(np.lib.stride_tricks.as_strided(loc_arr, shape=shape, strides=strides), -1)
def smooth_data():														#DONE

	global frames_time_smoothed
	global cluster_biggest_nb_smoothed
	global cluster_biggest_pc_smoothed	
	global cluster_biggest_size_smoothed
		
	#time
	frames_time_smoothed = rolling_avg(frames_time)
	
	#biggest cluster size
	cluster_biggest_nb_smoothed = rolling_avg(cluster_biggest_nb)
	cluster_biggest_pc_smoothed = rolling_avg(cluster_biggest_pc)
	cluster_biggest_size_smoothed = rolling_avg(cluster_biggest_size)
	
	#size
	for c_size in cluster_sizes_sampled:
		cluster_sizes_nb_smoothed[c_size] = rolling_avg(cluster_sizes_nb[c_size])
		cluster_sizes_pc_smoothed[c_size] = rolling_avg(cluster_sizes_pc[c_size])
	
	#groups
	if args.cluster_groups_file!="no":
		for g_index in cluster_groups_sampled:
			cluster_groups_nb_smoothed[g_index] = rolling_avg(cluster_groups_nb[g_index])
			cluster_groups_pc_smoothed[g_index] = rolling_avg(cluster_groups_pc[g_index])
		
	return

#=========================================================================================
# outputs
#=========================================================================================

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
	output_stat.write("Warning: a single cluster size (" + str(cluster_sizes_sampled[0]) + ") was detected throughout the trajectory. Check the -m, -c, -r or -n options (see cluster_prot -h).")
	output_stat.close()
	
	return

#sizes
def write_xvg_biggest():												#DONE
	filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_3_biggest/xvg/1_2_clusterprot_biggest.txt'
	filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_3_biggest/xvg/1_2_clusterprot_biggest.xvg'
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
	for f_index in range(0,nb_frames_to_process):
		results = str(frames_time[f_index]) + "	" + str(cluster_biggest_size[f_index]) + "	" + str(round(cluster_biggest_pc[f_index],2)) + "	" + str(cluster_biggest_nb[f_index])
		output_xvg.write(results + "\n")
	output_xvg.close()
	return
def graph_xvg_biggest():												#DONE
	#create filenames
	#----------------
	filename_png = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_3_biggest/png/1_2_clusterprot_biggest.png'
	filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_3_biggest/1_2_clusterprot_biggest.svg'

	#create figure
	#-------------
	fig=plt.figure(figsize=(8, 6.2))
	fig.suptitle("Evolution of the size of the biggest protein cluster")

	#plot data: %
	#------------
	ax1 = fig.add_subplot(211)
	p_upper = {}
	p_upper["size"] = plt.plot(frames_time, cluster_biggest_size, color = 'k', linewidth = 2.0, label = "size")
	fontP.set_size("small")
	ax1.legend(prop=fontP)
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('cluster size', fontsize="small")

	#plot data: nb #TO DO: make 2 y axis for the bottom bit, to include nb of clusters
	#-------------
	ax2 = fig.add_subplot(212)
	p_lower = {}
	p_lower["pc"] = plt.plot(frames_time, cluster_biggest_pc, color = 'c', linewidth = 2.0, label = "% of proteins")
	#p_lower["nb"]=plt.plot(time_sorted, cluster_biggest_nb_sorted, color='r', linewidth=2.0, label="nb of clusters")
	fontP.set_size("small")
	ax2.legend(prop=fontP)
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('% of proteins', fontsize="small")

	#save figure
	#-----------
	ax1.set_ylim(0, max(cluster_biggest_size)+2)
	ax2.set_ylim(0, 100)
	ax1.spines['top'].set_visible(False)
	ax1.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.spines['right'].set_visible(False)
	ax1.xaxis.set_ticks_position('bottom')
	ax1.yaxis.set_ticks_position('left')
	ax2.xaxis.set_ticks_position('bottom')
	ax2.yaxis.set_ticks_position('left')
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
def write_xvg_cluster_biggest_smoothed():								#DONE
	filename_txt=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_3_biggest/xvg/1_3_clusterprot_cluster_biggest_smooth.txt'
	filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_3_biggest/xvg/1_3_clusterprot_cluster_biggest_smooth.xvg'
	output_txt = open(filename_txt, 'w')
	output_txt.write("@[lipid tail order parameters statistics - written by order_param v" + str(version_nb) + "]\n")
	output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 1_3_clusterprot_cluster_biggest_smooth.xvg.\n")
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
	output_txt.write("1_3_clusterprot_cluster_biggest_smooth.xvg,1,size,k\n")
	output_txt.write("1_3_clusterprot_cluster_biggest_smooth.xvg,2,%,c\n")
	output_txt.write("1_3_clusterprot_cluster_biggest_smooth.xvg,3,nb,r\n")
	for f_index in range(0, frames_time_smoothed):
		results = str(frames_time_smoothed[f_index]) + "	" + str(cluster_biggest_size_smoothed[f_index]) + "	" + str(round(cluster_biggest_pc_smoothed[f_index],2)) + "	" + str(cluster_biggest_nb_smoothed[f_index])
		output_xvg.write(results + "\n")
	output_xvg.close()
	return
def graph_xvg_cluster_biggest_smoothed():								#DONE
	#create filenames
	#----------------
	filename_png = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_3_biggest/png/1_3_clusterprot_cluster_biggest_smoothed.png'
	filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_3_biggest/1_3_clusterprot_cluster_biggest_smoothed.svg'

	#create figure
	#-------------
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("Evolution of the size of the biggest protein cluster")
		
	#plot data: %
	#------------
	ax1 = fig.add_subplot(211)
	p_upper = {}
	p_upper["size"] = plt.plot(frames_time_smoothed, cluster_biggest_size_smoothed, color = 'k', linewidth = 2.0, label = "size")
	fontP.set_size("small")
	ax1.legend(prop=fontP)
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('cluster size', fontsize="small")

	#plot data: nb #TO DO: make 2 y axis for the bottom bit, to include nb of clusters
	#-------------
	ax2 = fig.add_subplot(212)
	p_lower = {}
	p_lower["pc"] = plt.plot(frames_time_smoothed, cluster_biggest_pc_smoothed, color = 'c', linewidth = 2.0, label = "% of proteins")
	#p_lower["nb"] = plt.plot(time_smoothed, cluster_biggest_nb_smoothed, color='r', linewidth=2.0, label="nb of clusters")
	fontP.set_size("small")
	ax2.legend(prop=fontP)
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('% of proteins', fontsize="small")

	#save figure
	#-----------
	ax1.set_ylim(0, np.floor(max(cluster_biggest_size)+2))
	ax2.set_ylim(0, 100)
	ax1.spines['top'].set_visible(False)
	ax1.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.spines['right'].set_visible(False)
	ax1.xaxis.set_ticks_position('bottom')
	ax1.yaxis.set_ticks_position('left')
	ax2.xaxis.set_ticks_position('bottom')
	ax2.yaxis.set_ticks_position('left')
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
def write_xvg_sizes_TM():												#DONE

	filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_2_plots_1D/xvg/1_2_clusterprot_1D_TM.txt'
	filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_2_plots_1D/xvg/1_2_clusterprot_1D_TM.xvg'
	output_txt = open(filename_txt, 'w')
	output_txt.write("# [protein aggregation statistics - written by cluster_prot v" + str(version_nb) + "]\n")
	output_txt.write("# Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 1_2_clusterprot_1D_TM.xvg.\n")
	output_xvg = open(filename_xvg, 'w')
	output_xvg.write("# [protein aggregation statistics - written by cluster_prot v" + str(version_nb) + "]\n")
	output_xvg.write("# - proteins detected: " + str(proteins_nb) + "\n")
	output_xvg.write("# - algorithm chosen: " + str(args.m_algorithm) + "\n")
	if args.m_algorithm == "density":
		output_xvg.write("# - search radius: " + str(args.dbscan_dist) + "\n")
		output_xvg.write("# - nb of neighbours: " + str(args.dbscan_nb) + "\n")
	else:
		output_xvg.write("# - cutoff for contact: " + str(args.cutoff_connect) + "\n")
	output_xvg.write("@ title \"Evolution of protein aggregation\"\n")
	output_xvg.write("@ xaxis  label \"time (ns)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	tmp_minus = 0
	output_xvg.write("@ legend length " + str(len(cluster_TM_sizes_sampled)*2) + "\n")
	#write caption: %
	for c_index in range(0,len(cluster_TM_sizes_sampled)):
		c_size = cluster_TM_sizes_sampled[c_index]
		output_xvg.write("@ s" + str(c_index) + " legend \"% " + str(c_size) + "\"\n")
		output_txt.write("1_2_clusterprot_1D.xvg," + str(c_index+1) + ",% " + str(c_size) + "," + mcolors.rgb2hex(mcolorconv.to_rgb(colours_sizes_dict[c_size])) + "\n")
	#write caption: nb
	for c_index in range(0,len(cluster_TM_sizes_sampled)):
		c_size = cluster_TM_sizes_sampled[c_index]
		output_xvg.write("@ s" + str(len(cluster_TM_sizes_sampled) + c_index) + " legend \"nb " + str(c_size) + "\"\n")
		output_txt.write("1_2_clusterprot_1D.xvg," + str(len(cluster_sizes_sampled) + c_index + 1) + ",nb " + str(c_size) + "," + mcolors.rgb2hex(mcolorconv.to_rgb(colours_sizes_dict[c_size])) + "\n")
	output_txt.close()
	#write results
	for f_index in range(0,nb_frames_to_process):
		results = str(frames_time[f_index])
		for c_size in cluster_sizes_sampled:
			results+="	" + str(round(cluster_sizes_pc[c_size][f_index],2))
		for c_size in cluster_sizes_sampled:
			results+="	" + str(round(cluster_sizes_nb[c_size][f_index],2))
		output_xvg.write(results + "\n")
	output_xvg.close()
	
	return
def graph_xvg_sizes_TM():												#DONE
	
	#create filenames
	#----------------
	filename_png = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_2_plots_1D/png/1_2_clusterprot_1D_TM.png'
	filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_2_plots_1D/1_2_clusterprot_1D_TM.svg'

	#create figure
	#-------------
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("Evolution of proteins distribution")
		
	#plot data: %
	#------------
	ax1 = fig.add_subplot(211)
	p_upper = {}
	for c_size in cluster_sizes_sampled:
		if c_size != -1 and c_size != proteins_nb+1:
			p_upper[c_size] = plt.plot(frames_time, cluster_sizes_pc[c_size], color = colours_sizes_dict[c_size], linewidth = 2.0, label = str(int(c_size)))
	fontP.set_size("small")
	ax1.legend(prop=fontP)
	plt.title("%", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('% of proteins', fontsize="small")

	#plot data: nb
	#-------------
	ax2 = fig.add_subplot(212)
	p_lower = {}
	for c_size in cluster_sizes_sampled:
		if c_size != -1 and c_size != proteins_nb+1:
			p_lower[c_size] = plt.plot(frames_time, cluster_sizes_nb[c_size], color = colours_sizes_dict[c_size], linewidth = 2.0, label=str(int(c_size)))
	fontP.set_size("small")
	ax2.legend(prop=fontP)
	plt.title("nb", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('nb of clusters', fontsize="small")

	#save figure
	#-----------
	ax1.set_ylim(0, 100)
	ax2.set_ylim(0, np.max(cluster_sizes_nb.values())+1)
	ax1.spines['top'].set_visible(False)
	ax1.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.spines['right'].set_visible(False)
	ax1.xaxis.set_ticks_position('bottom')
	ax1.yaxis.set_ticks_position('left')
	ax2.xaxis.set_ticks_position('bottom')
	ax2.yaxis.set_ticks_position('left')
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
def write_xvg_sizes_interfacial():										#DONE

	filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_2_plots_1D/xvg/1_2_clusterprot_1D_interfacial.txt'
	filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_2_plots_1D/xvg/1_2_clusterprot_1D_interfacial.xvg'
	output_txt = open(filename_txt, 'w')
	output_txt.write("# [protein aggregation statistics - written by cluster_prot v" + str(version_nb) + "]\n")
	output_txt.write("# Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 1_2_clusterprot_1D_interfacial.xvg.\n")
	output_xvg = open(filename_xvg, 'w')
	output_xvg.write("# [protein aggregation statistics - written by cluster_prot v" + str(version_nb) + "]\n")
	output_xvg.write("# - proteins detected: " + str(proteins_nb) + "\n")
	output_xvg.write("# - algorithm chosen: " + str(args.m_algorithm) + "\n")
	if args.m_algorithm == "density":
		output_xvg.write("# - search radius: " + str(args.dbscan_dist) + "\n")
		output_xvg.write("# - nb of neighbours: " + str(args.dbscan_nb) + "\n")
	else:
		output_xvg.write("# - cutoff for contact: " + str(args.cutoff_connect) + "\n")
	output_xvg.write("@ title \"Evolution of the number of interfacial proteins\"\n")
	output_xvg.write("@ xaxis  label \"time (ns)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length 4\n")
	output_xvg.write("@ s0 legend \"lower (%)\"\n")
	output_xvg.write("@ s1 legend \"upper (%)\"\n")
	output_xvg.write("@ s2 legend \"lower (nb)\"\n")
	output_xvg.write("@ s3 legend \"upper (nb)\"\n")
	output_txt.write("1_2_clusterprot_1D_interfacial.xvg,1, lower (%)," + mcolors.rgb2hex(mcolorconv.to_rgb(colours_sizes_dict[-1])) + "\n")
	output_txt.write("1_2_clusterprot_1D_interfacial.xvg,2, upper (%)," + mcolors.rgb2hex(mcolorconv.to_rgb(colours_sizes_dict[proteins_nb+1])) + "\n")
	output_txt.write("1_2_clusterprot_1D_interfacial.xvg,3, lower (nb)," + mcolors.rgb2hex(mcolorconv.to_rgb(colours_sizes_dict[-1])) + "\n")
	output_txt.write("1_2_clusterprot_1D_interfacial.xvg,4, upper (nb)," + mcolors.rgb2hex(mcolorconv.to_rgb(colours_sizes_dict[proteins_nb+1])) + "\n")
	output_txt.close()
	for f_index in range(0,nb_frames_to_process):
		results = str(frames_time[f_index]) + "	" + str(round(cluster_sizes_pc[-1][f_index],2)) + "	" + str(round(cluster_sizes_pc[proteins_nb+1][f_index],2)) + "	" + str(round(cluster_sizes_pc[-1][f_index],2)) + "	" + str(round(cluster_sizes_pc[proteins_nb+1][f_index],2))
		output_xvg.write(results + "\n")
	output_xvg.close()
	
	return
def graph_xvg_sizes_interfacial():										#DONE
	
	#create filenames
	#----------------
	filename_png = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_2_plots_1D/png/1_2_clusterprot_1D_interfacial.png'
	filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_2_plots_1D/1_2_clusterprot_1D_interfacial.svg'

	#create figure
	#-------------
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("Evolution of proteins distribution")
		
	#plot data: %
	#------------
	ax1 = fig.add_subplot(211)
	p_upper = {}
	p_upper[proteins_nb+1] = plt.plot(frames_time, cluster_sizes_pc[-1], color = colours_sizes_dict[proteins_nb+1], linewidth = 2.0, label = "upper leaflet")
	p_upper[-1] = plt.plot(frames_time, cluster_sizes_pc[-1], color = colours_sizes_dict[-1], linewidth = 2.0, label = "lower leaflet")
	fontP.set_size("small")
	ax1.legend(prop=fontP)
	plt.title("%", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('% of proteins', fontsize="small")

	#plot data: nb
	#-------------
	ax2 = fig.add_subplot(212)
	p_lower = {}
	p_lower[proteins_nb+1] = plt.plot(frames_time, cluster_sizes_nb[-1], color = colours_sizes_dict[proteins_nb+1], linewidth = 2.0, label = "upper leaflet")
	p_lower[-1] = plt.plot(frames_time, cluster_sizes_nb[-1], color = colours_sizes_dict[-1], linewidth = 2.0, label = "lower leaflet")
	fontP.set_size("small")
	ax2.legend(prop=fontP)
	plt.title("nb", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('nb of clusters', fontsize="small")

	#save figure
	#-----------
	ax1.set_ylim(0, 100)
	ax2.set_ylim(0, np.max(cluster_sizes_nb.values())+1)
	ax1.spines['top'].set_visible(False)
	ax1.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.spines['right'].set_visible(False)
	ax1.xaxis.set_ticks_position('bottom')
	ax1.yaxis.set_ticks_position('left')
	ax2.xaxis.set_ticks_position('bottom')
	ax2.yaxis.set_ticks_position('left')
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

def write_xvg_sizes_smoothed():											#DONE
	filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_4_plots_1D_smoothed/xvg/1_4_clusterprot_1D_smoothed.txt'
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
	output_xvg.write("@ legend length " + str(len(cluster_sizes_sampled)*2) + "\n")
	#write caption: %
	for c_index in range(0,len(cluster_sizes_sampled)):
		c_size = cluster_sizes_sampled[c_index]
		output_xvg.write("@ s" + str(c_index) + " legend \"% " + str(c_size) + "\"\n")
		output_txt.write("1_4_clusterprot_1D_smoothed.xvg," + str(c_index+1) + ",% " + str(c_size) + "," + mcolors.rgb2hex(mcolorconv.to_rgb(colours_sizes_dict[c_size])) + "\n")
	#write caption: nb
	for c_index in range(0,len(cluster_sizes_sampled)):
		c_size = cluster_sizes_sampled[c_index]
		output_xvg.write("@ s" + str(c) + " legend \"nb " + str(c_size) + "\"\n")
		output_txt.write("1_4_clusterprot_1D_smoothed.xvg," + str(len(cluster_sizes_sampled) + c_index + 1) + ",nb " + str(c_size) + "," + mcolors.rgb2hex(mcolorconv.to_rgb(colours_sizes_dict[c_size])) + "\n")
	output_txt.close()
	#write results
	for f_index in range(0, len(time_smoothed)):
		results = str(frames_time_smoothed[f_index])
		for c_size in cluster_sizes_sampled:
			results += "	" + str(round(cluster_sizes_pc_smoothed[c_size][f_index],2))
		for c_size in cluster_sizes_sampled:
			results += "	" + str(round(cluster_sizes_nb_smoothed[c_size][f_index],2))
		output_xvg.write(results + "\n")
	output_xvg.close()

	return
def graph_xvg_sizes_smoothed():											#DONE
	
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
	p_upper = {}
	for c_size in cluster_sizes_sampled:
		p_upper[c_size]=plt.plot(frames_time_smoothed, cluster_sizes_pc_smoothed[c_size], color = colours_sizes_dict[c_size], linewidth = 2.0, label = str(c_size))
	fontP.set_size("small")
	ax1.legend(prop=fontP)
	plt.title("%", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('% of proteins', fontsize="small")

	#plot data: nb
	#-------------
	ax2 = fig.add_subplot(212)
	p_lower = {}
	for c_size in cluster_sizes_sampled:
		p_lower[c_size] = plt.plot(frames_time_smoothed, cluster_sizes_nb_smoothed[c_size], color = colours_sizes_dict[c_size], linewidth = 2.0, label = str(c_size))
	fontP.set_size("small")
	ax2.legend(prop=fontP)
	plt.title("nb", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('nb of clusters', fontsize="small")

	#save figure
	#-----------
	ax1.set_ylim(0, 100)
	ax2.set_ylim(0, max(max(cluster_sizes_nb_smoothed.values()))+1)
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
def graph_aggregation_2D_sizes():										#TO CHECK

	#create filenames
	filename_png = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_1_plots_2D/png/1_1_clusterprot_2D.png'
	filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_1_plots_2D/1_1_clusterprot_2D.svg'

	#build color map
	color_map = mcolors.LinearSegmentedColormap.from_list('custom', colours_sizes_list, colours_sizes_nb)

	#set colours for surfacic peptide
	color_map.set_under(colour_leaflet_lower)
	color_map.set_over(colour_leaflet_upper)
	
	#determine nb of colours and their boundaries
	bounds = []
	cb_ticks_lab = []
	if args.cutoff_leaflet != "no":
		cb_ticks_lab.append("lower")
	for c in range(colours_sizes_range[0],colours_sizes_range[1]+1):
		bounds.append(c-0.5)
		cb_ticks_lab.append(str(c))
	if args.cutoff_leaflet != "no":
		cb_ticks_lab.append("upper")
	bounds.append(sorted(colours_sizes_dict.iterkeys())[-1]+0.5)	
	norm = mpl.colors.BoundaryNorm(bounds, color_map.N)
				
	#create figure ('norm' requires at least 2 elements to work)
	fig = plt.figure(figsize=(9, 8))
	ax_plot = fig.add_axes([0.10, 0.1, 0.75, 0.77])	
	ax_plot.matshow(proteins_cluster_status_sizes, origin = 'lower', interpolation = 'nearest', cmap = color_map, aspect = 'auto', norm = norm)

	#create color bar
	ax_cbar = fig.add_axes([0.88, 0.1, 0.025, 0.77])
	if args.cutoff_leaflet != "no":
		extend_val = "both"
		boundaries_val = np.concatenate([[bounds[0]-1], bounds, [bounds[-1]+1]])
	else:
		extend_val = "neither"
		boundaries_val = bounds
	cb = mpl.colorbar.ColorbarBase(ax_cbar, cmap = color_map, norm = norm, boundaries = boundaries_val, extend = extend_val)
	
	#position and label color bar ticks
	cb_ticks_pos = []
	if args.cutoff_leaflet != "no":
		cb_ticks_pos.append(bounds[0])
	for b in range(1,len(bounds)):
		cb_ticks_pos.append(bounds[b-1]+(bounds[b]-bounds[b-1])/2)
	if args.cutoff_leaflet != "no":
		cb_ticks_pos.append(bounds[-1])
	cb.set_ticks(cb_ticks_pos)
	cb.set_ticklabels(cb_ticks_lab)
	for t in cb.ax.get_yticklabels():
		t.set_fontsize('xx-small')
	
	#format axes
	ax_plot.xaxis.set_label_position('bottom') 
	ax_plot.xaxis.set_ticks_position('bottom')
	ax_plot.yaxis.set_ticks_position('left')
	ax_plot.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x,p: '{0:0g}'.format(x+1)))	#increase the index by 1 to get 1-based numbers
	ax_plot.set_xlabel("time (ns)", fontsize = "medium")
	ax_plot.set_ylabel("protein #", fontsize = "medium")
	ax_plot.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax_plot.yaxis.set_major_locator(MaxNLocator(prune = 'lower'))	
	ax_plot.set_title("Evolution of the cluster size in which proteins are involved", fontsize = "medium")
	plt.setp(ax_plot.xaxis.get_majorticklabels(), fontsize = "small" )
	plt.setp(ax_plot.yaxis.get_majorticklabels(), fontsize = "small" )
	ax_cbar.set_ylabel('cluster size', fontsize = 'small')
	
	#save figure
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
			
	return

#groups
def write_xvg_groups():													

	filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_2_plots_1D/xvg/2_2_clusterprot_1D.txt'
	filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_2_plots_1D/xvg/2_2_clusterprot_1D.xvg'
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
	output_xvg.write("@ legend length " + str(len(cluster_groups_sampled)*2) + "\n")
	#write caption: %
	for c_index in range(0,len(cluster_groups_sampled)):
		c_group = cluster_groups_sampled[c_index]
		output_xvg.write("@ s" + str(c_index) + " legend \"% " + str(c_group) + "\"\n")
		output_txt.write("1_2_clusterprot_1D.xvg," + str(c_index+1) + ",% " + str(c_group) + "," + mcolors.rgb2hex(mcolorconv.to_rgb(colours_groups_dict[c_group])) + "\n")
	#write caption: nb
	for c_index in range(0,len(cluster_groups_sampled)):
		c_group = cluster_groups_sampled[c_index]
		output_xvg.write("@ s" + str(c) + " legend \"nb " + str(c_group) + "\"\n")
		output_txt.write("1_2_clusterprot_1D.xvg," + str(len(cluster_groups_sampled) + c_index + 1) + ",nb " + str(c_group) + "," + mcolors.rgb2hex((colours_groups_dict[c_group])) + "\n")
	output_txt.close()
	#write results
	for f_index in range(0,nb_frames_to_process):
		results = str(frames_time[f_index])
		for c_group in cluster_groups_sampled:
			results+="	" + str(round(cluster_groups_pc[c_group][f_index],2))
		for c_group in cluster_groups_sampled:
			results+="	" + str(round(cluster_groups_nb[c_group][f_index],2))
		output_xvg.write(results + "\n")
	output_xvg.close()
	
	return
def graph_xvg_groups():													
	
	#create filenames
	#----------------
	filename_png = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_2_plots_1D/png/2_2_clusterprot_1D.png'
	filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_2_plots_1D/2_2_clusterprot_1D.svg'

	#create figure
	#-------------
	fig = plt.figure(figgroup=(8, 6.2))
	fig.suptitle("Evolution of proteins distribution")
		
	#plot data: %
	#------------
	ax1 = fig.add_subplot(211)
	p_upper = {}
	for c_group in cluster_groups_sampled:
		p_upper[c_group] = plt.plot(frames_time, cluster_groups_pc[c_group], color = colours_groups_dict[c_group], linewidth = 2.0, label = str(c_group))
	fontP.set_group("small")
	ax1.legend(prop=fontP)
	plt.title("%", fontgroup="small")
	plt.xlabel('time (ns)', fontgroup="small")
	plt.ylabel('% of proteins', fontgroup="small")

	#plot data: nb
	#-------------
	ax2 = fig.add_subplot(212)
	p_lower = {}
	for c_group in cluster_groups_sampled:
		p_lower[c_group] = plt.plot(frames_time, cluster_groups_nb[c_group], color = colours_groups_dict[c_group], linewidth = 2.0, label=str(c_group))
	fontP.set_group("small")
	ax2.legend(prop=fontP)
	plt.title("nb", fontgroup="small")
	plt.xlabel('time (ns)', fontgroup="small")
	plt.ylabel('nb of clusters', fontgroup="small")

	#save figure
	#-----------
	ax1.set_ylim(0, 100)
	ax2.set_ylim(0, max(max(cluster_groups_nb.values()))+1)
	ax1.spines['top'].set_visible(False)
	ax1.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.spines['right'].set_visible(False)
	ax1.xaxis.set_ticks_position('bottom')
	ax1.yaxis.set_ticks_position('left')
	ax2.xaxis.set_ticks_position('bottom')
	ax2.yaxis.set_ticks_position('left')
	ax1.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax1.yaxis.set_major_locator(MaxNLocator(nbins=5))
	ax2.xaxis.set_major_locator(MaxNLocator(nbins=5))
	ax2.yaxis.set_major_locator(MaxNLocator(nbins=5))
	plt.setp(ax1.xaxis.get_majorticklabels(), fontgroup="small" )
	plt.setp(ax1.yaxis.get_majorticklabels(), fontgroup="small" )
	plt.setp(ax2.xaxis.get_majorticklabels(), fontgroup="small" )
	plt.setp(ax2.yaxis.get_majorticklabels(), fontgroup="small" )	
	plt.subplots_adjust(top=0.9, bottom=0.07, hspace=0.37, left=0.09, right=0.96)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()

	return
def write_xvg_groups_smoothed():											
	filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_4_plots_1D_smoothed/xvg/2_4_clusterprot_1D_smoothed.txt'
	output_txt = open(filename_txt, 'w')
	output_txt.write("@[lipid tail order parameters statistics - written by order_param v" + str(version_nb) + "]\n")
	output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 2_4_clusterprot_1D_smoothed.xvg.\n")
	filename_xvg=os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_4_plots_1D_smoothed/xvg/2_4_clusterprot_1D_smoothed.xvg'
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
	output_xvg.write("@ legend length " + str(len(cluster_groups_sampled)*2) + "\n")
	#write caption: %
	for c_index in range(0,len(cluster_groups_sampled)):
		c_group = cluster_groups_sampled[c_index]
		output_xvg.write("@ s" + str(c_index) + " legend \"% " + str(c_group) + "\"\n")
		output_txt.write("2_4_clusterprot_1D_smoothed.xvg," + str(c_index+1) + ",% " + str(c_group) + "," + mcolors.rgb2hex((colours_groups_dict[c_group])) + "\n")
	#write caption: nb
	for c_index in range(0,len(cluster_groups_sampled)):
		c_group = cluster_groups_sampled[c_index]
		output_xvg.write("@ s" + str(c) + " legend \"nb " + str(c_group) + "\"\n")
		output_txt.write("2_4_clusterprot_1D_smoothed.xvg," + str(len(cluster_groups_sampled) + c_index + 1) + ",nb " + str(c_group) + "," + mcolors.rgb2hex(mcolorconv.to_rgb(colours_groups_dict[c_group])) + "\n")
	output_txt.close()
	#write results
	for f_index in range(0, len(time_smoothed)):
		results = str(frames_time_smoothed[f_index])
		for c_group in cluster_groups_sampled:
			results += "	" + str(round(cluster_groups_pc_smoothed[c_group][f_index],2))
		for c_group in cluster_groups_sampled:
			results += "	" + str(round(cluster_groups_nb_smoothed[c_group][f_index],2))
		output_xvg.write(results + "\n")
	output_xvg.close()

	return
def graph_xvg_groups_smoothed():											
	
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
	p_upper = {}
	for c_size in cluster_sizes_sampled:
		p_upper[c_size]=plt.plot(frames_time_smoothed, cluster_sizes_pc_smoothed[c_size], color = colours_sizes_dict[c_size], linewidth = 2.0, label = str(c_size))
	fontP.set_size("small")
	ax1.legend(prop=fontP)
	plt.title("%", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('% of proteins', fontsize="small")

	#plot data: nb
	#-------------
	ax2 = fig.add_subplot(212)
	p_lower = {}
	for c_size in cluster_sizes_sampled:
		p_lower[c_size] = plt.plot(frames_time_smoothed, cluster_sizes_nb_smoothed[c_size], color = colours_sizes_dict[c_size], linewidth = 2.0, label = str(c_size))
	fontP.set_size("small")
	ax2.legend(prop=fontP)
	plt.title("nb", fontsize="small")
	plt.xlabel('time (ns)', fontsize="small")
	plt.ylabel('nb of clusters', fontsize="small")

	#save figure
	#-----------
	ax1.set_ylim(0, 100)
	ax2.set_ylim(0, max(max(cluster_sizes_nb_smoothed.values()))+1)
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


def graph_aggregation_2D_groups():
	
	#create filenames
	filename_png=os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_1_plots_2D//png/2_1_clusterprot_2D.png'
	filename_svg=os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_1_plots_2D//2_1_clusterprot_2D.svg'

	#create data
	lip_2D_evolution=np.zeros((proteins_nb,len(time_stamp.keys())))
	for p_index in range(0,proteins_nb):
		lip_2D_evolution[p_index,:]=np.asarray(proteins_cluster_group[p_index])

	#build color map
	color_map=mcolors.LinearSegmentedColormap.from_list('custom', colours_groups_list, colours_groups_nb)
	
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
		tmp_xticks_lab.append('{0:0g}'.format(np.floor(t)))
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
def write_stability_groups():
	
	filename_details = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_0_clusterprot_stability.stat'
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
	
	#group info
	tmp_cap1 = ""
	tmp_cap2 = "-----"
	for g_index in cluster_TM_groups_sampled:
		if g_index == groups_number:
			tmp_cap1 += "	other"
		elif groups_boundaries[g_index][1]==100000:
			tmp_cap1 += "	>=" + str(groups_boundaries[g_index][0])
		else:
			tmp_cap1+="	" + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1])
		tmp_cap2 += "--------"

	#results
	output_stat.write("\n")
	output_stat.write(tmp_cap1 + "\n")
	output_stat.write(tmp_cap2 + "\n")
	results = str("")
	for g_index in cluster_groups_sampled:
		results += "	" + str(proteins_groups_stability[g_index])
	output_stat.write(results + "\n")

	output_stat.close()
	
	return

#annotations
def write_frame_stat(f_nb, f_index, f_t):								#DONE

	#case: gro file or xtc summary
	#-----------------------------
	if f_index == "all":
		#sizes
		#-----
		#create file
		if args.xtcfilename == "no":
			filename_details = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/' + args.grofilename[:-4] + '_annotated_clustprot_sizes_sampled.stat'		
		else:
			filename_details = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/' + args.xtcfilename[:-4] + '_annotated_clustprot_sizes_sampled.stat'		
		output_stat = open(filename_details, 'w')		
	
		#general info
		output_stat.write("[protein clustering statistics - written by cluster_prot v" + str(version_nb) + "]\n")
		output_stat.write("\n")
		output_stat.write("1. nb of proteins: " + str(proteins_nb) + "\n")
		output_stat.write("\n")
		output_stat.write("2. cluster detection Method:\n")
		if args.m_algorithm == "min":
			output_stat.write(" - connectivity based (min distances)\n")
			output_stat.write(" - contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
		elif args.m_algorithm == "cog":
			output_stat.write(" - connectivity based (cog distances)\n")
			output_stat.write(" - contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
		else:
			output_stat.write(" - density based (DBSCAN)\n")
			output_stat.write(" - search radius = " + str(args.dbscan_dist) + " Angstrom, nb of neighbours = " + str(args.dbscan_nb) + "\n")
		if args.xtcfilename != "no":
			output_stat.write("\n")
			output_stat.write("3. frame:	" + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")
		
		#data results
		output_stat.write("\n")
		output_stat.write("Sizes sampled: ")
		tmp_res = ""
		for c_size in cluster_sizes_sampled:
			tmp_res += "," + str(c_size)
		output_stat.write(tmp_res[1:] + "\n")
		output_stat.write("\n")
		output_stat.close()
		
		#groups
		#------
		if args.cluster_groups_file != "no":
			#create file
			if args.xtcfilename == "no":
				filename_details = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/' + args.grofilename[:-4] + '_annotated_clustprot_groups_sampled.stat'		
			else:
				filename_details = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/' + args.xtcfilename[:-4] + '_annotated_clustprot_groups_sampled.stat'		
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
				output_stat.write("3. nb frames processed:	" + str(nb_frames_to_process) + " (" + str(nb_frames_xtc) + " frames in xtc, step=" + str(args.frames_dt) + ")\n")

			#group index definition
			output_stat.write("\n")
			output_stat.write("Group size ranges:\n")
			for g_index in range(0,groups_number+3):
				output_stat.write(str(g_index) + "=" + str(groups_labels[g_index]) + "\n")

			#data results
			output_stat.write("\n")
			output_stat.write("Groups sampled: ")
			tmp_res = ""
			for g_index in cluster_groups_sampled:
				tmp_res += "," + str(g_index)
			output_stat.write(tmp_res[1:] + "\n")		
			output_stat.write("\n")
		output_stat.close()
		
	#case: xtc snapshot
	#==================
	else:
		#sizes
		#-----
		#create file
		if args.xtcfilename == "no":
			filename_details=os.getcwd() + '/' + str(args.output_folder) + '/3_snapshots/sizes/' + args.grofilename[:-4] + '_annotated_clusterprot_sizes.stat'
		else:
			filename_details=os.getcwd() + "/" + str(args.output_folder) + '/3_snapshots/sizes/' + args.xtcfilename[:-4] + '_annotated_clusterprot_sizes_' + str(int(f_t)).zfill(5) + 'ns.stat'
		output_stat = open(filename_details, 'w')		
	
		#general info
		output_stat.write("[protein clustering statistics - written by cluster_prot v" + str(version_nb) + "]\n")
		output_stat.write("\n")
		output_stat.write("1. nb of proteins: " + str(proteins_nb) + "\n")
		output_stat.write("2. cluster detection Method:\n")
		if args.m_algorithm == "min":
			output_stat.write(" - connectivity based (min distances)\n")
			output_stat.write(" - contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
		elif args.m_algorithm == "cog":
			output_stat.write(" - connectivity based (cog distances)\n")
			output_stat.write(" - contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
		else:
			output_stat.write(" - density based (DBSCAN)\n")
			output_stat.write(" - search radius = " + str(args.dbscan_dist) + " Angstrom, nb of neighbours = " + str(args.dbscan_nb) + "\n")
		if args.xtcfilename != "no":
			output_stat.write("\n")
			output_stat.write("3. time: " + str(f_t) + "ns (frame " + str(f_nb) + "/" + str(nb_frames_xtc) + ")\n")

		#what's in this file
		output_stat.write("\n")
		output_stat.write("Distribution of proteins by cluster size:\n")
		output_stat.write("\n")
		tmp_cap1 = "cluster_size"
		tmp_cap2 = "-----------"
		tmp_res = "% proteins"
		for c_size in cluster_sizes_sampled:
			tmp_cap1 += "	" + str(c_size)
			tmp_cap2 += "--------"
		output_stat.write(tmp_cap1 + "\n")
		output_stat.write(tmp_cap2 + "\n")
		for c_size in cluster_sizes_sampled:		
			tmp_res += "	" + str(round(cluster_sizes_pc[c_size][f_index], 1))
		output_stat.write(tmp_res + "\n")		
		output_stat.close()
		
		#groups
		#======
		if args.cluster_groups_file != "no":
			#create file
			if args.xtcfilename == "no":
				filename_details = os.getcwd() + '/' + str(args.output_folder) + '/3_snapshots/groups/' + args.grofilename[:-4] + '_annotated_clusterprot_groups.stat'
			else:
				filename_details = os.getcwd() + "/" + str(args.output_folder) + '/3_snapshots/groups/' + args.xtcfilename[:-4] + '_annotated_clusterprot_groups_' + str(int(f_t)).zfill(5) + 'ns.stat'
			output_stat = open(filename_details, 'w')		
		
			#general info
			output_stat.write("[protein clustering statistics - written by cluster_prot v" + str(version_nb) + "]\n")
			output_stat.write("\n")
			output_stat.write("1. nb of proteins: " + str(proteins_nb) + "\n")
			output_stat.write("2. cluster detection Method:\n")
			if args.m_algorithm == "min":
				output_stat.write(" - connectivity based (min distances)\n")
				output_stat.write(" - contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
			elif args.m_algorithm == "cog":
				output_stat.write(" - connectivity based (cog distances)\n")
				output_stat.write(" - contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
			else:
				output_stat.write(" - density based (DBSCAN)\n")
				output_stat.write(" - search radius = " + str(args.dbscan_dist) + " Angstrom, nb of neighbours = " + str(args.dbscan_nb) + "\n")
			if args.xtcfilename != "no":
				output_stat.write("\n")
				output_stat.write("3. time: " + str(f_t) + "ns (frame " + str(f_nb) + "/" + str(nb_frames_xtc) + ")\n")
		
			#group index definition
			output_stat.write("\n")
			output_stat.write("Group size ranges:\n")
			for g_index in range(0,groups_number+3):
				output_stat.write(str(g_index) + "=" + str(groups_labels[g_index]) + "+\n")

			#what's in this file
			output_stat.write("\n")
			output_stat.write("Distribution of proteins by cluster size group:\n")
			output_stat.write("\n")
			tmp_cap1 = "cluster_group"
			tmp_cap2 = "-----------"
			tmp_res = "% proteins"
			for c_group in cluster_groups_sampled:
				tmp_cap1 += "	" + str(c_group)
				tmp_cap2 += "--------"
			output_stat.write(tmp_cap1 + "\n")
			output_stat.write(tmp_cap2 + "\n")
			for c_group in cluster_groups_sampled:		
				tmp_res += "	" + str(round(groups_pc[c_group][f_index],1))
			output_stat.write(tmp_res + "\n")		
			output_stat.close()

	return
def write_frame_snapshot(f_index, f_t):									#DONE
	
	#sizes
	#-----
	#store cluster info in beta factor field
	for p_index in range(0,proteins_nb):
			proteins_sele[p_index].set_bfactor(proteins_cluster_status_sizes[p_index,f_index])
	
	#write annotated file
	if args.xtcfilename == "no":
		all_atoms.write(os.getcwd() + '/' + str(args.output_folder) + '/3_snapshots/sizes/' + args.grofilename[:-4] + '_annotated_clusterprot_sizes', format="PDB")
	else:
		tmp_name = os.getcwd() + "/" + str(args.output_folder) + '/3_snapshots/sizes/' + args.xtcfilename[:-4] + '_annotated_clusterprot_sizes_' + str(int(f_t)).zfill(5) + 'ns.pdb'
		W = Writer(tmp_name, nb_atoms)
		W.write(all_atoms)
	
	#groups
	#------
	if args.cluster_groups_file != "no":
		#store cluster info in beta factor field
		for p_index in range(0,proteins_nb):
				proteins_sele[p_index].set_bfactor(proteins_cluster_status_groups[p_index,f_index])
		
		#write annotated file
		if args.xtcfilename == "no":
			all_atoms.write(os.getcwd() + '/' + str(args.output_folder) + '/3_snapshots/groups/' + args.grofilename[:-4] + '_annotated_clusterprot_groups', format = "PDB")
		else:
			tmp_name=os.getcwd() + "/" + str(args.output_folder) + '/3_snapshots/groups/' + args.xtcfilename[:-4] + '_annotated_clusterprot_groups_' + str(int(f_t)).zfill(5) + 'ns.pdb'
			W = Writer(tmp_name, nb_atoms)
			W.write(all_atoms)
		
	return
def write_frame_annotation(f_index, f_t):								#DONE

	#sizes
	#-----
	#create file
	if args.xtcfilename == "no":
		filename_details=os.getcwd() + '/' + str(args.output_folder) + '/3_snapshots/sizes/' + args.grofilename[:-4] + '_annotated_clusterprot_sizes.txt'
	else:
		filename_details=os.getcwd() + "/" + str(args.output_folder) + '/3_snapshots/sizes/' + args.xtcfilename[:-4] + '_annotated_clusterprot_sizes_' + str(int(f_t)).zfill(5) + 'ns.txt'
	output_stat = open(filename_details, 'w')		

	#output VMD protein selection line
	tmp_prot_sele = ""
	for p_index in range(0,proteins_nb):
		tmp_prot_sele += "." + proteins_sele_string_VMD[p_index]
	output_stat.write(tmp_prot_sele[1:] + "\n")
				
	#ouput min and max size
	output_stat.write(str(np.min(proteins_cluster_status_sizes[:,f_index])) + "." + str(np.max(proteins_cluster_status_sizes[:,f_index])) + "\n")
	
	#ouptut cluster size for each protein
	tmp_sizes = "1"
	for p_index in range(0,proteins_nb):
		tmp_sizes += "." + str(proteins_cluster_status_sizes[p_index,f_index])
	output_stat.write(tmp_sizes + "\n")
	output_stat.close()
	
	#groups
	#------
	if args.cluster_groups_file != "no":
		#create file
		if args.xtcfilename == "no":
			filename_details = os.getcwd() + '/' + str(args.output_folder) + '/3_snapshots/groups/' + args.grofilename[:-4] + '_annotated_clusterprot_groups.txt'
		else:
			filename_details = os.getcwd() + "/" + str(args.output_folder) + '/3_snapshots/groups/' + args.xtcfilename[:-4] + '_annotated_clusterprot_groups_' + str(int(f_t)).zfill(5) + 'ns.txt'
		output_stat = open(filename_details, 'w')		
	
		#output VMD protein selection line
		tmp_prot_sele = ""
		for p_index in range(0,proteins_nb):
			tmp_prot_sele += "." + proteins_sele_string_VMD[p_index]
		output_stat.write(tmp_prot_sele[1:] + "\n")
		
		#ouput min and max size
		output_stat.write(str(np.min(proteins_cluster_status_groups[:,f_index])) + "." + str(np.max(proteins_cluster_status_groups[:,f_index])) + "\n")
		
		#ouptut cluster size for each protein
		tmp_groups = "1"
		for p_index in range(0,proteins_nb):
			tmp_groups += "." + str(proteins_cluster_status_groups[p_index,f_index])
		output_stat.write(tmp_groups + "\n")
		output_stat.close()
	
	return
def write_xtc_annotation(action):										#DONE
	
	global output_xtc_annotate_cluster_size
	
	#create file name
	#----------------
	if action == "initialise":
		output_xtc_annotate_cluster_size = os.getcwd() + '/' + str(args.output_folder) + '/4_VMD/' + args.xtcfilename[:-4] + '_annotated_clustprot_sizes.txt'
		if args.cluster_groups_file != "no":
			global output_xtc_annotate_cluster_group
			output_xtc_annotate_cluster_group = os.getcwd() + '/' + str(args.output_folder) + '/4_VMD/' + args.xtcfilename[:-4] + '_annotated_clustprot_groups.txt'

	#pre-prend info to file
	#----------------------
	elif action == "finish":
	
		#add remaining bit of info and read the whole file content before closing the file
		with open(output_xtc_annotate_cluster_size, 'a') as f:
			f.write(vmd_cluster_size)		
		#read content of file
		with open(output_xtc_annotate_cluster_size, 'r') as f:
			tmp_data = f.read()
		#pre-pend info to it
		with open(output_xtc_annotate_cluster_size, 'w') as f:
			#output VMD protein selection line
			tmp_prot_sele = ""
			for p_index in range(0,proteins_nb):
				tmp_prot_sele += "." + proteins_sele_string_VMD[p_index]
			f.write(tmp_prot_sele[1:] + "\n")
			#ouput min and max size
			f.write(str(min(cluster_sizes_sampled)) + "." + str(max(cluster_sizes_sampled)) + "\n")
			#write data
			f.write(tmp_data)
	
		if args.cluster_groups_file != "no":
			#add remaining bit of info and read the whole file content before closing the file
			with open(output_xtc_annotate_cluster_group, 'a') as f:
				f.write(vmd_cluster_group)
			#read content of file
			with open(output_xtc_annotate_cluster_group, 'r') as f:
				tmp_data = f.read()
			#pre-pend info to it
			with open(output_xtc_annotate_cluster_group, 'w') as f:
				#output VMD protein selection line
				for p_index in range(0,proteins_nb):
					tmp_prot_sele += "." + proteins_sele_string_VMD[p_index]
				output_stat.write(tmp_prot_sele[1:] + "\n")
				#ouput min and max size
				f.write(str(min(cluster_groups_sampled)) + "." + str(max(cluster_groups_sampled)) + "\n")
				#write data
				f.write(tmp_data)

	return

##########################################################################################
# ALGORITHM
##########################################################################################

#=========================================================================================
#process inputs
#=========================================================================================

set_lipids_beads()
load_MDA_universe()
if args.selection_file_ff != "no":
	identify_ff()
if args.cutoff_leaflet != "no":
	identify_leaflets()
identify_proteins()
if args.cluster_groups_file != "no":
	initialise_groups()

#=========================================================================================
# initialise data structures
#=========================================================================================
print "\nInitialising data structures..."
data_struct_time()
data_struct_stats()
data_struct_proteins()
if args.xtcfilename != "no" and args.buffer_size != -1:
	write_xtc_annotation("initialise")

#=========================================================================================
# generate data
#=========================================================================================
print "\nDetecting proteins clusters..."

#case: gro file
#==============
if args.xtcfilename == "no":
	frames_nb[0] = 1
	frames_time[0] = 0

	#detect clusters
	#---------------
	if args.m_algorithm != "density":
		clusters = detect_clusters_connectivity(get_distances(U.trajectory.ts.dimensions), U.dimensions)
	else:
		clusters = detect_clusters_density(get_distances(U.trajectory.ts.dimensions), U.dimensions)
	
	#assign current cluster status
	#-----------------------------
	if args.cutoff_leaflet == "no":
		process_clusters(clusters, 0)
	else:
		process_clusters_TM(clusters, 0)
			
#case: xtc file
#==============
else:
	for f_index in range(0,nb_frames_to_process):
		ts = U.trajectory[frames_to_process[f_index]]
		if ts.time/float(1000) > args.t_end:
			break
		progress = '\r -processing frame ' + str(ts.frame) + '/' + str(nb_frames_xtc) + '                      '  
		sys.stdout.flush()
		sys.stdout.write(progress)

		#frame properties
		f_time = ts.time/float(1000)
		f_nb = ts.frame
		frames_nb[f_index] = f_nb
		frames_time[f_index] = f_time
		f_write = frames_to_write[f_index]
		box_dim = U.trajectory.ts.dimensions

		#detect clusters
		#---------------
		if args.m_algorithm != "density":
			clusters = detect_clusters_connectivity(get_distances(U.trajectory.ts.dimensions), box_dim)
		else:
			clusters = detect_clusters_density(get_distances(U.trajectory.ts.dimensions), box_dim)
		
		#assign current cluster status
		#-----------------------------
		if args.cutoff_leaflet == "no":
			process_clusters(clusters, f_index)
		else:
			process_clusters_TM(clusters, f_index, box_dim)
			
		#output results
		#--------------
		if f_write:
			print "  (writing snapshot...)"
			write_frame_stat(f_nb, f_index, f_time)
			write_frame_snapshot(f_index, f_time)
			write_frame_annotation(f_index, f_time)

		#buffer counter for outputting xtc annotation files
		#--------------------------------------------------
		if args.buffer_size != -1:
			if vmd_counter == args.buffer_size:
				vmd_counter = 0
			else:
				vmd_counter += 1 


	print ""

#=========================================================================================
# process data
#=========================================================================================
print "\nCalculating statistics..."
get_sizes_sampled()
update_color_dict()
calculate_statistics()
if args.nb_smoothing > 1:
	smooth_data()

#=========================================================================================
# produce outputs
#=========================================================================================
print "\nWriting outputs..."

#case: gro file
#==============
if args.xtcfilename == "no":
	if len(cluster_sizes_sampled)>1:
		print " -writing statistics..."
		write_frame_stat("all", 0, 0)
		print " -writing annotated pdb..."
		write_frame_snapshot(0, 0)
		write_frame_annotation(0, 0)
	else:
		print "\n"
		print "Warning: a single cluster size (", str(cluster_sizes_sampled[0]), ") was detected throughout the trajectory. Check the --algorithm, --cutoff, --radius or --neighbours options (see cluster_prot -h)."
		write_warning()

#case: xtc file
#==============
else:
	if len(cluster_sizes_sampled)>1:
		#writing statistics
		print " -writing statistics..."
		write_frame_stat("all", 0, 0)
				
		#write annotation files for VMD
		if args.buffer_size != -1:
			write_xtc_annotation("finish")
		
		#write xvg and graphs
		print " -writing xvg and graphs..."
		graph_aggregation_2D_sizes()
		write_xvg_biggest()
		graph_xvg_biggest()
		write_xvg_sizes_TM()
		graph_xvg_sizes_TM()
		write_xvg_sizes_interfacial()
		graph_xvg_sizes_interfacial()
		if args.nb_smoothing > 1:
			write_xvg_sizes_smoothed()
			graph_xvg_sizes_smoothed()
			write_xvg_cluster_biggest_smoothed()
			graph_xvg_cluster_biggest_smoothed()
		if args.cluster_groups_file != "no":
			graph_aggregation_2D_groups()
			write_stability_groups()
			write_xvg_groups()
			graph_xvg_groups()
			if args.nb_smoothing > 1:
				write_xvg_groups_smoothed()
				graph_xvg_groups_smoothed()
	else:
		print "\n"
		print "Warning: a single cluster size (", str(cluster_sizes_sampled[0]), ") was detected throughout the trajectory, maybe check the cluster detection options (see cluster_prot -h)."
		write_warning()
		
#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check output in ./" + args.output_folder + "/"
print ""
sys.exit(0)
