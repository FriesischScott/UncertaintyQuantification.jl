# Reference: https://opensees.berkeley.edu/wiki/index.php?title=Time_History_Analysis_of_a_2D_Elastic_Cantilever_Column
# --------------------------------------------------------------------------------------------------
# Example 1. cantilever 2D
# EQ ground motion with gravity
# all units are in kip, inch, second
# elasticBeamColumn ELEMENT
#		Silvia Mazzoni & Frank McKenna, 2006
# 		Marius Bittner, 2024
#
#    ^Y
#	 |
#	 3	
#    |          
#    |          
#    |
#    2       __ 
#    |         | 
#    |         | 
#    |         | 
#  (1)      36'
#    |         | 
#    |         | 
#    |         | 
#  =1=    ----  -------->X
#

# SET UP ----------------------------------------------------------------------------
wipe;						       # clear opensees model
model basic -ndm 2 -ndf 3;	       # 2 dimensions, 3 dof per node

# define GEOMETRY -------------------------------------------------------------
# nodal coordinates:
node 1 0. 0.;					   # node#, X Y
node 2 0. 216.;
node 3 0. 432.;

# Single point constraints -- Boundary Conditions
fix 1 1 1 1; 			           # node DX DY RZ

# nodal masses:
mass 2 2.59 0. 0.;			   # node#, Mx My Mz, Mass=Weight/g.
mass 3 5.18 0. 0.;

# Define ELEMENTS -------------------------------------------------------------
# define geometric transformation: performs a linear geometric transformation of beam stiffness and resisting force from the basic system to the global-coordinate system
geomTransf Linear 1;  		       # associate a tag to transformation

# connectivity:
element elasticBeamColumn 1 1 2 3600 3225 1080000 1;	# element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
element elasticBeamColumn 2 2 3 3600 3225 1080000 1;

# Define RECORDERS -------------------------------------------------------------
recorder Node -file displacement.out -time -node 3 -dof 1 disp;			            # displacements of free nodes

# define GRAVITY -------------------------------------------------------------
timeSeries Linear 1
pattern Plain 1 1 {
   load 3 0. -2000. 0.;			    # node#, FX FY MZ --  superstructure-weight
}
constraints Plain;     				# how it handles boundary conditions
numberer Plain;					    # renumber dof's to minimize band-width (optimization), if you want to
system BandGeneral;				    # how to store and solve the system of equations in the analysis
algorithm Linear;                   # use Linear algorithm for linear analysis
integrator LoadControl 0.1;			# determine the next time step for an analysis, # apply gravity in 10 steps
analysis Static					    # define type of analysis static or transient
analyze 10;					        # perform gravity analysis
loadConst -time 0.0;				# hold gravity constant and restart time

# DYNAMIC ground-motion analysis -------------------------------------------------------------
# create load pattern
set G 386
# Time step of signal injected as parameter
timeSeries Path 2 -dt {{{ :dt }}} -filePath ground-motion.dat -factor $G; # define acceleration vector from file (dt is associated with the input file gm)
pattern UniformExcitation 2 1 -accel 2;		         # define where and how (pattern tag, dof) acceleration is applied

# set damping based on first eigen mode
set freq [expr [eigen -fullGenLapack 1]**0.5]
set dampRatio 0.05
rayleigh 0. 0. 0. [expr 2*$dampRatio/$freq]

# create the analysis
wipeAnalysis;					     # clear previously-define analysis parameters
constraints Plain;     				 # how it handles boundary conditions
numberer Plain;					     # renumber dof's to minimize band-width (optimization), if you want to
system BandGeneral;					 # how to store and solve the system of equations in the analysis
algorithm Linear					 # use Linear algorithm for linear analysis
integrator Newmark 0.5 0.25 ;	     # determine the next time step for an analysis
analysis Transient;					 # define type of analysis: time-dependent
analyze 3995 0.01;					 # apply 3995 0.01-sec time steps in analysis


puts "Done!"
wipe