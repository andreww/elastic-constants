#!/usr/bin/env python
"""
generate_strain_gulp.py

Generate strained gulp input files and 
.cijdat files for elastic constants calculation
for analysis with elastics.py

Copyright (c) 2011 Andrew Walker (a.walker@ucl.ac.uk)
All rights reserved.
"""

import sys
import numpy as np

import gulp
import generate_strain as gs

version = 0.1
			
def main(input_options, libmode=False):

	# deal with options
	options, arguments = gs.get_options(input_options, libmode)
	seedname = arguments[0]
	
	(cell,pointgroup,atoms) = gulp.parse_dotgot(seedname)

	# Re-align lattice vectors on cartisian system
	a, b, c, al, be, ga = gs.cellCART2cellABC(cell)
	cell = gs.cellABC2cellCART(a, b, c, al, be, ga)


	# Not sure why the lattice types are enumerated like this, but this is how .cijdat does it...
	latticeTypes = {0:"Unknown", 1:"Triclinic", 2:"Monoclinic", 3:"Orthorhombic", \
	                4:"Tetragonal", 5:"Cubic", 6:"Trigonal-low", 7:"Trigonal-high/Hexagonal"}
	maxstrain = options.strain
	if (maxstrain == None):
		maxstrain = 0.1
	numsteps = options.numsteps
	if (numsteps == None):
		numsteps = 3 
	# Which strain pattern to use?
	if (options.lattice == None):
		if (pointgroup == None):
			# Nothing from user and nothing from 
			# .got: we are in trouble
			print "No point group found in .got file so the strain pattern cannot be determined\n"
			print "A strain pattern can also be provided using the -l flag\n"
			sys.exit(1)
		else:
			# Use the value from .got
			latticeCode = gs.PointGroup2StrainPat(pointgroup)
	else:
		if (pointgroup == None):
			# Noting in .got - use users choice
			latticeCode = options.lattice
		else:
			# Use users choice, but check and warn
			latticeCode = options.lattice
			if (latticeCode != PointGroup2StrainPat(pointgroup)):
				print "WARNING: User supplied lattice code is inconsistant with the point group\n"
				print "         used in the bulk GULP run. Using user supplied lattice code.\n"
		
	print "Cell parameters: a = %f gamma = %f" % (a, al)
	print "                 b = %f beta  = %f" % (b, be)
	print "                 c = %f gamma = %f \n" % (c, ga)
	print "Lattce vectors:  %7f %7f %7f " % (cell[0][0], cell[0][1], cell[0][2])
	print "                 %7f %7f %7f " % (cell[1][0], cell[1][1], cell[1][2])
	print "                 %7f %7f %7f \n " % (cell[2][0], cell[2][1], cell[2][2])
	patterns = gs.GetStrainPatterns(latticeCode)
	numStrainPatterns = len(patterns)
	print "Lattice type is ", latticeTypes[latticeCode] 
	print "Number of patterns: "+ str(numStrainPatterns) +"\n"

	cijdat = open(seedname+".cijdat","w")
	print "Writing strain data to ", seedname+".cijdat"
	cijdat.write(str(latticeCode) + ' ' + str(numsteps*2) + ' 0 ' + '0 \n')
  	cijdat.write(str(maxstrain)+"\n")

	# The order of these three loops matters for the analysis code.
	for patt in range(numStrainPatterns):
		this_pat = patterns[patt]
		for a in range(0,numsteps):
			for neg in range(0,2):
				if (neg == 0):
					this_mag = (float(a)+1) / (float(numsteps)) * maxstrain
				else:
					this_mag = -1.0 * (float(a)+1) / (float(numsteps)) * maxstrain
				disps = [x * this_mag for x in patterns[patt]]
				# Build the strain tensor (IRE convention but 1 -> 0 etc.)
				this_strain = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
				# diagonal elements - strain is displacment / lattice vector length
				this_strain[0] = disps[0] / np.sqrt(cell[0][0]**2+cell[0][1]**2+cell[0][2]**2)
				this_strain[1] = disps[1] / np.sqrt(cell[1][0]**2+cell[1][1]**2+cell[1][2]**2)
				this_strain[2] = disps[2] / np.sqrt(cell[2][0]**2+cell[2][1]**2+cell[2][2]**2)
				# off diagonals - we only strain upper right corner of cell matrix, so strain is 1/2*du/dx...
				this_strain[3] = 0.5 * (disps[3] / np.sqrt(cell[1][0]**2+cell[1][1]**2+cell[1][2]**2))
				this_strain[4] = 0.5 * (disps[4] / np.sqrt(cell[0][0]**2+cell[0][1]**2+cell[0][2]**2))
				this_strain[5] = 0.5 * (disps[5] / np.sqrt(cell[0][0]**2+cell[0][1]**2+cell[0][2]**2))

				# Deform cell - only apply deformation to upper right corner
				defcell = [[cell[0][0]+disps[0], cell[0][1]+disps[5], cell[0][2]+disps[4]],
					   [cell[1][0],          cell[1][1]+disps[1], cell[1][2]+disps[3]],
					   [cell[2][0],          cell[2][1],          cell[2][2]+disps[2]]]

				pattern_name = seedname + "_cij__" + str(patt+1) + "__" + str((a*2)+1+neg)

				print "Pattern Name = ", pattern_name
				print "Pattern = ", this_pat
				print "Magnitude = ", this_mag
				cijdat.write(pattern_name+"\n")
				cijdat.write(str(this_strain[0]) + " " + str(this_strain[5]) + " " + str(this_strain[4]) + "\n")
				cijdat.write(str(this_strain[5]) + " " + str(this_strain[1]) + " " + str(this_strain[3]) + "\n")
				cijdat.write(str(this_strain[4]) + " " + str(this_strain[3]) + " " + str(this_strain[2]) + "\n")
				gulp.produce_dotgin(seedname, pattern_name+".gin", defcell, atoms)
	

	

if __name__ == "__main__":
	main(sys.argv[1:])
