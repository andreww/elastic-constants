#!/usr/bin/env python
"""
generate_strain.py

Generate strained castep input files and 
.cijdat files for elastic constants calculation
for analysis with elastics.py

Copyright (c) 2010 Andrew Walker (a.walker@ucl.ac.uk)
All rights reserved.
"""

import sys
import os
import optparse
import re
import numpy as np
import castep

version = 0.1

def PointGroup2StrainPat(pointgroup):
	"""
	Converts point group number (as ordered by CASTEP
	and in the International Tables) to a number representing 
	the needed strain pattern.
	"""
        supcode = 0
	if (pointgroup < 1):
		print "Point group number " + str(pointgroup) + " not recognized.\n"
		sys.exit(1)
	elif (pointgroup <= 2): 
		# Triclinic
		patt = 1
	elif (pointgroup <= 5):
		# Monoclinic
		patt = 2
	elif (pointgroup <= 8):
		# Orthorhombic
		patt = 3
	elif (pointgroup <= 15):
		# Tetragonal
		patt = 4
	elif (pointgroup <= 17):
		# Trigonal-Low
		patt = 6
	elif (pointgroup <= 20):
		# Trigonal-High
		patt = 7
                supcode = 1
	elif (pointgroup <= 27):
		# Hexagonal
		patt = 7
	elif (pointgroup <= 32):
		# Cubic
		patt = 5
	else:
		print "Point group number " + str(pointgroup) + " not recognized.\n"
		sys.exit(1)
	return patt, supcode

def GetStrainPatterns(code, supcode):
	"""
	Given a code number for the crystal symmetry, 
	returns a list of strain patterns needed for
	the calculation of the elastic constants tensor.
	Each pattern is a 6 element list, the subscript
	reflects the strain in IRE notation

	Supported Strain Patterns
	-------------------------

	5 Cubic: e1+e4
	7 Hexagonal: e3 and e1+e4
	7 Trigonal-High (32, 3m, -3m): e1 and e3+e4
	6 Trigonal-Low (3, -3): e1 and e3+e4
	4 Tetragonal: e1+e4 and e3+e6
	3 Orthorhombic: e1+e4 and e2+e5 and e3+e6
	2 Monoclinic: e1+e4 and e3+e6 and e2 and e5
	1 Triclinic: e1 to e6 separately
	0 Unknown...
	"""

	if (code == 1):
		pattern = [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
		           [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
		           [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
		           [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
		           [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
		           [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]
	elif (code == 2):
		pattern = [[1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
		           [0.0, 0.0, 1.0, 0.0, 0.0, 1.0],
		           [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
		           [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]]
	elif (code == 3):
		pattern = [[1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
		           [0.0, 1.0, 0.0, 0.0, 1.0, 0.0],
		           [0.0, 0.0, 1.0, 0.0, 0.0, 1.0]]
	elif (code == 4):
		pattern = [[1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
		           [0.0, 0.0, 1.0, 0.0, 0.0, 1.0]]
	elif (code == 5):
		pattern = [[1.0, 0.0, 0.0, 1.0, 0.0, 0.0]]
	elif (code == 6):
		pattern = [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
		           [0.0, 0.0, 1.0, 1.0, 0.0, 0.0]]
	elif (code == 7):
		# NB is this correct for hex and hig trig? - see missmatch above/
		# I suspect I have to rotate lattice for trig high?
		pattern = [[0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
		           [1.0, 0.0, 0.0, 1.0, 0.0, 0.0]]
		if supcode == 1:
			pattern = [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
			           [0.0, 0.0, 1.0, 1.0, 0.0, 0.0]]

	return pattern	

def get_options(input_options, libmode):
	"""
	Just extracts the command line arguments into an options object
	"""
	if not libmode:
		p = optparse.OptionParser(usage="%prog [options] seedname\nGenerate CASTEP input for elastic constants calculation", \
		    version="%prog "+str(version))
		p.add_option('--debug', '-d', action='store_true', \
		              help="Debug mode (output to stdout rather than file)")
		p.add_option('--steps', '-n', action='store', type='int', dest="numsteps", \
	   	              help='Number of positive strain magnitudes to impose (defaults to 3)')
		p.add_option('--strain', '-s', action='store', type='float', dest="strain", \
		              help='Maximum magnitude of deformation to produced strained cells (defaults to 0.1)')
		p.add_option('--lattice', '-l',  action='store', type='int', dest="lattice", \
		              help='Lattice type to set pattern of deformation (extracted from .castep file)')
		options,arguments = p.parse_args(args=input_options)

	return options, arguments


def cellABC2cellCART (a, b, c, alp, bet, gam, Convention=1):
	"""
	Given three lattice vector lengths and angles, returns
	three vectors (as list of lists:
	[[a_x, a_y, a_z], [b_x, b_y, b_z], [c_x, c_y, c_z]]) representing
	the vectors on a cartisian frame. 
	"""
	# Get lattice vetors on cart frame from a, b, c and angles
	# For monoclinic, b // Y and c // Z
	if (alp == 90.0):
                sina = 1.0
      		cosa = 0.0
	else:
                sina = np.sin(np.radians(alp))
      		cosa = np.cos(np.radians(alp))
	if (bet == 90.0):
		sinb = 1.0
     		cosb = 0.0
	else:
		sinb = np.sin(np.radians(bet))
		cosb = np.cos(np.radians(bet))
	if (gam == 90.0):
		sing = 1.0
		cosg = 0.0
	else:
		sing = np.sin(np.radians(gam))
		cosg = np.cos(np.radians(gam))
        c_x = 0.0
        c_y = 0.0
        c_z = c

	b_x = 0.0
	b_y = b*sina
	b_z = b*cosa

        a_z = a*cosb
        a_y = a*(cosg - cosa*cosb)/sina
	trm1 = a_y/a
	a_x = a*np.sqrt(1.0 - cosb**2 - trm1**2)

	return [[a_x, a_y, a_z], [b_x, b_y, b_z], [c_x, c_y, c_z]]

def cellCART2cellABC (lattice):
	"""
	Given three latice vectors (with three componets each) return 
	the lattice vector lengths and angles between them. Input argument
	should be [[a_x, a_y, a_z], [b_x, b_y, bz], [c_x, c_y, c_z]]. Angles
	returned in degrees.
	"""
	# Does not care about orentation...
	a = np.sqrt(lattice[0][0]**2 + lattice[0][1]**2 + lattice[0][2]**2)
	b = np.sqrt(lattice[1][0]**2 + lattice[1][1]**2 + lattice[1][2]**2)
	c = np.sqrt(lattice[2][0]**2 + lattice[2][1]**2 + lattice[2][2]**2)
	gam = np.arccos(np.dot(lattice[0],lattice[1]) / (a*b))
	bet = np.arccos(np.dot(lattice[0],lattice[2]) / (a*c))
	alp = np.arccos(np.dot(lattice[1],lattice[2]) / (b*c))
	return a, b, c, np.degrees(alp), np.degrees(bet), np.degrees(gam)

			
def main(input_options, libmode=False):

	# deal with options
	options, arguments = get_options(input_options, libmode)
	seedname = arguments[0]
	
	(cell,pointgroup,atoms) = castep.parse_dotcastep(seedname)
	# Re-align lattice vectors on cartisian system
	a, b, c, al, be, ga = cellCART2cellABC(cell)
	cell = cellABC2cellCART(a, b, c, al, be, ga)


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
			# .castep: we are in trouble
			print "No point group found in .castep file so the strain pattern cannot be determined\n"
			print "A strain pattern can also be provided using the -l flag\n"
			sys.exit(1)
		else:
			# Use the value from .castep
			latticeCode, supcode = PointGroup2StrainPat(pointgroup)
	else:
		if (pointgroup == None):
			# Noting in .castep - use users choice
			latticeCode = options.lattice
		else:
			# Use users choice, but check and warn
			latticeCode = options.lattice
			if (latticeCode != PointGroup2StrainPat(pointgroup)[0]):
				print "WARNING: User supplied lattice code is inconsistant with the point group\n"
				print "         found by CASTEP. Using user supplied lattice code.\n"
		
	print "Cell parameters: a = %f gamma = %f" % (a, al)
	print "                 b = %f beta  = %f" % (b, be)
	print "                 c = %f gamma = %f \n" % (c, ga)
	print "Lattce vectors:  %7f %7f %7f " % (cell[0][0], cell[0][1], cell[0][2])
	print "                 %7f %7f %7f " % (cell[1][0], cell[1][1], cell[1][2])
	print "                 %7f %7f %7f \n " % (cell[2][0], cell[2][1], cell[2][2])
	patterns = GetStrainPatterns(latticeCode, supcode)
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
				castep.produce_dotcell(seedname, pattern_name+".cell", defcell, atoms)
				os.symlink(seedname+".param", pattern_name+".param")
	

	

if __name__ == "__main__":
	main(sys.argv[1:])
