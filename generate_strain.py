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
import math

def GetStrainPatterns(code):
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

	return pattern	

def main(input_options, libmode=False):

	# deal with options
	def get_options():
		if not libmode:
			p = optparse.OptionParser()
			p.add_option('--debug', '-d', action='store_true', \
			              help="Debug mode (output to stdout rather than file)")
			p.add_option('--num_steps', '-n', action='store', type='int', dest="numsteps")
			p.add_option('--strain_mag', '-s', action='store', type='float', dest="strain")
			p.add_option('--lattice_type', '-l', action='store', type='int', dest="lattice")
			options,arguments = p.parse_args(args=input_options)

		return options, arguments

	def parse_dotcastep(seedname):
		# regular expression to match the whole of the final cell from a .castep file
		m_BFGS_cell = re.compile("\sBFGS:\sFinal\sConfiguration:\s*\n=+\s*\n\s*\n\s+\-+\s*\n\s+Unit\sCell\s*\n\s+\-+\s*\n\s+Real\sLattice\(A\)\s+Reciprocal\sLattice\(1/A\)\s*\n\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s*\n\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s*\n\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s*\n")
		dotCastep = open(seedname+".castep","r")
		latticeblock = m_BFGS_cell.findall(dotCastep.read())[-1] # Get the last block - handle concat restarts
                dotCastep.close()
		lattice = []
		lattice.append([float(latticeblock[0]), float(latticeblock[1]), float(latticeblock[2])])
		lattice.append([float(latticeblock[6]), float(latticeblock[7]), float(latticeblock[8])])
		lattice.append([float(latticeblock[12]), float(latticeblock[13]), float(latticeblock[14])])
		
		return lattice

	def produce_dotcell(seedname, filename, newlattice):
		lattice_start_RE = re.compile("%BLOCK\sLATTICE_CART",re.IGNORECASE)
		lattice_end_RE = re.compile("%ENDBLOCK\sLATTICE_CART",re.IGNORECASE)
		in_lattice = False
		inputfile = open(seedname+".cell", "r")
		outputfile = open(filename, "w")
		for line in inputfile:
			if (lattice_end_RE.search(line) and in_lattice):
				in_lattice = False
			elif (lattice_start_RE.search(line) and not in_lattice):
				outputfile.write("%block LATTICE_CART\n")
				outputfile.write(str(defcell[0][0]) + " " + str(defcell[0][1]) + " " + str(defcell[0][2]) + "\n")
				outputfile.write(str(defcell[1][0]) + " " + str(defcell[1][1]) + " " + str(defcell[1][2]) + "\n")
				outputfile.write(str(defcell[2][0]) + " " + str(defcell[2][1]) + " " + str(defcell[2][2]) + "\n")
				outputfile.write("%endblock LATTICE_CART\n")
				outputfile.write("FIX_ALL_CELL true\n")
				in_lattice = True
			elif(not in_lattice):
				outputfile.write(line)
		inputfile.close
		outputfile.close
			
	options, arguments = get_options()
	seedname = arguments[0]
	
	cell = parse_dotcastep(seedname)

	cijdat = open(seedname+".cijdat","w")
	print "\nWriting strain data to ", seedname+".cijdat\n"

	# Not sure why the lattice types are enumerated like this, but this is how .cijdat does it...
	latticeTypes = {0:"Unknown", 1:"Triclinic", 2:"Monoclinic", 3:"Orthorhombic", \
	                4:"Tetragonal", 5:"Cubic", 6:"Trigonal-low", 7:"Trigonal-high/Hexagonal"}
	maxstrain = options.strain
	numsteps = options.numsteps
	latticeCode = options.lattice
	patterns = GetStrainPatterns(latticeCode)
	numStrainPatterns = len(patterns)
	print "\nLattice type is ", latticeTypes[latticeCode] +"\n"
	print "Number of patterns: "+ str(numStrainPatterns) +"\n"
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
				this_strain[0] = disps[0] / math.sqrt(cell[0][0]**2+cell[0][1]**2+cell[0][2]**2)
				this_strain[1] = disps[1] / math.sqrt(cell[1][0]**2+cell[1][1]**2+cell[1][2]**2)
				this_strain[2] = disps[2] / math.sqrt(cell[2][0]**2+cell[2][1]**2+cell[2][2]**2)
				# off diagonals - we only strain upper right corner of cell matrix, so strain is 1/2*du/dx...
				this_strain[3] = 0.5 * (disps[3] / math.sqrt(cell[0][0]**2+cell[0][1]**2+cell[0][2]**2))
				this_strain[4] = 0.5 * (disps[4] / math.sqrt(cell[1][0]**2+cell[1][1]**2+cell[1][2]**2))
				this_strain[5] = 0.5 * (disps[5] / math.sqrt(cell[2][0]**2+cell[2][1]**2+cell[2][2]**2))

				# Deform cell - only apply deformation to upper right corner
				defcell = [[cell[0][0]+disps[0], cell[0][1]+disps[5], cell[0][2]+disps[4]],
					   [cell[1][0],          cell[1][1]+disps[1], cell[1][2]+disps[3]],
					   [cell[2][0],          cell[2][1],          cell[2][2]+disps[2]]]

				pattern_name = seedname + "_cij__" + str(patt+1) + "__" + str((a*2)+1+neg)

				print "\nPattern Name = ", pattern_name
				print "Pattern = ", this_pat
				print "Magnitude = ", this_mag
				cijdat.write(pattern_name+"\n")
				cijdat.write(str(this_strain[0]) + " " + str(this_strain[5]) + " " + str(this_strain[4]) + "\n")
				cijdat.write(str(this_strain[5]) + " " + str(this_strain[1]) + " " + str(this_strain[3]) + "\n")
				cijdat.write(str(this_strain[4]) + " " + str(this_strain[3]) + " " + str(this_strain[2]) + "\n")
				produce_dotcell(seedname, pattern_name+".cell", defcell)
				os.symlink(seedname+".param", pattern_name+".param")
	

	

if __name__ == "__main__":
	main(sys.argv[1:])
