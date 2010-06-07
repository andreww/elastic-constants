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



	# regular expression to match the whole of the final cell from a .castep file
	dotcastep_latt_RE = re.compile("""\sBFGS:\sFinal\sConfiguration:\s*\n
	            =+\s*\n\s*\n\s+\-+\s*\n\s+Unit\sCell\s*\n\s+\-+\s*\n
		    \s+Real\sLattice\(A\)\s+Reciprocal\sLattice\(1/A\)\s*\n
	            \s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s*\n
	            \s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s*\n
	            \s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s*\n""",
	          re.VERBOSE)
	# Start of the 'final configuration'
	dotcastep_infinal_RE = re.compile("BFGS: Final Configuration:")
	# Once inside final configuation, this should only match a line with atoms
	dotcastep_atomline_RE = re.compile("x\s+(\w+)\s+\d+\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+x")
	def parse_dotcastep(seedname):
		dotCastep = open(seedname+".castep","r")
		# Find the lattice
		latticeblock = dotcastep_latt_RE.findall(dotCastep.read())[-1] # Get the last block - handle concat restarts
		lattice = []
		lattice.append([float(latticeblock[0]), float(latticeblock[1]), float(latticeblock[2])])
		lattice.append([float(latticeblock[6]), float(latticeblock[7]), float(latticeblock[8])])
		lattice.append([float(latticeblock[12]), float(latticeblock[13]), float(latticeblock[14])])
		# rewind and search for and final atomic positions (these will be absent if e.g. they are all on symmetry positions)
		dotCastep.seek(0)
		in_atoms = False
		atoms = []
		for line in dotCastep:
			if (in_atoms and dotcastep_atomline_RE.search(line)):
				atom_line = dotcastep_atomline_RE.search(line)
				atoms.append([atom_line.group(1), float(atom_line.group(2)), float(atom_line.group(3)), float(atom_line.group(4))])
			elif ((not in_atoms) and (dotcastep_infinal_RE.search(line))):
				in_atoms = True
			
                dotCastep.close()
		
		return (lattice, atoms)

	# Regular expressions to match a lattice block in a castep .cell file. Note that these
	# can be of the form %block lattice_abc or %block lattice_cart and are case insensitive
	dotcell_lattice_start_RE = re.compile("^\s*%BLOCK\s+LATTICE_(?:CART|ABC)",re.IGNORECASE)
	dotcell_lattice_end_RE = re.compile("^\s*%ENDBLOCK\s+LATTICE_(?:CART|ABC)",re.IGNORECASE)
	dotcell_atoms_start_RE = re.compile("^\s*%BLOCK\s+POSITIONS_(?:FRAC|ABS)", re.IGNORECASE)
	dotcell_atoms_end_RE = re.compile("^\s*%ENDBLOCK\s+POSITIONS_(?:FRAC|ABS)", re.IGNORECASE)
	def produce_dotcell(seedname, filename, newlattice, atoms):
		"""
		produce_dotcell: reads <seedname>.cell (CASTEP cell file
		and writes a new .cell file to <filename> replacing the 
		lattice block with a new crystalographic lattice <newlattice> 
		(which should be supplied as a list of three lists, each with 
		three elements). Also adds command to fix cell during optimization.
		"""
		in_lattice = False
		in_atoms = False
		have_atoms = (atoms != []) # If we have an empty list, no atoms were optimized so just leave them in the .cell file.
		inputfile = open(seedname+".cell", "r")
		outputfile = open(filename, "w")
		for line in inputfile:
			if (dotcell_lattice_end_RE.search(line) and in_lattice):
				in_lattice = False
			elif (dotcell_lattice_start_RE.search(line) and not in_lattice):
				outputfile.write("%block LATTICE_CART\n")
				outputfile.write(str(defcell[0][0]) + " " + str(defcell[0][1]) + " " + str(defcell[0][2]) + "\n")
				outputfile.write(str(defcell[1][0]) + " " + str(defcell[1][1]) + " " + str(defcell[1][2]) + "\n")
				outputfile.write(str(defcell[2][0]) + " " + str(defcell[2][1]) + " " + str(defcell[2][2]) + "\n")
				outputfile.write("%endblock LATTICE_CART\n")
				outputfile.write("FIX_ALL_CELL true\n")
				in_lattice = True
			elif (dotcell_atoms_end_RE.search(line) and in_atoms and have_atoms):
				in_atoms = False
			elif ((dotcell_atoms_start_RE.search(line)) and (not in_atoms) and have_atoms):
				outputfile.write("%block POSITIONS_FRAC\n")
				for atom in atoms:
					outputfile.write("  " + atom[0] + "  " + str(atom[1]) + "  " + str(atom[2]) + "  " + str(atom[3]) + "\n")
				outputfile.write("%endblock POSITIONS_FRAC\n")
				in_atoms = True
			elif(not (in_lattice or in_atoms)):
				outputfile.write(line)
		inputfile.close
		outputfile.close
			
	options, arguments = get_options()
	seedname = arguments[0]
	
	(cell,atoms) = parse_dotcastep(seedname)

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
				produce_dotcell(seedname, pattern_name+".cell", defcell, atoms)
				os.symlink(seedname+".param", pattern_name+".param")
	

	

if __name__ == "__main__":
	main(sys.argv[1:])
