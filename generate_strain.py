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

def PointGroup2StrainPat(pointgroup):
	"""
	Converts point group number (as ordered by CASTEP
	and in the International Tables) to a number representing 
	the needed strain pattern.
	"""
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
	elif (pointgroup <= 27):
		# Hexagonal
		patt = 7
	elif (pointgroup <= 32):
		# Cubic
		patt = 5
	else:
		print "Point group number " + str(pointgroup) + " not recognized.\n"
		sys.exit(1)
	return patt

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

def get_options(input_options, libmode):
	"""
	Just extracts the command line arguments into an options object
	"""
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
dotcastep_poinggroup_RE = re.compile("^\s+Point group of crystal =\s+([\+\-]?\d+):")

def parse_dotcastep(seedname):
	"""
	Extract lattice and atom positions from a .castep
	file. List of atoms may be empty (e.g. MgO)
	"""
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
	pointgroup = None
	atoms = []
	for line in dotCastep:
		sym_line = dotcastep_poinggroup_RE.search(line)
		atom_line = dotcastep_atomline_RE.search(line)
		if (in_atoms and atom_line):
			atoms.append([atom_line.group(1), float(atom_line.group(2)), \
			              float(atom_line.group(3)), float(atom_line.group(4))])
		elif ((not in_atoms) and (dotcastep_infinal_RE.search(line))):
			in_atoms = True
		elif (sym_line):
			pointgroup = int(sym_line.group(1))
		
	dotCastep.close()
	print lattice
	a, b, c, al, be, ga = cellCART2cellABC(lattice)
	print a, b, c, al, be, ga
	lattice = cellABC2cellCART(a, b, c, al, be, ga)
	print lattice
	
	return (lattice, pointgroup, atoms)

def cellABC2cellCART (a, b, c, alp, bet, gam):
	"""
	Given three lattice vector lengths and angles, returns
	three vectors (as list of lists:
	[[a_x, a_y, a_z], [b_x, b_y, b_z], [c_x, c_y, c_z]]) representing
	the vectors on a cartisian frame. 
	"""
	# Get lattice vetors on cart frame from a, b, c and angles
	# For monoclinic, b // Y and c // Z
	if (alp == 90.0):
      		cosa = 0.0;
	else:
      		cosa = np.cos(np.radians(alp))
	if (bet == 90.0):
     		cosb = 0.0
	else:
		cosb = np.cos(np.radians(bet))
	if (gam == 90.0):
		sing = 1.0
		cosg = 0.0
	else:
		sing = np.sin(np.radians(gam))
		cosg = np.cos(np.radians(gam))
	rv21 = 0.0
	rv31 = 0.0
	rv32 = 0.0
	rv11 = a
	rv12 = b*cosg
	rv22 = b*sing
	rv13 = c*cosb
	rv23 = c*(cosa - cosg*cosb)/sing
	trm1 = rv23/c
	rv33 = c*np.sqrt(1.0 - cosb**2 - trm1**2)
	return [[rv11, rv12, rv13], [rv21, rv22, rv23], [rv31, rv32, rv33]]

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


# Regular expressions to match a lattice block in a castep .cell file. Note that these
# can be of the form %block lattice_abc or %block lattice_cart and are case insensitive
dotcell_lattice_start_RE = re.compile("^\s*%BLOCK\s+LATTICE_(?:CART|ABC)",re.IGNORECASE)
dotcell_lattice_end_RE = re.compile("^\s*%ENDBLOCK\s+LATTICE_(?:CART|ABC)",re.IGNORECASE)
dotcell_atoms_start_RE = re.compile("^\s*%BLOCK\s+POSITIONS_(?:FRAC|ABS)", re.IGNORECASE)
dotcell_atoms_end_RE = re.compile("^\s*%ENDBLOCK\s+POSITIONS_(?:FRAC|ABS)", re.IGNORECASE)
def produce_dotcell(seedname, filename, defcell, atoms):
	"""
	produce_dotcell: reads <seedname>.cell (CASTEP cell file
	and writes a new .cell file to <filename> replacing the 
	lattice block with a new crystalographic lattice <defcell> 
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
	return()
			
def main(input_options, libmode=False):

	# deal with options
	options, arguments = get_options(input_options, libmode)
	seedname = arguments[0]
	
	(cell,pointgroup,atoms) = parse_dotcastep(seedname)

	cijdat = open(seedname+".cijdat","w")
	print "\nWriting strain data to ", seedname+".cijdat\n"

	# Not sure why the lattice types are enumerated like this, but this is how .cijdat does it...
	latticeTypes = {0:"Unknown", 1:"Triclinic", 2:"Monoclinic", 3:"Orthorhombic", \
	                4:"Tetragonal", 5:"Cubic", 6:"Trigonal-low", 7:"Trigonal-high/Hexagonal"}
	maxstrain = options.strain
	numsteps = options.numsteps
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
			latticeCode = PointGroup2StrainPat(pointgroup)
	else:
		if (pointgroup == None):
			# Noting in .castep - use users choice
			latticeCode = options.lattice
		else:
			# Use users choice, but check and warn
			latticeCode = options.lattice
			if (latticeCode != PointGroup2StrainPat(pointgroup)):
				print "WARNING: User supplied lattice code is inconsistant with the point group\n"
				print "         found by CASTEP. Using user supplied lattice code.\n"
		
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
				this_strain[0] = disps[0] / np.sqrt(cell[0][0]**2+cell[0][1]**2+cell[0][2]**2)
				this_strain[1] = disps[1] / np.sqrt(cell[1][0]**2+cell[1][1]**2+cell[1][2]**2)
				this_strain[2] = disps[2] / np.sqrt(cell[2][0]**2+cell[2][1]**2+cell[2][2]**2)
				# off diagonals - we only strain upper right corner of cell matrix, so strain is 1/2*du/dx...
				this_strain[3] = 0.5 * (disps[3] / np.sqrt(cell[0][0]**2+cell[0][1]**2+cell[0][2]**2))
				this_strain[4] = 0.5 * (disps[4] / np.sqrt(cell[1][0]**2+cell[1][1]**2+cell[1][2]**2))
				this_strain[5] = 0.5 * (disps[5] / np.sqrt(cell[2][0]**2+cell[2][1]**2+cell[2][2]**2))

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
