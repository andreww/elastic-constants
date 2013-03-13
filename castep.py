#!/usr/bin/env python
"""
castep.py

Various bits of python to read/write castep files

Copyright (c) 2010 Andrew Walker (a.walker@ucl.ac.uk)
All rights reserved.
"""

import re
import scipy as S

version = 0.1

# regular expression to match the whole of the final cell from a .castep file
dotcastep_latt_RE = re.compile("""\sBFGS\s*:\sFinal\sConfiguration:\s*\n
            =+\s*\n\s*\n\s+\-+\s*\n\s+Unit\sCell\s*\n\s+\-+\s*\n
	    \s+Real\sLattice\(A\)\s+Reciprocal\sLattice\(1/A\)\s*\n
            \s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s*\n
            \s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s*\n
            \s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s*\n""",
          re.VERBOSE)

# Start of the 'final configuration'
dotcastep_infinal_RE = re.compile("BFGS\s*: Final Configuration:")

# Once inside final configuation, this should only match a line with atoms
dotcastep_atomline_RE = re.compile("x\s+(\w+)\s+\d+\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+x")

# Get the point group number
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
	return (lattice, pointgroup, atoms)


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

# regular expression which matches the whole stress tensor block from a .castep file
stressRE = re.compile("\s\*+\s(?:Symmetrised\s)?Stress\sTensor\s\*+\n.+\n.+?\((\w+)\).+\n.+\n.+\n.+\n\s\*\s+x\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+\*\n\s\*\s+y\s+[\+\-]?\d+.\d+\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+\*\n\s\*\s+z\s+[\+\-]?\d+.\d+\s+[\+\-]?\d+.\d+\s+([\+\-]?\d+.\d+)\s+\*\n")

def get_stress_dotcastep(filename):
	"""Extract the stress tensor from a .castep file
	
	   Returns a tuple of (<units>, <stress>) where <units>
	   is a string representing the stress units and 
	   <stress> is a numpy vector of the elements of the 
	   stress tensor in the order s(1,1), s(2,2), s(3,3)
	   s(3,2), s(3,1), s(2,1).
	"""
	dotCastep = open(filename,"r")
	stressData = stressRE.findall(dotCastep.read())[0]
	dotCastep.close()
	units = stressData[0]
	stress = S.array([float(stressData[1]),float(stressData[4]),
                          float(stressData[6]),float(stressData[5]),
                          float(stressData[3]),float(stressData[2])])
	return(units, stress)

