#!/usr/bin/env python
"""
castep.py

Various bits of python to read/write gulp files

Copyright (c) 2011 Andrew Walker (a.walker@ucl.ac.uk)
All rights reserved.
"""

import re

version = 0.1

# regular expression to match the whole of the final cell from a .castep file
dotgulp_latt_RE = re.compile("""\s+Final\sCartesian\slattice\svectors\s\(Angstroms\)\s:\s*\n
            \s*\n
            \s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s*\n
            \s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s*\n
            \s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s*\n""",
          re.VERBOSE)

# Start of the 'final configuration'
dotgulp_infinal_RE = re.compile("\s+Final asymmetric unit coordinates :")

# Once inside final configuation, this should only match a line with atoms
dotgulp_atomline_RE = re.compile("\s+\d+\s+(\w+)\s+(c|s|bs)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)")

dotgulp_atomsend_RE = re.compile("------------")

# Get the point group number
dotgulp_poinggroup_RE = re.compile("^\s+Point group of crystal =\s+([\+\-]?\d+):")

def parse_dotgot(seedname):
	"""
	Extract lattice and atom positions from a .castep
	file. List of atoms may be empty (e.g. MgO)
	"""
	dotGot = open(seedname+".got","r")
	# Find the lattice
	latticeblock = dotgulp_latt_RE.findall(dotGot.read())[-1] # Get the last block - handle concat restarts
	lattice = []
	lattice.append([float(latticeblock[0]), float(latticeblock[1]), float(latticeblock[2])])
	lattice.append([float(latticeblock[3]), float(latticeblock[4]), float(latticeblock[5])])
	lattice.append([float(latticeblock[6]), float(latticeblock[7]), float(latticeblock[8])])
	# rewind and search for and final atomic positions (these will be absent if e.g. they are all on symmetry positions)
	dotGot.seek(0)
	in_atoms = False
        at_sepcounter = 3
	pointgroup = None
	atoms = []
	for line in dotGot:
		sym_line = dotgulp_poinggroup_RE.search(line)
		atom_line = dotgulp_atomline_RE.search(line)
		if (in_atoms and atom_line):
			atoms.append([atom_line.group(1), atom_line.group(2), \
			          float(atom_line.group(3)), float(atom_line.group(4)), \
			          float(atom_line.group(5)), float(atom_line.group(6))])
		elif ((not in_atoms) and (dotgulp_infinal_RE.search(line))):
			in_atoms = True
                elif ((in_atoms) and (dotgulp_atomsend_RE.search(line))):
                    at_sepcounter = at_sepcounter - 1
                    if (at_sepcounter == 0): in_atoms = False
		elif (sym_line):
			pointgroup = int(sym_line.group(1))
		
	dotGot.close()
	return (lattice, pointgroup, atoms)


# Regular expressions to match a lattice block in a gulp .gin file. Note that this
# parser if far from foolproof.
dotgulp_lattice_start_RE = re.compile("\s*(?:CELL|VECT)",re.IGNORECASE)
dotgulp_lattice_end_RE = re.compile("\s*(?:\d+\.?\d*\s+\d+\.?\d*\s+\d+\.?\d*|(?:[01]\s+){6}|\d+\.?\d*\s+\d+\.?\d*\s+\d+\.?\d*\s+\d+\.?\d*\s+\d+\.?\d*\s+\d+\.?\d*\s*(?:[01]\s*){6})",re.IGNORECASE) # Note - this is a backwards match/
#dotgulp_atoms_start_RE = re.compile("\s*(?:FRAC|CART)", re.IGNORECASE)
dotgulp_atoms_start_RE = re.compile("fractional", re.IGNORECASE)
dotgulp_atoms_end_RE = re.compile("^\s*\w\w?\s+(?:core|shel|bshe)?\s*\d+\.?\d*", re.IGNORECASE) # Also nve search

def produce_dotgin(seedname, filename, defcell, atoms):
	"""
	produce_dotcell: reads <seedname>.gin (gulp input file
	and writes a new .gin file to <filename> replacing the 
	lattice block with a new crystalographic lattice <defcell> 
	(which should be supplied as a list of three lists, each with 
	three elements). Also adds keyword fix cell during optimization.
	"""
	done_keywords = False
        in_lattice = False
	in_atoms = False
	have_atoms = (atoms != []) # If we have an empty list, no atoms were optimized so just leave them in the .cell file.
	inputfile = open(seedname+".gin", "r")
	outputfile = open(filename, "w")
	for line in inputfile:
                # NB this is a serise of if ... continue blocks not a lot of 
                # if ... elif ... as we need to carry with the termnating line
                # for some cases
                if (re.match('\s*#',line)):
		    outputfile.write(line) # Just a comment - skip (may be before keywords)
                    continue

                if (not done_keywords):
                    #line = re.sub('conp', 'conv', line, 1, flags=re.IGNORECASE) # python 2.7
                    line = re.sub('conp', 'conv', line, 1)
	            done_keywords = True
		    outputfile.write(line)
                    continue

		if ((not dotgulp_lattice_end_RE.match(line)) and in_lattice):
			in_lattice = False

		if (dotgulp_lattice_start_RE.match(line) and not in_lattice):
			outputfile.write("vectors\n")
			outputfile.write(str(defcell[0][0]) + " " + str(defcell[0][1]) + " " + str(defcell[0][2]) + "\n")
			outputfile.write(str(defcell[1][0]) + " " + str(defcell[1][1]) + " " + str(defcell[1][2]) + "\n")
			outputfile.write(str(defcell[2][0]) + " " + str(defcell[2][1]) + " " + str(defcell[2][2]) + "\n")
			in_lattice = True
                        continue

		if ((not dotgulp_atoms_end_RE.match(line)) and in_atoms and have_atoms):
			in_atoms = False
  
		if ((dotgulp_atoms_start_RE.match(line)) and (not in_atoms) and have_atoms):
			outputfile.write("fractional\n")
			for atom in atoms:
                                if (atom[1] == 'c'):
                                    attype = 'core'
                                elif(atom[1] == 's'):
                                    attype = 'shel'
                                elif(atom[1] == 'bs'):
                                    attype = 'bshe'
                                else:
                                    print "Gulp parse error - wrong atom type"
				outputfile.write("  " + atom[0] + "  " + attype + "  " + str(atom[2]) + "  " + str(atom[3]) + "  " + str(atom[4]) + "\n")
			in_atoms = True
                        continue
		
		if(not (in_lattice or in_atoms)):
			outputfile.write(line)

	inputfile.close
	outputfile.close
	return()
			
