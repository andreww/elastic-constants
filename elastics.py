#!/usr/bin/env python
# encoding: utf-8
"""
elastics.py

Part of the MaterialsGrid project
Copyright (c) 2007-2008 Dan Wilson. All rights reserved.

"""

import sys
import os
import re
import optparse
import scipy as S

version = 1
	
def main(input_options, libmode=False):
	
	def dumpXML(seedname):
		
		import time
		
		
		S.set_printoptions(threshold=S.nan)  #don't skip printing parts of the array
		
		if os.path.isfile(seedname+"-geomopt1.cml"):
			# we'll inject the data into this file, then
			basefile = seedname+"-geomopt1.cml"
    		
			print "Inserting CML into: "+ basefile
			f = open(basefile, "r")
			tree = etree.parse(f)
			f.close()
			
			fp = d["{http://www.castep.org/cml/dictionary/}finalProperties"]
			finalproperties = fp.findin_context(tree)
			assert len(finalproperties) == 1
		    
			startnode = finalproperties[0]
		elif os.path.isfile(seedname+"-energy.cml"):
			# we'll inject the data into this file, then
			basefile = seedname+"-energy.cml"

			print "Inserting CML into: "+ basefile
			f = open(basefile, "r")
			tree = etree.parse(f)
			f.close()

			fp = d["{http://www.castep.org/cml/dictionary/}finalProperties"]
			finalproperties = fp.findin_context(tree)
			assert len(finalproperties) == 1

			startnode = finalproperties[0]
		else:
			
			
			basefile = seedname+"_cij_analysis.cml"	
			
			print "Outputting CML into: "+ basefile
			
			# create a header
			CML_NAMESPACE = "http://www.xml-cml.org/schema"
			DC_NAMESPACE ="http://purl.org/dc/elements/1.1/"		
			CML = "{%s}" % CML_NAMESPACE
			NSMAP = {None : CML_NAMESPACE, "dc" : DC_NAMESPACE}
				
			startnode = etree.Element(CML +"cml", nsmap=NSMAP)
		
			metadata = etree.Element("metadata")
			metadata.attrib["name"] = "dc:date"
			metadata.attrib["content"] = time.strftime("%Y-%m-%dT%H:%M:%S")
			startnode.append(metadata)
		
			metadata = etree.Element("metadata")
			metadata.attrib["name"] = "dc:creator"
			metadata.attrib["content"] = "MGElastics"
			startnode.append(metadata)
		
			metadata = etree.Element("metadata")
			metadata.attrib["name"] = "dc:hasVersion"
			metadata.attrib["content"] = str(version)
			startnode.append(metadata)
		
			metadata = etree.Element("metadata")
			metadata.attrib["name"] = "dc:subject"
			metadata.attrib["content"] = "Analysis file for Cij calculations"
			startnode.append(metadata)
			
			tree = etree.ElementTree(startnode)
			
			
		# now append Cij data
		cijnode = etree.Element("propertyList")
		cijnode.attrib["title"] = "MaterialsGrid: Elastic Properties"
		cijnode.attrib["dictRef"] = "castep:elastic_properties"  
		startnode.append(cijnode) 
		
		cijProp = etree.Element("property")
		cijProp.attrib["dictRef"] = "castep:elastic_stiffness_constants"
		cijProp.attrib["title"] = "Elastic Stiffness Constants"
		cijnode.append(cijProp)
		
		cijTensor = etree.Element("matrix")
		cijTensor.attrib["rows"] = "6"
		cijTensor.attrib["columns"] = "6"
		cijTensor.attrib["dataType"] = "xsd:double"
		cijTensor.attrib["units"] = "castepunits:"+ units.lower()
		cijTensor.text = S.array2string(finalCijMatrix.reshape(1,36)[0],max_line_width=1000,suppress_small=True).strip("[] ") 
		cijProp.append(cijTensor)
		
		sijProp = etree.Element("property")
		sijProp.attrib["dictRef"] = "castep:elastic_compliance_constants"
		sijProp.attrib["title"] = "Elastic Compliance Constants"
		cijnode.append(sijProp)
		
		sijTensor = etree.Element("matrix")
		sijTensor.attrib["rows"] = "6"
		sijTensor.attrib["columns"] = "6"
		sijTensor.attrib["dataType"] = "xsd:double"
		sijTensor.attrib["units"] = "castepunits:"+ units.lower()+"-1"  #this has to be changed - i don't think the '/' is allowed
		sijTensor.text = S.array2string(sij.reshape(1,36)[0],max_line_width=1000,suppress_small=True).strip("[] ") 
		sijProp.append(sijTensor)
		
		#### Young's Moduli ####
		youngsMod = etree.Element("propertyList")
		youngsMod.attrib["dictRef"] = "castep:youngs_moduli"
		youngsMod.attrib["title"] = "Young's Moduli"
		cijnode.append(youngsMod)
		
		# X
		youngsModValX = etree.Element("property")
		youngsModValX.attrib["title"] = "Young's Modulus X"
		youngsModValX.attrib["dictRef"] = "castep:young_x"
		youngsMod.append(youngsModValX)
		
		youngsModValXScalar = etree.Element("scalar")
		youngsModValXScalar.attrib["dataType"] = "xsd:double"
		youngsModValXScalar.attrib["units"] = "castepunits:"+ units.lower()
		youngsModValXScalar.text = str(youngX)
		youngsModValX.append(youngsModValXScalar)
		
		# Y
		youngsModValY = etree.Element("property")
		youngsModValY.attrib["title"] = "Young's Modulus Y"
		youngsModValY.attrib["dictRef"] = "castep:young_y"
		youngsMod.append(youngsModValY)
		
		youngsModValYScalar = etree.Element("scalar")
		youngsModValYScalar.attrib["dataType"] = "xsd:double"
		youngsModValYScalar.attrib["units"] = "castepunits:"+ units.lower()
		youngsModValYScalar.text = str(youngY)
		youngsModValY.append(youngsModValYScalar)
		
		# Z
		youngsModValZ = etree.Element("property")
		youngsModValZ.attrib["title"] = "Young's Modulus Z"
		youngsModValZ.attrib["dictRef"] = "castep:young_z"
		youngsMod.append(youngsModValZ)
		
		youngsModValZScalar = etree.Element("scalar")
		youngsModValZScalar.attrib["dataType"] = "xsd:double"
		youngsModValZScalar.attrib["units"] = "castepunits:"+ units.lower()
		youngsModValZScalar.text = str(youngZ)
		youngsModValZ.append(youngsModValZScalar)
		
		
		
		#### Poisson's Ratio ####
		poissonMod = etree.Element("propertyList")
		poissonMod.attrib["dictRef"] = "castep:poisson_ratio"
		poissonMod.attrib["title"] = "Poisson Ratio"
		cijnode.append(poissonMod)
		
		# XY
		poissonModValXY = etree.Element("property")
		poissonModValXY.attrib["title"] = "Poisson Ratio XY"
		poissonModValXY.attrib["dictRef"] = "castep:poisson_xy"
		poissonModValXY.attrib["units"] = "castepunits:dimensionless"
		
		poissonMod.append(poissonModValXY)
		
		poissonModValXYScalar = etree.Element("scalar")
		poissonModValXYScalar.attrib["dataType"] = "xsd:double"
		poissonModValXYScalar.attrib["units"] = "castepunits:dimensionless"
		poissonModValXYScalar.text = str(poissonXY)
		poissonModValXY.append(poissonModValXYScalar)
		
		# XZ
		poissonModValXZ = etree.Element("property")
		poissonModValXZ.attrib["title"] = "Poisson Ratio XZ"
		poissonModValXZ.attrib["dictRef"] = "castep:poisson_xz"
		poissonMod.append(poissonModValXZ)
		
		poissonModValXZScalar = etree.Element("scalar")
		poissonModValXZScalar.attrib["dataType"] = "xsd:double"
		poissonModValXZScalar.attrib["units"] = "castepunits:dimensionless"
		poissonModValXZScalar.text = str(poissonXZ)
		poissonModValXZ.append(poissonModValXZScalar)
		
		# YX
		poissonModValYX = etree.Element("property")
		poissonModValYX.attrib["title"] = "Poisson Ratio YX"
		poissonModValYX.attrib["dictRef"] = "castep:poisson_yx"
		poissonMod.append(poissonModValYX)
		
		poissonModValYXScalar = etree.Element("scalar")
		poissonModValYXScalar.attrib["dataType"] = "xsd:double"
		poissonModValYXScalar.attrib["units"] = "castepunits:dimensionless"
		poissonModValYXScalar.text = str(poissonYX)
		poissonModValYX.append(poissonModValYXScalar)
		
		# YZ
		poissonModValYZ = etree.Element("property")
		poissonModValYZ.attrib["title"] = "Poisson Ratio YZ"
		poissonModValYZ.attrib["dictRef"] = "castep:poisson_yz"
		poissonMod.append(poissonModValYZ)
		
		poissonModValYZScalar = etree.Element("scalar")
		poissonModValYZScalar.attrib["dataType"] = "xsd:double"
		poissonModValYZScalar.attrib["units"] = "castepunits:dimensionless"
		poissonModValYZScalar.text = str(poissonYZ)
		poissonModValYZ.append(poissonModValYZScalar)
		
		# ZX
		poissonModValZX = etree.Element("property")
		poissonModValZX.attrib["title"] = "Poisson Ratio ZX"
		poissonModValZX.attrib["dictRef"] = "castep:poisson_zx"
		poissonMod.append(poissonModValZX)
		
		poissonModValZXScalar = etree.Element("scalar")
		poissonModValZXScalar.attrib["dataType"] = "xsd:double"
		poissonModValZXScalar.attrib["units"] = "castepunits:dimensionless"
		poissonModValZXScalar.text = str(poissonZX)
		poissonModValZX.append(poissonModValZXScalar)
		
		# ZY
		poissonModValZY = etree.Element("property")
		poissonModValZY.attrib["title"] = "Poisson Ratio ZY"
		poissonModValZY.attrib["dictRef"] = "castep:poisson_zy"
		poissonMod.append(poissonModValZY)
		
		poissonModValZYScalar = etree.Element("scalar")
		poissonModValZYScalar.attrib["dataType"] = "xsd:double"
		poissonModValZYScalar.attrib["units"] = "castepunits:dimensionless"
		poissonModValZYScalar.text = str(poissonZY)
		poissonModValZY.append(poissonModValZYScalar)	
		
		
		#### Polycrystalline Results ####
		
		#bulk moduli
		bulkModPL = etree.Element("propertyList")
		bulkModPL.attrib["dictRef"] = "castep:polycrystalline_bulk_moduli"
		bulkModPL.attrib["title"] = "Polycrystalline Bulk Moduli"
		cijnode.append(bulkModPL)
		
		bulkModVoigt = etree.Element("property")
		bulkModVoigt.attrib["title"] = "Voigt"
		bulkModVoigt.attrib["dictRef"] = "castep:bulk_modulus_voigt"
		bulkModPL.append(bulkModVoigt)
		
		bulkModReuss = etree.Element("property")
		bulkModReuss.attrib["title"] = "Reuss"
		bulkModReuss.attrib["dictRef"] = "castep:bulk_modulus_reuss"
		bulkModPL.append(bulkModReuss)
		
		bulkModHill = etree.Element("property")
		bulkModHill.attrib["title"] = "Hill"
		bulkModHill.attrib["dictRef"] = "castep:bulk_modulus_hill"
		bulkModPL.append(bulkModHill)
		
		bulkModVoigtScalar = etree.Element("scalar")
		bulkModVoigtScalar.attrib["dataType"] = "xsd:double"
		bulkModVoigtScalar.attrib["units"] = "castepunits:"+ units.lower()
		bulkModVoigtScalar.text = str(voigtB)
		bulkModVoigt.append(bulkModVoigtScalar)
		
		bulkModReussScalar = etree.Element("scalar")
		bulkModReussScalar.attrib["dataType"] = "xsd:double"
		bulkModReussScalar.attrib["units"] = "castepunits:"+ units.lower()
		bulkModReussScalar.text = str(reussB)
		bulkModReuss.append(bulkModReussScalar)
		
		bulkModHillScalar = etree.Element("scalar")
		bulkModHillScalar.attrib["dataType"] = "xsd:double"
		bulkModHillScalar.attrib["units"] = "castepunits:"+ units.lower()
		bulkModHillScalar.text = str((voigtB+reussB)/2)
		bulkModHill.append(bulkModHillScalar)
		
		
		#shear moduli
		shearModPL = etree.Element("propertyList")
		shearModPL.attrib["dictRef"] = "castep:polycrystalline_shear_moduli"
		shearModPL.attrib["title"] = "Polycrystalline Shear Moduli"
		cijnode.append(shearModPL)
		
		shearModVoigt = etree.Element("property")
		shearModVoigt.attrib["title"] = "Voigt"
		shearModVoigt.attrib["dictRef"] = "castep:shear_modulus_voigt"
		shearModPL.append(shearModVoigt)
		
		shearModReuss = etree.Element("property")
		shearModReuss.attrib["title"] = "Reuss"
		shearModReuss.attrib["dictRef"] = "castep:shear_modulus_reuss"
		shearModPL.append(shearModReuss)
		
		shearModHill = etree.Element("property")
		shearModHill.attrib["title"] = "Hill"
		shearModHill.attrib["dictRef"] = "castep:shear_modulus_hill"
		shearModPL.append(shearModHill)
		
		shearModVoigtScalar = etree.Element("scalar")
		shearModVoigtScalar.attrib["dataType"] = "xsd:double"
		shearModVoigtScalar.attrib["units"] = "castepunits:"+ units.lower()
		shearModVoigtScalar.text = str(voigtG)
		shearModVoigt.append(shearModVoigtScalar)
		
		shearModReussScalar = etree.Element("scalar")
		shearModReussScalar.attrib["dataType"] = "xsd:double"
		shearModReussScalar.attrib["units"] = "castepunits:"+ units.lower()
		shearModReussScalar.text = str(reussG)
		shearModReuss.append(shearModReussScalar)
		
		shearModHillScalar = etree.Element("scalar")
		shearModHillScalar.attrib["dataType"] = "xsd:double"
		shearModHillScalar.attrib["units"] = "castepunits:"+ units.lower()
		shearModHillScalar.text = str((voigtG+reussG)/2)
		shearModHill.append(shearModHillScalar)
		
		if(options.debug):
			print etree.tostring(startnode,pretty_print=True)
		else:
			# wrap it in an ElementTree instance, and save as XML
			tree.write(basefile,pretty_print=True)
			
		return
		
	def analysePatterns(strain):

		# these are the IRE conventions, except that integers are 0->5 rather than 1->6
		strainDict = {0:"xx",1:"yy",2:"zz",3:"yz", 4:"zx", 5:"xy"}

		strainsUsed = S.zeros((6,1))

		for a in range(0,S.size(strain)):
			if strain[a] != 0.0:
				print strainDict[a], "component is non-zero"
				strainsUsed[a] = 1
			else:
				strainsUsed[a] = 0

		return strainsUsed
	

	def getCMLStress(file):

		def __castep(dict, term):
			return dict["{http://www.castep.org/cml/dictionary/}%s" % term]
		
            
		stress = __castep(d, "stress")

		tree = etree.parse(file)
		stress_xml = stress.findin(tree)
		assert len(stress_xml) == 1
		stressTensor = stress.getvalue(stress_xml[0])
		return stressTensor, stressTensor.unit
	

	def cMatrix(symmetryType,TetrHigh):
		if symmetryType == "Cubic" :
			return S.matrix([[1, 7, 7, 0, 0, 0],
						[7, 1, 7, 0, 0, 0],
						[7, 7, 1, 0, 0, 0],
						[0, 0, 0, 4, 0, 0],
						[0, 0, 0, 0, 4, 0],
						[0, 0, 0, 0, 0, 4]])

		elif symmetryType == "Trigonal-high/Hexagonal":
			return S.matrix([[1, 7, 8, 9, 10, 0],
						[7, 1, 8, 0,-9, 0],
						[8, 8, 3, 0, 0, 0],
						[9, -9, 0, 4, 0, 0],
						[10, 0, 0, 0, 4, 0],
						[0, 0, 0, 0, 0, 6]])

		elif symmetryType == "Trigonal-low":
			return S.matrix([[1, 7, 8, 9, 10, 0],
						[7, 1, 8, -9, -10, 0],
						[8, 8, 3, 0, 0, 0],
						[9, -9, 0, 4, 0, -10],
						[10,-10, 0, 0, 4, 9],
						[0, 0, 0, -10, 9, 6]])

		elif symmetryType == "Tetragonal":
			if TetrHigh == "-1":
				print "Higher-symmetry tetragonal (422,4mm,4-2m,4/mmm)"
				return S.matrix([[1, 7, 8, 0, 0, 0],
							[7, 1, 8, 0, 0, 0],
							[8, 8, 3, 0, 0, 0],
							[0, 0, 0, 4, 0, 0],
							[0, 0, 0, 0, 4, 0],
							[0, 0, 0, 0, 0, 6]])
			else:
				print "Lower-symmetry tetragonal (4,-4,4/m)"
				return S.matrix([[1, 7, 8, 0, 0, 11],
							[7, 1, 8, 0, 0, -11],
							[8, 8, 3, 0, 0, 0],
							[0, 0, 0, 4, 0, 0],
							[0, 0, 0, 0, 4, 0],
							[11, -11, 0, 0, 0, 6]])

		elif symmetryType == "Orthorhombic":
			return S.matrix([[ 1,  7,  8,  0,  0,  0],
						[ 7,  2, 12,  0,  0,  0],
						[ 8, 12,  3,  0,  0,  0],
						[ 0,  0,  0,  4,  0,  0],
						[ 0,  0,  0,  0,  5,  0],
						[ 0,  0,  0,  0,  0,  6]])

		elif symmetryType == "Monoclinic":
			return S.matrix([[ 1,  7,  8,  0,  10,  0],
						[ 7,  2, 12,  0, 14,  0],
						[ 8, 12,  3,  0, 17,  0],
						[ 0,  0,  0,  4,  0,  20],
						[10, 14, 17,  0,  5,  0],
						[ 0,  0,  0, 20,  0,  6]])

		elif symmetryType == "Triclinic":
			return S.matrix([[ 1,  7,  8,  9,  10, 11],
						[ 7,  2, 12,  13, 14,  15],
						[ 8, 12,  3,  16, 17,  18],
						[ 9, 13, 16,  4,  19,  20],
						[10, 14, 17, 19,  5,  21],
						[11, 15, 18, 20,  21,  6]])	
    

	def get_options():
		# deal with options
		if not libmode:
			p = optparse.OptionParser()
			p.add_option('--xml', '-x', action='store_true',  help="CML mode")
			p.add_option('--castep', '-c', action='store_true', help="CASTEP mode")
			p.add_option('--force-cml-output','-f', action='store_true', help="Force CML output",dest="force")
			p.add_option('--graphics', '-g', action='store_true', help="Show graphics (requires matplotlib)")
			p.add_option('--debug', '-d', action='store_true', help="Debug mode (output to stdout rather than file)")

			options,arguments = p.parse_args(args=input_options)
		else:
			class PotM_Options:
				xml = True
				graphics = True
				castep = False
				debug = False
			
			options = PotM_Options()
			
			taskRe = re.compile(r"(.+)-\w+\.cml")
			arguments = taskRe.findall(input_options[1][0])
			global outfile
			outfile = input_options[0]
		
		if options.castep and options.xml:
			p.error("options -x and -c are mutually exclusive")
		elif not options.castep and not options.xml:
			p.error("one of -x or -c is required")
		
		
		# For some reason, golem/lxml needs to be imported first    
		if options.xml or options.force:
			try:
				# import golem, for use globally
				global golem
				import golem
				
				# import etree, for use globally
				global etree
				from lxml import etree
				
				# set dictionary, for use globally
				global d
				d = golem.loadDictionary("castepDict.xml")
				
			except ImportError:
				print >> sys.stderr, "You need to have golem and lxml installed for the --xml option"
				sys.exit(1)
				
		if options.graphics:
			try:
				global P
				import pylab as P
			except ImportError:
				print >> sys.stderr, "You need to have matplotlib installed for the --graphics option"
				sys.exit(1)
				
		if(options.castep):
			print "Reading stress data from the .castep files"
		elif(options.xml):
			print "Reading stress data from the .xml files"
		else:
			print "Don't know where to get the stress data, please use flags --castep or --xml"
			sys.exit(1)
		
		return options, arguments
	
	options, arguments = get_options()

		
	# Not sure why the lattice types are enumerated like this, but this is how .cijdat does it...
	latticeTypes = {0:"Unknown", 1:"Triclinic", 2:"Monoclinic", 3:"Orthorhombic", \
		4:"Tetragonal", 5:"Cubic", 6:"Trigonal-low", 7:"Trigonal-high/Hexagonal"}

	# regular expression which matches the whole stress tensor block from a .castep file
	stressRE = re.compile("\s\*+\sSymmetrised\sStress\sTensor\s\*+\n.+\n.+?\((\w+)\).+\n.+\n.+\n.+\n\s\*\s+x\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+\*\n\s\*\s+y\s+[\+\-]?\d+.\d+\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+\*\n\s\*\s+z\s+[\+\-]?\d+.\d+\s+[\+\-]?\d+.\d+\s+([\+\-]?\d+.\d+)\s+\*\n")

	# Get strain tensors
	seedname = arguments[0]

	cijdat = open(seedname+".cijdat","r")
	print "\nReading strain data from ", seedname+".cijdat\n"

	numStrainPatterns = (len(cijdat.readlines())-2)/4 #total for all strain patterns

	#rewind
	cijdat.seek(0)

	# deal with those first four integers
	latticeType,numsteps,TetrHigh,TrigHigh = cijdat.readline().split()
	numsteps = int(numsteps)

	symmetryType = latticeTypes[int(latticeType)]
	print "System is", symmetryType,"\n"
	
	# get maximum magnitude of strains
	magnitude = float(cijdat.readline())
	print numsteps, "steps of maximum magnitude",magnitude
	
	# if using graphics, do some initial set-up
	if options.graphics:	
		fig = P.figure(num=1, figsize=(9.5,8),facecolor='white')
		fig.subplots_adjust(left=0.07,right=0.97,top=0.97,bottom=0.07,wspace=0.5,hspace=0.5)
		colourDict = {0: '#BAD0EF', 1:'#FFCECE', 2:'#BDF4CB', 3:'#EEF093',4:'#FFA4FF',5:'#75ECFD'}

		for index1 in range(6):
		    for index2 in range(6):
				# position this plot in a 6x6 grid
				sp = P.subplot(6,6,6*(index1)+index2+1)
				sp.set_axis_off()
				# change the labels on the axes
				# xlabels = sp.get_xticklabels()
				# P.setp(xlabels,'rotation',90,fontsize=7)
				# ylabels = sp.get_yticklabels()
				# P.setp(ylabels,fontsize=7)
				P.text(0.4,0.4, "n/a")
			
	print "\n<>---------------------------- ANALYSIS ---------------------------------<>"		
	
	# initialise 1d array to store all 21 unique elastic constants - will be transformed into 6x6 matrix later
	finalCijs = S.zeros((21,1))
	
	for patt in range(numStrainPatterns/numsteps):
		
		print "\nAnalysing pattern", patt+1, ":"
		
		for a in range(0,numsteps):  
		
			pattern = cijdat.readline()
	
			# grab the strain data from the .cijdat file
			line1 = cijdat.readline().split()
			line2 = cijdat.readline().split()
			line3 = cijdat.readline().split()
		
			# only take from the top right triangle
			# numbering according to IRE conventions (Proc IRE, 1949)
		
			if a == 0:
				strain = S.array([float(line1[0]),float(line2[1]),float(line3[2]),2*float(line2[2]),2*float(line1[2]),2*float(line1[1])])
			else:
				strain = S.row_stack((strain,S.array([float(line1[0]),float(line2[1]),float(line3[2]),2*float(line2[2]),2*float(line1[2]),2*float(line1[1])])))
		
			if(options.castep):
				# now get corresponding stress data from .castep
				dotCastep = open(seedname+"_cij__"+str(patt+1)+"__"+str(a+1)+".castep","r")	
				stressData = stressRE.findall(dotCastep.read())[0]
				dotCastep.close()
			
				units = stressData[0]
		
				# again, top right triangle
				if a == 0:
					stress = S.array([float(stressData[1]),float(stressData[4]),float(stressData[6]),float(stressData[5]),float(stressData[3]),float(stressData[2])])
				else:
					stress = S.row_stack((stress,S.array([float(stressData[1]),float(stressData[4]),float(stressData[6]),float(stressData[5]),float(stressData[3]),float(stressData[2])])))		
			
			elif(options.xml):
				# now get corresponding stress data from .xml (using Golem)
				dotXMLfilename = seedname+"_cij__"+str(patt+1)+"__"+str(a+1)+"-geomopt1.cml"
				cmlStressData, units = getCMLStress(dotXMLfilename)
			
				if a == 0:
					stress = S.array([float(cmlStressData[0][0]),float(cmlStressData[1][1]),float(cmlStressData[2][2]),float(cmlStressData[1][2]),float(cmlStressData[0][2]),float(cmlStressData[0][1])])
				else: 
					stress = S.row_stack((stress,S.array([float(cmlStressData[0][0]),float(cmlStressData[1][1]),float(cmlStressData[2][2]),float(cmlStressData[1][2]),float(cmlStressData[0][2]),float(cmlStressData[0][1])]) ))	
				

	

		"""
		Both the stress and strain matrices use the IRE conventions to reduce the
		3x3 matrices to 1x6 arrays. These 1D arrays are then stacked to form a 
		Nx6 array, where N=number of steps.
	
		Note that strain and stress arrays are numbered 0->5 rather than 1->6
		"""
		
		def __fit(index1, index2):
			from scipy import stats, sqrt, square
			
			# do the fit
			(cijFitted,intercept,r,tt,stderr) = stats.linregress(strain[:,index2-1],stress[:,index1-1])

			# correct for scipy weirdness - see http://www.scipy.org/scipy/scipy/ticket/8
			stderr = S.sqrt((numsteps * stderr**2)/(numsteps-2))
			error  = stderr/sqrt(sum(square(strain[:,index2-1])))
			
			# print info about the fit
			print '\n'
			print     'Cij (gradient)          :    ',cijFitted
			print     'Error in Cij            :    ', error
			if abs(r) > 0.9:
				print 'Correlation coefficient :    ',r
			else:
				print 'Correlation coefficient :    ',r, '     <----- WARNING'
			
			# if using graphics, add a subplot
			if options.graphics:
					
				# position this plot in a 6x6 grid
				sp = P.subplot(6,6,6*(index1-1)+index2)
				sp.set_axis_on()
				
				# change the labels on the axes
				xlabels = sp.get_xticklabels()
				P.setp(xlabels,'rotation',90,fontsize=7)
				ylabels = sp.get_yticklabels()
				P.setp(ylabels,fontsize=7)
			
				# colour the plot depending on the strain pattern
				sp.set_axis_bgcolor(colourDict[patt])

				# plot the data
				P.plot([strain[0,index2-1],strain[numsteps-1,index2-1]],[cijFitted*strain[0,index2-1]+intercept,cijFitted*strain[numsteps-1,index2-1]+intercept])
				P.plot(strain[:,index2-1],stress[:,index1-1],'ro')
			
			return cijFitted
		
			
		def __appendOrReplace(valList,val):
			try:
				valList.append(val)
				return sum(valList)/len(valList)
			except NameError:
				return val
		
				
		def __createListAndAppend(val):
			newList = []
			newList.append(val)
			return val, newList
		

		cij = S.zeros(21)
		
		# Analyse the patterns to see which strains were applied
		strainsUsed = analysePatterns(strain[0,:])

		# should check strains are as expected

		if symmetryType == "Cubic":
			
			if S.all(strainsUsed.transpose() == S.array([[1.0, 0.0, 0.0, 1.0, 0.0, 0.0]])):	# strain pattern e1+e4		

				finalCijs[0] = __fit(1,1)                    # fit C11
				finalCijs[6] = (__fit(2,1) + __fit(3,1))/2   # fit C21+C31		
				finalCijs[3] = __fit(4,4)                    # fit C44
				
			else:
				print "Unsupported strain pattern"
				sys.exit(1)
				
		elif symmetryType == "Trigonal-high/Hexagonal":
			if S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]])):	# strain pattern e3 (hexagonal)

					# fit C13 + C23, and add to list (more values coming...)
					cij13 = []
					cij13.append(__fit(1,3))
					cij13.append(__fit(2,3))
										
					finalCijs[2] = __fit(3,3)                # fit C33
					
			elif S.all(strainsUsed.transpose() == S.array([[1.0, 0.0, 0.0, 1.0, 0.0, 0.0]])):	# strain pattern e1+e4 (hexagonal)

					finalCijs[0] = __fit(1,1)                          # fit C11
					finalCijs[6] = __fit(2,1)                          # fit C21
					finalCijs[7] = __appendOrReplace(cij13,__fit(3,1)) # fit C31
					finalCijs[3] = __fit(4,4)                          # fit C44

			elif S.all(strainsUsed.transpose() == S.array([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]])):	
				
					# strain pattern e1 (trigonal-high)

					finalCijs[0] = __fit(1,1)                # fit C11
					finalCijs[6] = __fit(2,1)                # fit C21
					finalCijs[7] = __fit(3,1)                # fit C31
					finalCijs[8] = __fit(4,1)                # fit C41
					finalCijs[9] = __fit(5,1)                # fit C51
					
					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 1.0, 1.0, 0.0, 0.0]])):	
					
					# strain pattern e3+e4 (trigonal-high)
					# could recalculate C13/C14/C23/C24/C46 here, but won't just now
														
					finalCijs[2] = __fit(3,3)                # fit C33
					finalCijs[3] = __fit(4,4)                # fit C44
					
			else:
				print "Unsupported strain pattern"
				sys.exit(1)
				
		elif symmetryType == "Trigonal-low":
			if S.all(strainsUsed.transpose() == S.array([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]])):	
				
					# strain pattern e1 

					finalCijs[0] = __fit(1,1)                # fit C11
					finalCijs[6] = __fit(2,1)                # fit C21
					finalCijs[7] = __fit(3,1)                # fit C31
					finalCijs[8] = __fit(4,1)                # fit C41
					finalCijs[9] = __fit(5,1)                # fit C51
					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 1.0, 1.0, 0.0, 0.0]])):	
				
					# strain pattern e3+e4
					# could recalculate C13/C14/C23/C24/C46 here, but won't just now

					finalCijs[2] = __fit(3,3)                # fit C33
					finalCijs[3] = __fit(4,4)                # fit C44
					
			else:
				print "Unsupported strain pattern"
				sys.exit(1)
		
		elif symmetryType == "Tetragonal":			
			if S.all(strainsUsed.transpose() == S.array([[1.0, 0.0, 0.0, 1.0, 0.0, 0.0]])):	# strain pattern e1+e4 

					finalCijs[0]  = __fit(1,1)               # fit C11
					finalCijs[6]  = __fit(2,1)               # fit C21
					finalCijs[7]  = __fit(3,1)               # fit C31
					finalCijs[10] = __fit(6,1)               # fit C61
					finalCijs[3]  = __fit(4,4)               # fit C44

					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 1.0, 0.0, 0.0, 1.0]])):	# strain pattern e3+e6
					
					finalCijs[2] = __fit(3,3)                # fit C33
					finalCijs[5] = __fit(6,6)                # fit C66
					
			else:
				print "Unsupported strain pattern"
				sys.exit(1)
									
		elif symmetryType == "Orthorhombic":			
			if S.all(strainsUsed.transpose() == S.array([[1.0, 0.0, 0.0, 1.0, 0.0, 0.0]])):	# strain pattern e1+e4 

					finalCijs[0] = __fit(1,1)                                # fit C11
					finalCijs[6], cij12 = __createListAndAppend(__fit(2,1))  # fit C21
					finalCijs[7], cij13 = __createListAndAppend(__fit(3,1))  # fit C31
					finalCijs[3] = __fit(4,4)                                # fit C44
					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 1.0, 0.0, 0.0, 1.0, 0.0]])):	# strain pattern e2+e5 

					
					finalCijs[6] = __appendOrReplace(cij12,__fit(1,2))       # fit C12	
					finalCijs[1] = __fit(2,2)                                # fit C22
					finalCijs[11], cij23 = __createListAndAppend(__fit(3,2)) # fit C32						
					finalCijs[4] = __fit(5,5)                                # fit C55
					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 1.0, 0.0, 0.0, 1.0]])):	# strain pattern e3+e6 

					finalCijs[7]  = __appendOrReplace(cij13,__fit(1,3))      # fit C13
					finalCijs[11] = __appendOrReplace(cij23,__fit(2,3))      # fit C23
					finalCijs[2]  = __fit(3,3)                               # fit C33
					finalCijs[5]  = __fit(6,6)                               # fit C66
					
			else:
				print "Unsupported strain pattern"
				sys.exit(1)
				
		elif symmetryType == "Monoclinic":			
			if S.all(strainsUsed.transpose() == S.array([[1.0, 0.0, 0.0, 1.0, 0.0, 0.0]])):	# strain pattern e1+e4 

					finalCijs[0] = __fit(1,1)                                # fit C11
					finalCijs[6], cij12 = __createListAndAppend(__fit(2,1))  # fit C21
					finalCijs[7], cij13 = __createListAndAppend(__fit(3,1))  # fit C31
					finalCijs[3] = __fit(4,4)                                # fit C44				
					finalCijs[9], cij51 = __createListAndAppend(__fit(5,1))  # fit C51	
					finalCijs[19], cij64 = __createListAndAppend(__fit(6,4)) # fit C64

					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 1.0, 0.0, 0.0, 1.0]])):	# strain pattern e3+e6 

					finalCijs[7] = __appendOrReplace(cij13,__fit(1,3))       # fit C13
					finalCijs[11], cij23 = __createListAndAppend(__fit(2,3)) # fit C23
					finalCijs[2] = __fit(3,3)                                # fit C33
					finalCijs[16], cij53 = __createListAndAppend(__fit(5,3)) # fit C53
					finalCijs[19] = __appendOrReplace(cij64,__fit(4,6))      # fit C46
					finalCijs[5] = __fit(6,6)                                # fit C66
					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 1.0, 0.0, 0.0, 0.0, 0.0]])):	# strain pattern e2

					finalCijs[6]  = __appendOrReplace(cij12,__fit(1,2))      # fit C12
					finalCijs[1]  = __fit(2,2)                               # fit C22
					finalCijs[11] = __appendOrReplace(cij23,__fit(3,2))      # fit C32					
					finalCijs[13], cij52 = __createListAndAppend(__fit(5,2)) # fit C52

					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 0.0, 0.0, 1.0, 0.0]])):	# strain pattern e5

					finalCijs[9]  = __appendOrReplace(cij51,__fit(1,5))      # fit C15
					finalCijs[13] = __appendOrReplace(cij52,__fit(2,5))      # fit C25
					finalCijs[16] = __appendOrReplace(cij53,__fit(3,5))      # fit C35
					finalCijs[4]  = __fit(5,5)                               # fit C55
			else:
				print "Unsupported strain pattern"
				sys.exit(1)
				
		elif symmetryType == "Triclinic":
			
			if S.all(strainsUsed.transpose() == S.array([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]])):	# strain pattern e1

					finalCijs[0]  = __fit(1,1)                               # fit C11
					finalCijs[6], cij12 = __createListAndAppend(__fit(2,1))  # fit C21	
					finalCijs[7], cij13 = __createListAndAppend(__fit(3,1))  # fit C31					
					finalCijs[8], cij14 = __createListAndAppend(__fit(4,1))  # fit C41					
					finalCijs[9], cij15 = __createListAndAppend(__fit(5,1))  # fit C51					
					finalCijs[10],cij16 = __createListAndAppend(__fit(6,1))  # fit C61					
					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 1.0, 0.0, 0.0, 0.0, 0.0]])):	# strain pattern e2

					finalCijs[6]  = __appendOrReplace(cij12,__fit(1,2))       # fit C12
					finalCijs[1]  = __fit(2,2)                                # fit C22
					finalCijs[11], cij23 = __createListAndAppend(__fit(3,2))  # fit C32	
					finalCijs[12], cij24 = __createListAndAppend(__fit(4,2))  # fit C42	
					finalCijs[13], cij25 = __createListAndAppend(__fit(5,2))  # fit C52	
					finalCijs[14], cij26 = __createListAndAppend(__fit(6,2))  # fit C62		
					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]])):	# strain pattern e3

					finalCijs[7]  = __appendOrReplace(cij13,__fit(1,3))       # fit C13
					finalCijs[11] = __appendOrReplace(cij23,__fit(2,3))       # fit C23					
					finalCijs[2]  = __fit(3,3)                                # fit C33
					finalCijs[15], cij34 = __createListAndAppend(__fit(4,3))  # fit C43	
					finalCijs[16], cij35 = __createListAndAppend(__fit(5,3))  # fit C53	
					finalCijs[17], cij36 = __createListAndAppend(__fit(6,3))  # fit C63	
					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]])):	# strain pattern e4

					finalCijs[8]   = __appendOrReplace(cij14,__fit(1,4))      # fit C14
					finalCijs[12]  = __appendOrReplace(cij24,__fit(2,4))      # fit C24
					finalCijs[15]  = __appendOrReplace(cij34,__fit(3,4))      # fit C34
					finalCijs[3]   = __fit(4,4)                               # fit C44
					finalCijs[18], cij45 = __createListAndAppend(__fit(5,4))  # fit C54	
					finalCijs[19], cij46 = __createListAndAppend(__fit(6,4))  # fit C64		

			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 0.0, 0.0, 1.0, 0.0]])):	# strain pattern e5
			
					finalCijs[9]    = __appendOrReplace(cij15,__fit(1,5))     # fit C15
					finalCijs[13]   = __appendOrReplace(cij25,__fit(2,5))     # fit C25
					finalCijs[16]   = __appendOrReplace(cij35,__fit(3,5))     # fit C35
					finalCijs[18]   = __appendOrReplace(cij45,__fit(4,5))     # fit C45
					finalCijs[4]    = __fit(5,5)                              # fit C55
					finalCijs[20], cij56 = __createListAndAppend(__fit(6,5))  # fit C65	
					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])):	# strain pattern e6
			
					finalCijs[10]  = __appendOrReplace(cij16,__fit(1,6))      # fit C16
					finalCijs[14]  = __appendOrReplace(cij26,__fit(2,6))      # fit C26
					finalCijs[17]  = __appendOrReplace(cij36,__fit(3,6))      # fit C36
					finalCijs[19]  = __appendOrReplace(cij46,__fit(4,6))      # fit C46
					finalCijs[20]  = __appendOrReplace(cij56,__fit(5,6))      # fit C56
					finalCijs[5]   = __fit(6,6)                               # fit C66		
					
			else:
				print "Unsupported strain pattern"
				sys.exit(1)
		else:
			print "Unsupported symmetry type. Exiting"
			sys.exit(1)
	
	if options.graphics:
		P.savefig(os.path.basename(seedname)+'_fits')
		
	cijdat.close()
	
	if symmetryType == "Trigonal-high/Hexagonal" or symmetryType == "Trigonal-low":
		# for these systems, C66 is calculated as a combination of the other Cijs.
		finalCijs[5] = 0.5*(finalCijs[0]-finalCijs[6])
	
	c = cMatrix(symmetryType,TetrHigh)
	
	# Generate the 6x6 matrix of elastic constants 
	# - negative values signify a symmetry relation
	finalCijMatrix = S.zeros((6,6))	
	for i in range(0,6):
		for j in range(0,6):
			index = int(c[i,j])
			if index > 0:	
				finalCijMatrix[i,j] = finalCijs[index-1]
			elif index < 0:
				finalCijMatrix[i,j] = -finalCijs[-index-1]
				
	# Tests
	if symmetryType == "Cubic":
		if finalCijs[3] <= 0:
			print "\n *** WARNING: C44 is less than or equal to zero ***\n"
		if finalCijs[0] <= abs(finalCijs[6]):
			print "\n *** WARNING: C11 is less than or equal to |C12| ***\n"
		if (finalCijs[0]+2*finalCijs[6]) <= 0:
			print "\n *** WARNING: C11+2C12 is less than or equal to zero ***\n"
			
	
	print "\n<>---------------------------- RESULTS ----------------------------------<>\n"		
	print "Final Cij matrix ("+units+"):"
	print S.array2string(finalCijMatrix,max_line_width=130,suppress_small=True)
	
	
	sij = S.linalg.inv(finalCijMatrix)
	
	print "\nFinal Sij matrix ("+units+"-1):"
	print S.array2string(sij,max_line_width=130,suppress_small=True)
	
	# bulkModulus = (finalCijs[1]+2*finalCijs[7])/3  #cubic only
	youngX = 1/sij[0,0]
	youngY = 1/sij[1,1]
	youngZ = 1/sij[2,2]
	
	
	format = "%18s : %11.5f %8s"
	
	print ""
	# print format % ("Bulk Modulus", bulkModulus[0], units)
	# print format % ("Compressibility", 1/bulkModulus[0], "1/"+ units)
	
	print "\n                          x           y           z"
	print "%18s : %11.5f %11.5f %11.5f %6s" % ("Young's Modulus", youngX, youngY, youngZ, units)

	poissonXY = -1*sij[0,1]*youngX
	poissonXZ = -1*sij[0,2]*youngX
	poissonYX = -1*sij[1,0]*youngY
	poissonYZ = -1*sij[1,2]*youngY
	poissonZX = -1*sij[2,0]*youngZ
	poissonZY = -1*sij[2,1]*youngZ
	
	
	print "\n                        xy       xz       yx       yz       zx       zy"
	format = "%18s :  %6.5f  %6.5f  %6.5f  %6.5f  %6.5f  %6.5f"
	print format % ("Poisson's Ratios", poissonXY, poissonXZ, poissonYX, poissonYZ, poissonZX, poissonZY)
	
	
	print "\n<>--------------------- POLYCRYSTALLINE RESULTS -------------------------<>\n"		


	# These equations might be valid only for orthorhombic systems - check!
	
	voigtB = (1.0/9)*(finalCijMatrix[0,0] + finalCijMatrix[1,1] + finalCijMatrix[2,2] ) \
		+ (2.0/9)*(finalCijMatrix[0,1] + finalCijMatrix[0,2] + finalCijMatrix[1,2])
	
	reussB = 1.0/((sij[0,0]+sij[1,1]+sij[2,2]) + 2*(sij[0,1]+sij[0,2]+sij[1,2]))
	
	voigtG = (1.0/15)*(finalCijMatrix[0,0] + finalCijMatrix[1,1] + finalCijMatrix[2,2] - finalCijMatrix[0,1] - finalCijMatrix[0,2] - finalCijMatrix[1,2]) \
		+ (1.0/5)*(finalCijMatrix[3,3] + finalCijMatrix[4,4] + finalCijMatrix[5,5])
	
	reussG = 15.0/(4*(sij[0,0]+sij[1,1]+sij[2,2]) - 4*(sij[0,1]+sij[0,2]+sij[1,2]) + 3*(sij[3,3]+sij[4,4]+sij[5,5]))
	
	format = "%16s : %11.5f %11.5f %11.5f %6s"
	print "                      Voigt       Reuss       Hill"
	print format % ("Bulk Modulus", voigtB, reussB, (voigtB+reussB)/2, units)
	print format % ("Shear Modulus", voigtG, reussG, (voigtG+reussG)/2, units)
	
	print "\n<>-----------------------------------------------------------------------<>\n"		
	
	
	"""
	Here on in is just the XML output
	"""
	if(options.xml or options.force):
		dumpXML(seedname)

	# 
	
def calculate(outfile, files, params=None, paramfile=None):
    return main([outfile,files], libmode=True)

if __name__ == '__main__':
	main(sys.argv[1:])

