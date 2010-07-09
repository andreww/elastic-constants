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
import CijUtil

version = 1
	
def main(input_options, libmode=False):
	
		
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
				
		if options.graphics:
			try:
				global P
				import pylab as P
			except ImportError:
				print >> sys.stderr, "You need to have matplotlib installed for the --graphics option"
				sys.exit(1)
				
		return options, arguments
	
	options, arguments = get_options()

		
	# Not sure why the lattice types are enumerated like this, but this is how .cijdat does it...
	latticeTypes = {0:"Unknown", 1:"Triclinic", 2:"Monoclinic", 3:"Orthorhombic", \
		4:"Tetragonal", 5:"Cubic", 6:"Trigonal-low", 7:"Trigonal-high/Hexagonal"}

	# regular expression which matches the whole stress tensor block from a .castep file
	stressRE = re.compile("\s\*+\s(?:Symmetrised\s)?Stress\sTensor\s\*+\n.+\n.+?\((\w+)\).+\n.+\n.+\n.+\n\s\*\s+x\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+\*\n\s\*\s+y\s+[\+\-]?\d+.\d+\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+\*\n\s\*\s+z\s+[\+\-]?\d+.\d+\s+[\+\-]?\d+.\d+\s+([\+\-]?\d+.\d+)\s+\*\n")

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
	errors = S.zeros((21,1))
	
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

			if (S.__version__ < '0.7.0'):
				# correct for scipy weirdness - see http://www.scipy.org/scipy/scipy/ticket/8
				# This was fixed before 0.7.0 release. Maybe in some versions of 0.6.x too - 
				# will report huge errors if the check is wrong
				stderr = S.sqrt((numsteps * stderr**2)/(numsteps-2))
				error  = stderr/sqrt(sum(square(strain[:,index2-1])))
			else:
				# Work out the error ourselves as I cannot get it from
				# stderr and this has been checked with gnuplot's fitter
				fit_str = ((strain[:,index2-1] * cijFitted) + intercept)
				error = sqrt((sum(square(stress[:,index1-1] - fit_str)) / \
				             (numsteps-2))/(sum(square(strain[:,index2-1]))))
			
			# print info about the fit
			print '\n'
			print     'Cij (gradient)          :    ', cijFitted
			print     'Error in Cij            :    ', error
			print     'Intercept               :    ', intercept
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
			
			return cijFitted, error
		
			
		def __appendOrReplace(valList,erList,val):
			try:
				valList.append(val[0])
				erList.append(val[1])
				return (sum(valList)/len(valList)), (S.sqrt(sum([x**2 for x in erList])/len(erList)**2))
			except NameError:
				return val[0], val[1]
		
				
		def __createListAndAppend(val):
			newList = []
			newList.append(val[0])
			errorList = []
			errorList.append(val[1])
			return val[0], newList, val[1], errorList
		

		cij = S.zeros(21)
		
		# Analyse the patterns to see which strains were applied
		strainsUsed = analysePatterns(strain[0,:])

		# should check strains are as expected

		if symmetryType == "Cubic":
			
			if S.all(strainsUsed.transpose() == S.array([[1.0, 0.0, 0.0, 1.0, 0.0, 0.0]])):	# strain pattern e1+e4		

				finalCijs[0], errors[0] = __fit(1,1)                    # fit C11
				fit_21, fit_21_error = __fit(2,1)
				fit_31, fit_31_error = __fit(3,1)
				finalCijs[6] = (fit_21 + fit_31)/2   # fit C21+C31		
				errors[6] = S.sqrt((fit_21_error**2)/4 + (fit_31_error**2)/4)
				finalCijs[3], errors[3] = __fit(4,4)                    # fit C44
				
			else:
				print "Unsupported strain pattern"
				sys.exit(1)
				
		elif symmetryType == "Trigonal-high/Hexagonal":
			if S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]])):	# strain pattern e3 (hexagonal)

					# fit C13 + C23, and add to list (more values coming...)
					finalCijs[7], cij13, errors[7], er13 = __createListAndAppend(__fit(1,3))
					finalCijs[7], cij13, errors[7], er13 = __createListAndAppend(__fit(2,3))
										
					finalCijs[2], errors[2] = __fit(3,3)                # fit C33
					
			elif S.all(strainsUsed.transpose() == S.array([[1.0, 0.0, 0.0, 1.0, 0.0, 0.0]])):	# strain pattern e1+e4 (hexagonal)

					finalCijs[0], errors[0] = __fit(1,1)                          # fit C11
					finalCijs[6], errors[6] = __fit(2,1)                          # fit C21
					finalCijs[7], errors[7] = __appendOrReplace(cij13,er13,__fit(3,1)) # fit C31
					finalCijs[3], errors[3] = __fit(4,4)                          # fit C44

			elif S.all(strainsUsed.transpose() == S.array([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]])):	
				
					# strain pattern e1 (trigonal-high)

					finalCijs[0], errors[0] = __fit(1,1)                # fit C11
					finalCijs[6], errors[6] = __fit(2,1)                # fit C21
					finalCijs[7], errors[7] = __fit(3,1)                # fit C31
					finalCijs[8], errors[8] = __fit(4,1)                # fit C41
					finalCijs[9], errors[9] = __fit(5,1)                # fit C51
					
					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 1.0, 1.0, 0.0, 0.0]])):	
					
					# strain pattern e3+e4 (trigonal-high)
					# could recalculate C13/C14/C23/C24/C46 here, but won't just now
														
					finalCijs[2], errors[2] = __fit(3,3)                # fit C33
					finalCijs[3], errors[3] = __fit(4,4)                # fit C44
					
			else:
				print "Unsupported strain pattern"
				sys.exit(1)
				
		elif symmetryType == "Trigonal-low":
			if S.all(strainsUsed.transpose() == S.array([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]])):	
				
					# strain pattern e1 

					finalCijs[0], errors[0] = __fit(1,1)                # fit C11
					finalCijs[6], errors[6] = __fit(2,1)                # fit C21
					finalCijs[7], errors[7] = __fit(3,1)                # fit C31
					finalCijs[8], errors[8] = __fit(4,1)                # fit C41
					finalCijs[9], errors[9] = __fit(5,1)                # fit C51
					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 1.0, 1.0, 0.0, 0.0]])):	
				
					# strain pattern e3+e4
					# could recalculate C13/C14/C23/C24/C46 here, but won't just now

					finalCijs[2], errors[2] = __fit(3,3)                # fit C33
					finalCijs[3], errors[3] = __fit(4,4)                # fit C44
					
			else:
				print "Unsupported strain pattern"
				sys.exit(1)
		
		elif symmetryType == "Tetragonal":			
			if S.all(strainsUsed.transpose() == S.array([[1.0, 0.0, 0.0, 1.0, 0.0, 0.0]])):	# strain pattern e1+e4 

					finalCijs[0],  errors[0]  = __fit(1,1)               # fit C11
					finalCijs[6],  errors[6]  = __fit(2,1)               # fit C21
					finalCijs[7],  errors[7]  = __fit(3,1)               # fit C31
					finalCijs[10], errors[10] = __fit(6,1)               # fit C61
					finalCijs[3],  errors[3]  = __fit(4,4)               # fit C44

					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 1.0, 0.0, 0.0, 1.0]])):	# strain pattern e3+e6
					
					finalCijs[2], errors[2] = __fit(3,3)                # fit C33
					finalCijs[5], errors[5] = __fit(6,6)                # fit C66
					
			else:
				print "Unsupported strain pattern"
				sys.exit(1)
									
		elif symmetryType == "Orthorhombic":			
			if S.all(strainsUsed.transpose() == S.array([[1.0, 0.0, 0.0, 1.0, 0.0, 0.0]])):	# strain pattern e1+e4 

					finalCijs[0], errors[0] = __fit(1,1)                                # fit C11
					finalCijs[6], cij12, errors[6], er12 = __createListAndAppend(__fit(2,1))  # fit C21
					finalCijs[7], cij13, errors[7], er13 = __createListAndAppend(__fit(3,1))  # fit C31
					finalCijs[3], errors[3] = __fit(4,4)                                # fit C44
					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 1.0, 0.0, 0.0, 1.0, 0.0]])):	# strain pattern e2+e5 

					
					finalCijs[6], errors[6] = __appendOrReplace(cij12,er12,__fit(1,2))       # fit C12	
					finalCijs[1], errors[1] = __fit(2,2)                                # fit C22
					finalCijs[11], cij23, errors[11], er23 = __createListAndAppend(__fit(3,2)) # fit C32						
					finalCijs[4], errors[4] = __fit(5,5)                                # fit C55
					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 1.0, 0.0, 0.0, 1.0]])):	# strain pattern e3+e6 

					finalCijs[7],  errors[7]  = __appendOrReplace(cij13,er13,__fit(1,3))      # fit C13
					finalCijs[11], errors[11] = __appendOrReplace(cij23,er23,__fit(2,3))      # fit C23
					finalCijs[2], errors[2]  = __fit(3,3)                               # fit C33
					finalCijs[5], errors[5]  = __fit(6,6)                               # fit C66
					
			else:
				print "Unsupported strain pattern"
				sys.exit(1)
				
		elif symmetryType == "Monoclinic":			
			if S.all(strainsUsed.transpose() == S.array([[1.0, 0.0, 0.0, 1.0, 0.0, 0.0]])):	# strain pattern e1+e4 

					finalCijs[0], errors[0] = __fit(1,1)                                # fit C11
					finalCijs[6], cij12, errors[6], er12 = __createListAndAppend(__fit(2,1))  # fit C21
					finalCijs[7], cij13, errors[7], er13 = __createListAndAppend(__fit(3,1))  # fit C31
					finalCijs[3], errors[3] = __fit(4,4)                                # fit C44				
					finalCijs[9], cij51, errors[9], er51 = __createListAndAppend(__fit(5,1))  # fit C51	
					finalCijs[19], cij64, errors[19], er64 = __createListAndAppend(__fit(6,4)) # fit C64

					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 1.0, 0.0, 0.0, 1.0]])):	# strain pattern e3+e6 

					finalCijs[7], errors[7] = __appendOrReplace(cij13,er13,__fit(1,3))       # fit C13
					finalCijs[11], cij23, errors[11], er23 = __createListAndAppend(__fit(2,3)) # fit C23
					finalCijs[2], errors[2] = __fit(3,3)                                # fit C33
					finalCijs[16], cij53, errors[16], er53 = __createListAndAppend(__fit(5,3)) # fit C53
					finalCijs[19], errors[19] = __appendOrReplace(cij64,er64,__fit(4,6))      # fit C46
					finalCijs[5], errors[5] = __fit(6,6)                                # fit C66
					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 1.0, 0.0, 0.0, 0.0, 0.0]])):	# strain pattern e2

					finalCijs[6], errors[6]  = __appendOrReplace(cij12,er12,__fit(1,2))      # fit C12
					finalCijs[1], errors[1]  = __fit(2,2)                               # fit C22
					finalCijs[11],errors[11] = __appendOrReplace(cij23,er23,__fit(3,2))      # fit C32					
					finalCijs[13], cij52, errors[13], er52 = __createListAndAppend(__fit(5,2)) # fit C52

					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 0.0, 0.0, 1.0, 0.0]])):	# strain pattern e5

					finalCijs[9], errors[9]  = __appendOrReplace(cij51,er51,__fit(1,5))      # fit C15
					finalCijs[13],errors[13] = __appendOrReplace(cij52,er52,__fit(2,5))      # fit C25
					finalCijs[16],errors[16] = __appendOrReplace(cij53,er53,__fit(3,5))      # fit C35
					finalCijs[4], errors[4]  = __fit(5,5)                               # fit C55
			else:
				print "Unsupported strain pattern"
				sys.exit(1)
				
		elif symmetryType == "Triclinic":
			
			if S.all(strainsUsed.transpose() == S.array([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]])):	# strain pattern e1

					finalCijs[0], errors[0]  = __fit(1,1)                               # fit C11
					finalCijs[6], cij12, errors[6], er12 = __createListAndAppend(__fit(2,1))  # fit C21	
					finalCijs[7], cij13, errors[7], er13 = __createListAndAppend(__fit(3,1))  # fit C31					
					finalCijs[8], cij14, errors[8], er14 = __createListAndAppend(__fit(4,1))  # fit C41					
					finalCijs[9], cij15, errors[9], er15 = __createListAndAppend(__fit(5,1))  # fit C51					
					finalCijs[10],cij16, errors[10],er16 = __createListAndAppend(__fit(6,1))  # fit C61					
					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 1.0, 0.0, 0.0, 0.0, 0.0]])):	# strain pattern e2

					finalCijs[6], errors[6]  = __appendOrReplace(cij12,er12,__fit(1,2))       # fit C12
					finalCijs[1], errors[1]  = __fit(2,2)                                # fit C22
					finalCijs[11], cij23, errors[11], er23 = __createListAndAppend(__fit(3,2))  # fit C32	
					finalCijs[12], cij24, errors[12], er24 = __createListAndAppend(__fit(4,2))  # fit C42	
					finalCijs[13], cij25, errors[13], er25 = __createListAndAppend(__fit(5,2))  # fit C52	
					finalCijs[14], cij26, errors[14], er26 = __createListAndAppend(__fit(6,2))  # fit C62		
					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]])):	# strain pattern e3

					finalCijs[7],  errors[7]  = __appendOrReplace(cij13,er13,__fit(1,3))       # fit C13
					finalCijs[11], errors[11] = __appendOrReplace(cij23,er23,__fit(2,3))       # fit C23					
					finalCijs[2],  errors[2]  = __fit(3,3)                                # fit C33
					finalCijs[15], cij34, errors[15], er34 = __createListAndAppend(__fit(4,3))  # fit C43	
					finalCijs[16], cij35, errors[16], er35 = __createListAndAppend(__fit(5,3))  # fit C53	
					finalCijs[17], cij36, errors[17], er36 = __createListAndAppend(__fit(6,3))  # fit C63	
					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]])):	# strain pattern e4

					finalCijs[8],  errors[8]  = __appendOrReplace(cij14,er14,__fit(1,4))      # fit C14
					finalCijs[12], errors[12] = __appendOrReplace(cij24,er24,__fit(2,4))      # fit C24
					finalCijs[15], errors[15] = __appendOrReplace(cij34,er34,__fit(3,4))      # fit C34
					finalCijs[3],  errors[3]  = __fit(4,4)                               # fit C44
					finalCijs[18], cij45, errors[18], er45 = __createListAndAppend(__fit(5,4))  # fit C54	
					finalCijs[19], cij46, errors[19], er46 = __createListAndAppend(__fit(6,4))  # fit C64		

			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 0.0, 0.0, 1.0, 0.0]])):	# strain pattern e5
			
					finalCijs[9],  errors[9]  = __appendOrReplace(cij15,er15,__fit(1,5))     # fit C15
					finalCijs[13], errors[13] = __appendOrReplace(cij25,er25,__fit(2,5))     # fit C25
					finalCijs[16], errors[16] = __appendOrReplace(cij35,er35,__fit(3,5))     # fit C35
					finalCijs[18], errors[18] = __appendOrReplace(cij45,er45,__fit(4,5))     # fit C45
					finalCijs[4],  errors[4]  = __fit(5,5)                              # fit C55
					finalCijs[20], cij56, errors[20], er56 = __createListAndAppend(__fit(6,5))  # fit C65	
					
			elif S.all(strainsUsed.transpose() == S.array([[0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])):	# strain pattern e6
			
					finalCijs[10], errors[10]  = __appendOrReplace(cij16,er16,__fit(1,6))      # fit C16
					finalCijs[14], errors[14]  = __appendOrReplace(cij26,er26,__fit(2,6))      # fit C26
					finalCijs[17], errors[17]  = __appendOrReplace(cij36,er36,__fit(3,6))      # fit C36
					finalCijs[19], errors[19]  = __appendOrReplace(cij46,er46,__fit(4,6))      # fit C46
					finalCijs[20], errors[20]  = __appendOrReplace(cij56,er56,__fit(5,6))      # fit C56
					finalCijs[5],  errors[5]   = __fit(6,6)                               # fit C66		
					
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
		errors[5] = S.sqrt(0.25*(errors[0]**2+errors[6]**2))
	
	c = cMatrix(symmetryType,TetrHigh)
	
	# Generate the 6x6 matrix of elastic constants 
	# - negative values signify a symmetry relation
	finalCijMatrix = S.zeros((6,6))	
	finalErrors = S.zeros((6,6))
	for i in range(0,6):
		for j in range(0,6):
			index = int(c[i,j])
			if index > 0:	
				finalCijMatrix[i,j] = finalCijs[index-1]
				finalErrors[i,j] = errors[index-1]
			elif index < 0:
				finalCijMatrix[i,j] = -finalCijs[-index-1]
				finalErrors[i,j] = -errors[-index-1]
				
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
	print "\nErrors on Cij matrix ("+units+"):"
	print S.array2string(finalErrors,max_line_width=130,suppress_small=True)

	(sij, esij) = CijUtil.invertCij(finalCijMatrix,finalErrors)	
	
	print "\nFinal Sij matrix ("+units+"-1):"
	print S.array2string(sij,max_line_width=130,suppress_small=True)
	print "\nErrors on Sij matrix ("+units+"-1):"
	print S.array2string(esij,max_line_width=130,suppress_small=True)

	print"\n<>----------------------------------------------------------------------<>\n"	
	if symmetryType == "Cubic":
		print "  Zener anisotropy index     : %6.5f" % (CijUtil.zenerAniso(finalCijMatrix))
	print "  Universal anisotropy index : %6.5f" % (CijUtil.uAniso(finalCijMatrix))
	print "  (Rangnthn and Ostoja-Starzewski, PRL 101, 055504)\n"

	# bulkModulus = (finalCijs[1]+2*finalCijs[7])/3  #cubic only
	youngX = 1/sij[0,0]
	youngY = 1/sij[1,1]
	youngZ = 1/sij[2,2]
	
	
	format = "%18s : %11.5f %8s"
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
	(voigtB, reussB, voigtG, reussG, hillB, hillG) = CijUtil.polyCij(finalCijMatrix)
	format = "%16s : %11.5f %11.5f %11.5f %6s"
	print "                      Voigt       Reuss       Hill"
	print format % ("Bulk Modulus", voigtB, reussB, hillB, units)
	print format % ("Shear Modulus", voigtG, reussG, hillG, units)
	
	print "\n<>-----------------------------------------------------------------------<>\n"		
	
	S.savetxt(seedname + '_cij.txt', finalCijMatrix)	
	
	
def calculate(outfile, files, params=None, paramfile=None):
    return main([outfile,files], libmode=True)

if __name__ == '__main__':
	main(sys.argv[1:])

