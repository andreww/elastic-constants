#!/usr/bin/env python
# encoding: utf-8
"""
rotCij.py

Tool to rotate the 6*6 matrix representation of an 
elastic constants tensor. 

Copyright (c) 2010 Andrew Walker. All rights reserved.

"""

import numpy as np

version = 0.1

def _rotMat6(axis, theta):
	# Obviously, this could be heavily optimized to reduce duplicate trig...
	# Rotation matrix is most zero
	# see http://www.engin.brown.edu/courses/en224/anis_general/anis_general.htm
	rotMat = np.zeros((6,6))
	# start of leading diag is cos**2(theta) other than one 1.
	rotMat[0,0] = (np.cos(theta))**2
	rotMat[1,1] = (np.cos(theta))**2
	rotMat[2,2] = (np.cos(theta))**2
	rotMat[axis,axis] = 1.0
	# end of leading diag is cos(theta) other than one (cos2 * sin2)
	rotMat[3,3] = np.cos(theta)
	rotMat[4,4] = np.cos(theta)
	rotMat[5,5] = np.cos(theta)
	rotMat[axis+3,axis+3] = (np.cos(theta))**2 - (np.sin(theta))**2
	# The rest as if's
	if axis == 0:
		rotMat[2,1] = (np.sin(theta))**2
		rotMat[1,2] = (np.sin(theta))**2
		rotMat[1,3] = 2.0*np.cos(theta)*np.sin(theta)
		rotMat[2,3] = -2.0*np.cos(theta)*np.sin(theta)
		rotMat[3,1] = -1.0*np.cos(theta)*np.sin(theta)
		rotMat[3,2] = np.cos(theta)*np.sin(theta)
		rotMat[4,5] = -1.0*np.sin(theta)
		rotMat[5,4] = np.sin(theta)
	elif axis == 1:
		rotMat[2,0] = (np.sin(theta))**2
		rotMat[0,2] = (np.sin(theta))**2
		rotMat[0,4] = 2.0*np.cos(theta)*np.sin(theta)
		rotMat[2,4] = -2.0*np.cos(theta)*np.sin(theta)
		rotMat[4,0] = -1.0*np.cos(theta)*np.sin(theta)
		rotMat[4,2] = np.cos(theta)*np.sin(theta)
		rotMat[3,5] = -1.0*np.sin(theta)
		rotMat[5,3] = np.sin(theta)
	elif axis == 2:
		rotMat[1,0] = (np.sin(theta))**2
		rotMat[0,1] = (np.sin(theta))**2
		rotMat[0,5] = 2.0*np.cos(theta)*np.sin(theta)
		rotMat[1,5] = -2.0*np.cos(theta)*np.sin(theta)
		rotMat[5,0] = -1.0*np.cos(theta)*np.sin(theta)
		rotMat[5,1] = np.cos(theta)*np.sin(theta)
		rotMat[4,3] = -1.0*np.sin(theta)
		rotMat[3,4] = np.sin(theta)

	return rotMat


def rotCij (inCij, axis, theta):
	"""
	Rotates an elastic constants tensor around
	an axis (0, 1 or 2) by an angle (in rads).
	Input and output tensor is a 6*6 scipy 
	array object. Does not change input matrix.
	"""
	rotMat = _rotMat6(axis, theta)
	rotMat_T = np.array(np.transpose(np.matrix(rotMat)))
	# Must be matrix for mat mult, or we'll do it by element.
	outCij = np.array(np.matrix(rotMat) * np.matrix(inCij) * np.matrix(rotMat_T))
	return outCij

if __name__ == '__main__':
	import sys
	# This is about as basic as it can be, but 
	# seems OK for testing.
	axis = int(sys.argv[1])
	theta = float(sys.argv[2])
	inFile = file(sys.argv[3], 'r')
	Cij_in = np.loadtxt(inFile)
	print "Rotating matrix by ", theta, " degrees around axis ", axis
	print "Input matrix: "
	print  np.array2string(Cij_in,max_line_width=130,suppress_small=True)
	Cij_out = rotCij(Cij_in, axis, np.radians(theta))
	print "Output matrix: "
	print np.array2string(Cij_out,max_line_width=130,suppress_small=True)
	if (len(sys.argv) > 4):
		np.savetxt(sys.argv[4], Cij_out)
		
	

