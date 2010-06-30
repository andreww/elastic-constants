#!/usr/bin/env python

import numpy as np
from rotCij import rotCij
from generate_strain import cellABC2cellCART

def _rotMat3(axis, theta):
	cos_theta = np.cos(theta)
	sin_theta = np.sin(theta)
	rotMat = np.zeros((3,3))
	rotMat[0,0] = cos_theta
	rotMat[1,1] = cos_theta
	rotMat[2,2] = cos_theta
	rotMat[axis,axis] = 1
	if axis == 0:
		axis_name = 'x'
		rotMat[2,1] = -1 * sin_theta
		rotMat[1,2] = sin_theta
	elif axis == 1:
		axis_name = 'y'
		rotMat[0,2] = -1 * sin_theta
		rotMat[2,0] = sin_theta
	elif axis == 2:
		axis_name = 'z'
		rotMat[1,0] = -1 * sin_theta
		rotMat[0,1] = sin_theta
	print "Rotation matrix for rotation of " + str(theta) + "radians around " + axis_name + "axis:"	
	print rotMat

	return rotMat

def rotVec(inVec, axis, theta):
	"""
	Rotates a 3-element vector aroound an 
	axis (0, 1 or 2) by an axis (in rads).
	Input and output is something that can 
	become a scipy matrix object. Does not
	modify input. Output is a np.array.
	"""
	rotMat = np.matrix(_rotMat3(axis, theta))
	outVec = np.transpose(rotMat * np.transpose(np.matrix(inVec)))
	outVec = [outVec[0,0], outVec[0,1], outVec[0,2]]
	return outVec

if __name__ == '__main__':
	import sys
	a = float(sys.argv[1])
	b = float(sys.argv[2])
	c = float(sys.argv[3])
	al = float(sys.argv[4])
	be = float(sys.argv[5])
	ga = float(sys.argv[6])
	inFile = file(sys.argv[7], 'r')	
	Cij_in = np.loadtxt(inFile)
	# Get cart vectors for lattice
	(lat_a, lat_b, lat_c) = cellABC2cellCART(a, b, c, al, be, ga)
	print "initial lattice vectors:"
	print lat_a
	print lat_b
	print lat_c
	# rotate by 90 - gamma around z to put b on y axis
	angle = -1.0 * np.radians(90 - ga) 
	lat_a = rotVec(lat_a, 2, angle)
	lat_b = rotVec(lat_b, 2, angle)
	lat_c = rotVec(lat_c, 2, angle)
	Cij = rotCij(Cij_in, 2, angle)
	print "after rotating around z:"
	print lat_a
	print lat_b
	print lat_c
	print Cij

	# calculate angle between projection
	# of c on x-z plane and z axis
	c_proj = np.cross(np.transpose([0.0,1.0,0.0]),np.cross(np.transpose(lat_c),np.transpose([0.0,1.0,0.0])))
	print c_proj
	angle = -1 * np.arccos(np.dot(np.transpose(c_proj), np.transpose([0.0, 0.0, 1.0]))/np.linalg.norm(c_proj))
	print angle
#	angle = np.degrees(np.arctan2(lat_c[0],lat_c[2]))
	#angle = angle - 0.003555
	

	# rotate by that angle around y to put c on yz plane
	lat_a = rotVec(lat_a, 1, angle)
	lat_b = rotVec(lat_b, 1, angle)
	lat_c = rotVec(lat_c, 1, angle)
	c_proj = rotVec(c_proj, 1, angle)
	Cij = rotCij(Cij, 1, angle)
	print "after rotating around y:"
	print c_proj
	print lat_a
	print lat_b
	print lat_c
	print Cij
	
