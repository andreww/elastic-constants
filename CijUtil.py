#!/usr/bin/env python
# encoding: utf-8
"""
CijUtil.py

Bits to help with elastic constants manipulation.

polyCij(Cij): Given an elastic constants matrix for a single crystal,
              calculate the Voight and Reuss bounds on the bulk and 
              shear moduli for a random polycrystal.

Copyright (c) 2010 Andrew Walker. All rights reserved.

"""

import numpy as np

version = 0.1


def polyCij(Cij):

	# Need compliances too:
	sij = np.linalg.inv(Cij)

	# These equations are valid for all crystal systems (only 9 of the 21 elastic constants 
	# ever have to be used, e.g. see Anderson theory of the Earth, pg. 122).

	voigtB = (1.0/9)*(Cij[0,0] + Cij[1,1] + Cij[2,2] ) \
	       + (2.0/9)*(Cij[0,1] + Cij[0,2] + Cij[1,2])

	reussB = 1.0/((sij[0,0]+sij[1,1]+sij[2,2]) + 2*(sij[0,1]+sij[0,2]+sij[1,2]))

	voigtG = (1.0/15)*(Cij[0,0] + Cij[1,1] + Cij[2,2] - Cij[0,1] - Cij[0,2] - Cij[1,2]) \
	       + (1.0/5)*(Cij[3,3] + Cij[4,4] + Cij[5,5])

	reussG = 15.0/(4*(sij[0,0]+sij[1,1]+sij[2,2]) - 4*(sij[0,1]+sij[0,2]+sij[1,2]) + 3*(sij[3,3]+sij[4,4]+sij[5,5]))

	return (voigtB, reussB, voigtG, reussG, ((voigtB+reussB)/2.0), ((voigtG+reussG)/2.0))

if __name__ == '__main__':
	import sys
	inFile = file(sys.argv[1], 'r')
	Cij_in = np.loadtxt(inFile)
	print "Input matrix:"
	print  np.array2string(Cij_in,max_line_width=130,suppress_small=True)
	(voigtB, reussB, voigtG, reussG, hillB, hillG) = polyCij(Cij_in)
        format = "%16s : %11.5f %11.5f %11.5f"
        print "\n                      Voigt       Reuss       Hill"
        print format % ("Bulk Modulus", voigtB, reussB, hillB)
        print format % ("Shear Modulus", voigtG, reussG, hillG)
