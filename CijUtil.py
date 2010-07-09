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

def invertCij(Cij, eCij):
        """Given a square matrix and a square matrix of the errors
        on each element, return the inverse of the matrix and the 
        propogated errors on the inverse.

        We use numpy for the inversion and eq.10 of Lefebvre, 
        Keeler, Sobie and White ('Propagation of errors for 
        matrix inversion' Nuclear Instruments and Methods in 
        Physics Research A 451 pp.520-528; 2000) to calculate 
        the errors. We don't bother building the whole matrix 
        of coveriance matricies but we could from eq.9.

        Tested with the matrix:
                0.700(7) 0.200(2)
                0.400(4) 0.600(6)               
        which gives back the inverse and squared errors reported
        in Table 1 of the above reference (keeping in mind that 
        we only calculate the cov(eps[i,j],eps[i,j]) elements).
        """

        # Assuming we have a rank 2 square array
        # of the same size for input array. 
        # FIXME: add check

        # Need compliances 
        Sij = np.linalg.inv(Cij)

        eSij = np.zeros_like(eCij)
        array_size = eSij[0].size
        for a in range (array_size):
                for b in range (array_size):
                        for i in range (array_size):
                                for j in range (array_size):
                                        eSij[a,b] = eSij[a,b] + ((Sij[a,i]**2) * (eCij[i,j]**2) * (Sij[j,b]**2))
                        eSij[a,b] = np.sqrt(eSij[a,b])

        return (Sij, eSij)


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

def zenerAniso(Cij):
	"""Returns Zener anisotropy index, A, defined as
	2C44/(C11-C12). This is unity for an isotropic crystal 
	and, for a cubic crystal C44 and 1/2(C11-C12) are shear 
	strains accross the (100) and (110) planes, respectivly.
	See Zener, Elasticity and Anelasticity of Metals, 1948
	or doi:10.1103/PhysRevLett.101.055504 (c.f. uAniso).
	Note that we don't check that the crystal is cubic!"""
	return ((Cij[3,3]*2)/(Cij[0,0]-Cij[0,1]))

def uAniso(Cij):
	"""Returns the Universal elastic anisotropy index defined 
	by Ranganathan and Ostoja-Starzewski (PRL 101, 05504; 2008
	doi:10.1103/PhysRevLett.101.055504 ). Valid for all systems."""
	(voigtB, reussB, voigtG, reussG, hillB, hillG) = polyCij(Cij)
	return ((5*(voigtG/reussG))+(voigtB/reussB)-6)

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
