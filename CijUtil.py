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

def latexCij(Cij, eCij, outputfile, nt=False):

	f = open(outputfile,"w")
	for i in range(6):
		for j in range(i,6):
			if ((Cij[i,j] != 0.0) and (eCij[i,j] != 0.0)):
				if (nt):
					f.write("c$_{{{0}{1}}}$ & {2:5.1f}$\pm${3:3.1f}  \n".format(i+1,j+1,Cij[i,j],eCij[i,j]))
				else:
					f.write("c$_{{{0}{1}}}$ & {2:5.1f}$\pm${3:3.1f} \\\\ \n".format(i+1,j+1,Cij[i,j],eCij[i,j]))
	f.write(" & \\\\ \n")

	(vB, rB, vG, rG, hB, hG, evB, erB, evG, erG, ehB, ehG) = polyCij(Cij, eCij)
	if (nt):
		f.write("B^v & {0:5.1f}$\pm${1:3.1f}   \n".format(vB, evB))
		f.write("B^r & {0:5.1f}$\pm${1:3.1f}   \n".format(rB, erB))
		f.write("B^{{vrh}} & {0:5.1f}$\pm${1:3.1f}   \n".format(hB, ehB))
		f.write("G^v & {0:5.1f}$\pm${1:3.1f}    \n".format(vG, evG))
		f.write("G^r & {0:5.1f}$\pm${1:3.1f}  \n".format(rG, erG))
		f.write("G^{{vrh}} & {0:5.1f}$\pm${1:3.1f}   \n".format(hG, ehG))
		f.write(" & \\\\ \n")
		f.write("A$_U$ & {0:5.4f}  \n".format(uAniso(Cij)))
	else:
		f.write("B^v & {0:5.1f}$\pm${1:3.1f} \\\\ \n".format(vB, evB))
		f.write("B^r & {0:5.1f}$\pm${1:3.1f} \\\\ \n".format(rB, erB))
		f.write("B^{{vrh}} & {0:5.1f}$\pm${1:3.1f} \\\\ \n".format(hB, ehB))
		f.write("G^v & {0:5.1f}$\pm${1:3.1f} \\\\ \n".format(vG, evG))
		f.write("G^r & {0:5.1f}$\pm${1:3.1f} \\\\ \n".format(rG, erG))
		f.write("G^{{vrh}} & {0:5.1f}$\pm${1:3.1f} \\\\ \n".format(hG, ehG))
		f.write(" & \\\\ \n")
		f.write("A$_U$ & {0:5.4f} \\\\ \n".format(uAniso(Cij)))
	f.close()

	return None
	

def txtCij(Cij, filename):

	f = open(filename, "a")
	for i in range(6):
		for j in range(i,6):
			f.write("{0:5.1f} ".format(Cij[i,j]))
        f.write("\n")

def CijStability(Cij):
	"""Check that the elastic constants matrix is positive
	definite - i.e. that the structure is stable to small 
	strains. This is done by finding the eigenvalues by 
	diagonalization and checking that they are all positive.
	See Born & Huang, "Dynamical Theory of Crystal Lattices"
	(1954) page 141."""

	stable = False
	(eigenvalues, eigenvectors) = np.linalg.eig(Cij)
        if (np.amin(eigenvalues) > 0.0):
		stable = True
	else:
		print "Crystal not stable to small strains"
		print "(Cij not positive definite)"
		print "Eigenvalues: " + str(eigenvalues)

	return stable
    	

def invertCij(Cij, eCij):
        """Given a square matrix and a square matrix of the errors
        on each element, return the inverse of the matrix and the 
        propogated errors on the inverse.

        We use numpy for the inversion and eq.10 of Lefebvre, 
        Keeler, Sobie and White ('Propagation of errors for 
        matrix inversion' Nuclear Instruments and Methods in 
        Physics Research A 451 pp.520-528; 2000) to calculate 
        the errors. The errors can be reported directly as the 
	errors on the inverse matrix but to do useful further 
	propogation we need to report the covar matrix too.
        This is calculated from eq.9 and we then extract the 
	diagonal elements to get the errors (rather than implementing
	eq.10 too).

        Tested with the matrix:
                0.700(7) 0.200(2)
                0.400(4) 0.600(6)               
        which gives back the inverse and squared errors reported
        in Table 1 of the above reference.
        """

        # Assuming we have a rank 2 square array
        # of the same size for input array. 
        # FIXME: add check

        # Calculate the inverse using numpy
        Sij = np.linalg.inv(Cij)

	# Set up output arrays (init as zeros)
        eSij = np.zeros_like(eCij)
        array_size = eSij[0].size
	vcovSij = np.zeros((array_size,array_size,array_size,array_size),dtype=type(eSij))

	# Build covariance arrays (i.e COV(C^-1[a,b],S^-1[b,c] - a 4d array).
	# This is an implementation of eq.9 of Lefebvre et al.
        for a in range (array_size):
                for b in range (array_size):
			for c in range (array_size):
				for d in range (array_size):
                        		for i in range (array_size):
                                		for j in range (array_size):
                                        		vcovSij[a,b,c,d] = vcovSij[a,b,c,d] + \
							      ((Sij[a,i]*Sij[c,i]) * (eCij[i,j]**2) * (Sij[j,b]*Sij[j,d]))

	# Extrct the "diagonal" terms, which are
	# the errors on the elements of the inverse
	# and could also be calculated using eq.10
        for a in range (array_size):
		for b in range (array_size):                
			eSij[a,b] = np.sqrt(vcovSij[a,b,a,b])

        return (Sij, eSij, vcovSij)


def polyCij(Cij, eCij=np.zeros((6,6))):
	"""Returns voight-reuss-hill average of elastic constants tensor
           and propogated error given the 6*6 matrix of elastic constants
           and the 6*6 matrix of errors. The errors are optional. Assumes
           no co-varance between the errors on the elastic constants but 
           does include the co-varenance on the (calculated) compliance 
           matrix."""

	#FIXME: Need to check all the cov terms are correct in the addition 
        #       of errors.

	# Need compliances too:
	(sij, eSij, covSij) = invertCij(Cij, eCij)

	# These equations are valid for all crystal systems (only 9 of the 21 elastic constants 
	# ever have to be used, e.g. see Anderson theory of the Earth, pg. 122).

	voigtB = (1.0/9)*(Cij[0,0] + Cij[1,1] + Cij[2,2] ) \
	       + (2.0/9)*(Cij[0,1] + Cij[0,2] + Cij[1,2])

	evB = np.sqrt( (1.0/81)*(eCij[0,0]**2 + eCij[1,1]**2 + eCij[2,2]**2) \
	              +(2.0/81)*(eCij[0,1]**2 + eCij[0,2]**2 + eCij[1,2]**2) )

	reussB = 1.0/((sij[0,0]+sij[1,1]+sij[2,2]) + 2*(sij[0,1]+sij[0,2]+sij[1,2]))

	# Note that COV(X+Z,Y) = COV(X,Y)+COV(Z,Y) and COV(SUM(Xi),SUM(Yj)) = SUM(SUM(COV(Xi,Yj)
        # c.f. http://mathworld.wolfram.com/Covariance.html
	erB = (np.sqrt(    eSij[0,0]**2 + eSij[1,1]**2 + eSij[2,2]**2  \
	               + 4*eSij[0,1]**2 + 4*eSij[0,2]**2 + 4*eSij[1,2]**2  \
                       + 2*covSij[0,0,1,1] + 2*covSij[0,0,2,2] + 2*covSij[1,1,2,2] \
                       + 4*covSij[0,0,0,1] + 4*covSij[0,0,0,2] + 4*covSij[0,0,1,2] \
                       + 4*covSij[1,1,0,1] + 4*covSij[1,1,0,2] + 4*covSij[1,1,1,2] \
                       + 4*covSij[2,2,0,1] + 4*covSij[2,2,0,2] + 4*covSij[2,2,1,2] \
                       + 8*covSij[0,1,0,2] + 8*covSij[0,1,1,2] + 8*covSij[0,2,1,2] )) \
		* reussB**2

	voigtG = (1.0/15)*(Cij[0,0] + Cij[1,1] + Cij[2,2] - Cij[0,1] - Cij[0,2] - Cij[1,2]) \
	       + (1.0/5)*(Cij[3,3] + Cij[4,4] + Cij[5,5])

	evG = np.sqrt( (1.0/225)*(eCij[0,0]**2 + eCij[1,1]**2 + eCij[2,2]**2 + eCij[0,1]**2 + eCij[0,2]**2 + eCij[1,2]**2) \
	              +(1.0/25)*(eCij[3,3]**2 + eCij[4,4]**2 + eCij[5,5]**2) )

	reussG = 15.0/(4*(sij[0,0]+sij[1,1]+sij[2,2]) - 4*(sij[0,1]+sij[0,2]+sij[1,2]) + 3*(sij[3,3]+sij[4,4]+sij[5,5]))

	erG = np.sqrt( \
	              16*(eSij[0,0]**2 + eSij[1,1]**2 + eSij[2,2]**2) \
	            + 16*(eSij[0,1]**2 + eSij[0,2]**2 + eSij[1,2]**2) \
	            +  9*(eSij[3,3]**2 + eSij[4,4]**2 + eSij[5,5]**2) \
                    + 32*covSij[0,0,1,1] + 32*covSij[0,0,2,2] + 32*covSij[1,1,2,2] \
                    + 32*covSij[0,0,0,1] + 32*covSij[0,0,0,2] + 32*covSij[0,0,1,2] \
                    + 32*covSij[1,1,0,1] + 32*covSij[1,1,0,2] + 32*covSij[1,1,1,2] \
                    + 32*covSij[2,2,0,1] + 32*covSij[2,2,0,2] + 32*covSij[2,2,1,2] \
                    + 32*covSij[0,1,0,2] + 32*covSij[0,1,1,2] + 32*covSij[0,2,1,2] \
                    + 24*covSij[0,0,3,3] + 24*covSij[0,0,4,4] + 24*covSij[0,0,5,5] \
                    + 24*covSij[1,1,3,3] + 24*covSij[1,1,4,4] + 24*covSij[1,1,5,5] \
                    + 24*covSij[2,2,3,3] + 24*covSij[2,2,4,4] + 24*covSij[2,2,5,5] \
                    + 24*covSij[0,1,3,3] + 24*covSij[0,1,4,4] + 24*covSij[0,1,5,5] \
                    + 24*covSij[0,2,3,3] + 24*covSij[0,2,4,4] + 24*covSij[0,2,5,5] \
                    + 24*covSij[1,2,3,3] + 24*covSij[1,2,4,4] + 24*covSij[1,2,5,5] \
                    + 18*covSij[3,3,4,4] + 18*covSij[3,3,5,5] + 18*covSij[4,4,5,5] \
	             ) * (reussG**2 / 15)

	return (voigtB, reussB, voigtG, reussG, ((voigtB+reussB)/2.0), ((voigtG+reussG)/2.0), evB, erB, evG, erG, ((evB+erB)/2), ((evG+erG)/2))

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
	(voigtB, reussB, voigtG, reussG, hillB, hillG, evB, erB, evG, erG, ehB, ehG) = polyCij(Cij)
	return ((5*(voigtG/reussG))+(voigtB/reussB)-6) 

if __name__ == '__main__':
	import sys
	inFile = file(sys.argv[1], 'r')
	Cij_in = np.loadtxt(inFile)
	print "Input matrix:"
	print  np.array2string(Cij_in,max_line_width=130,suppress_small=True)
	(voigtB, reussB, voigtG, reussG, hillB, hillG, dum, dum, dum, dum, dum, dum) = polyCij(Cij_in)
        format = "%16s : %11.5f %11.5f %11.5f"
        print "\n                      Voigt       Reuss       Hill"
        print format % ("Bulk Modulus", voigtB, reussB, hillB)
        print format % ("Shear Modulus", voigtG, reussG, hillG)
        
