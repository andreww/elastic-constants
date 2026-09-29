#!/usr/bin/env python
# encoding: utf-8
"""
debye.py

Calculate (estimated) sound velocities and Debye temperatures according
to form proposed by Anderson [J. Phys. Chem. Solids 24, 909 (1963)].

Christopher D. Woodgate, 2024

"""

import numpy as np
import CijUtil

def AndersonDebye(rho, cell_vol, n_particles, Cij, eCij, nt=False):
    """Use calculated Hill bulk and shear moduli and return
       estimates of the polycrystalline longitudinal and transverse
       sound velocities and (consequently) the Debye temperature."""

    # TODO
    # 1. Think about error propogation

    # Planck constant
    h = 6.62607015e-34
    #Boltzmann constant
    k_B = 1.380649e-23

    (voigtB, reussB, voigtG, reussG, hillB, hillG, evB, erB, evG, erG, ehB, ehG) = CijUtil.polyCij(Cij, eCij)

    # Transverse and longitudinal sound velocities
    v_s = np.sqrt(hillG*1e9/rho)

    v_l = np.sqrt((hillB*1e9 + 4.0/3.0*hillG*1e9)/rho)

    # Averaged sound velocity according to
    # O. L. Anderson, J. Phys. Chem. Solids 24, 909-917 (1963)
    v_m = 1.0/np.cbrt(1.0/3.0*(2.0/v_s**3.0 + 1.0/v_l**3.0))

    # Anderson's formula (although not explicit in text)
    # works with the volume per ion of a material
    vol_per_ion = cell_vol/float(n_particles)

    # Anderson's formula for Debye temperature
    theta = h / k_B*v_m*np.cbrt((3.0)/(4.0*np.pi*vol_per_ion))

    return(v_s, v_l, v_m, theta)
    
