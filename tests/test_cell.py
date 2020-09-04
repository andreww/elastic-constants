#!/usr/bin/env python

import numpy as np
import numpy.testing as npt

import generate_strain as gs

def test_cell():
    # Should not end up rotates
    incell = [[9.1209595 ,  0.        , -2.65670807],
              [ 0.       ,  8.6       ,  0.   ],
              [ 0.       ,  0.        ,  5.2  ]]
    a, b, c, al, be, ga = gs.cellCART2cellABC (incell)
    outcell = gs.cellABC2cellCART (a, b, c, al, be, ga)
    npt.assert_almost_equal(incell, outcell)

def rot_cell(incell):
    a, b, c, al, be, ga = gs.cellCART2cellABC (incell)
    outcell = gs.cellABC2cellCART (a, b, c, al, be, ga)
    return outcell

def test_abc():
    # Should get out params back
    a = 9.5
    b = 8.6
    c = 5.2
    alp = 90.000
    bet = 106.239525
    gam = 90.000
    tmpcell = gs.cellABC2cellCART (a, b, c, alp, bet, gam)
    a_out = np.linalg.norm(tmpcell[0])
    b_out = np.linalg.norm(tmpcell[1])
    c_out = np.linalg.norm(tmpcell[2])
    alp_out = np.degrees(np.arccos(np.dot(tmpcell[1],tmpcell[2])/(np.linalg.norm(tmpcell[1])*np.linalg.norm(tmpcell[2]))))
    bet_out = np.degrees(np.arccos(np.dot(tmpcell[0],tmpcell[2])/(np.linalg.norm(tmpcell[0])*np.linalg.norm(tmpcell[2]))))
    gam_out = np.degrees(np.arccos(np.dot(tmpcell[0],tmpcell[1])/(np.linalg.norm(tmpcell[0])*np.linalg.norm(tmpcell[1]))))
    npt.assert_almost_equal(a_out,a)
    npt.assert_almost_equal(b_out,b)
    npt.assert_almost_equal(c_out,c)
    npt.assert_almost_equal(alp_out,alp)
    npt.assert_almost_equal(bet_out,bet)
    npt.assert_almost_equal(gam_out,gam)

def assert_cells_eq(cell_in, cell_out):
    a_out = np.linalg.norm(cell_out[0])
    b_out = np.linalg.norm(cell_out[1])
    c_out = np.linalg.norm(cell_out[2])
    alp_out = np.degrees(np.arccos(np.dot(cell_out[1],cell_out[2])/(np.linalg.norm(cell_out[1])*np.linalg.norm(cell_out[2]))))
    bet_out = np.degrees(np.arccos(np.dot(cell_out[0],cell_out[2])/(np.linalg.norm(cell_out[0])*np.linalg.norm(cell_out[2]))))
    gam_out = np.degrees(np.arccos(np.dot(cell_out[0],cell_out[1])/(np.linalg.norm(cell_out[0])*np.linalg.norm(cell_out[1]))))
    a = np.linalg.norm(cell_in[0])
    b = np.linalg.norm(cell_in[1])
    c = np.linalg.norm(cell_in[2])
    alp = np.degrees(np.arccos(np.dot(cell_in[1],cell_in[2])/(np.linalg.norm(cell_in[1])*np.linalg.norm(cell_in[2]))))
    bet = np.degrees(np.arccos(np.dot(cell_in[0],cell_in[2])/(np.linalg.norm(cell_in[0])*np.linalg.norm(cell_in[2]))))
    gam = np.degrees(np.arccos(np.dot(cell_in[0],cell_in[1])/(np.linalg.norm(cell_in[0])*np.linalg.norm(cell_in[1]))))
    npt.assert_almost_equal(a_out,a)
    npt.assert_almost_equal(b_out,b)
    npt.assert_almost_equal(c_out,c)
    npt.assert_almost_equal(alp_out,alp)
    npt.assert_almost_equal(bet_out,bet)

def test_jade():
    # Jadeite (the castep convention)
    jadeite_in = [[9.5000000,   0.0000000,   0.0000000],
                  [0.0000000,   8.6000000,   0.0000000],
                  [-1.4541981,   0.0000000,   4.9925252]]
    # Jadeite (c // z - the Cij convention)
    jadeite_out = [[9.1209595 ,  0.        , -2.65670807],
                   [ 0.       ,  8.6       ,  0.   ],
                   [ 0.       ,  0.        ,  5.2  ]]

    npt.assert_almost_equal(rot_cell(jadeite_in), jadeite_out)


def test_triclin():
    triclin_in = [[9.1209595 ,  0.        , 0.],
                  [ 0.5      ,  8.6       , 0 ],
                  [ 0.7      ,  0.2       ,  5.2]]
    triclin_out = [[ 9.0270942 ,  0.47424434,  1.21596251],
                   [ 0.        ,  8.60549714,  0.39423208],
                   [ 0.        ,  0.        ,  5.25071424]]
    assert_cells_eq(rot_cell(triclin_in), triclin_in)
    npt.assert_almost_equal(rot_cell(triclin_in), triclin_out)
    assert_cells_eq(triclin_in, triclin_out)

