# Test the elastics analysis command line code using pytest
# this runs in a tempdir and uses numpy to check the calculated
# elastic constants are close to a reference set.
# Data comes from my diopside & jadeite paper

import os
import shutil
import subprocess

import numpy as np

def test_jadite_analysis(tmpdir):
    tmpdir = str(tmpdir) # We may run in python 2 with old os
    # copy input files
    print(tmpdir)
    print(os.listdir(tmpdir))
    shutil.copytree(os.path.dirname(os.path.realpath(__file__)) + "/jadite_files",tmpdir+"/jadeite_files")
    print(os.listdir(tmpdir))
    print(os.listdir(tmpdir+"/jadeite_files"))
    # Run the analysis code
    return_code = subprocess.call(['./elastics.py', tmpdir+"/jadeite_files/jadeite_4_30_GPa"])
    print(os.listdir(tmpdir+"/jadeite_files"))

    # check the results
    result = np.loadtxt(tmpdir+"/jadeite_files/jadeite_4_30_GPa_cij.txt")
    reference = np.array([
        [2.834355856208732121e+02, 9.333739462764282280e+01, 7.252504292812315612e+01,
         0.000000000000000000e+00, 4.431718910590255689e+00, 0.000000000000000000e+00],
        [9.333739462764282280e+01, 2.525910644867334156e+02, 8.182799549565147856e+01,
         0.000000000000000000e+00, 1.404182374948849699e+01, 0.000000000000000000e+00],
        [7.252504292812315612e+01, 8.182799549565147856e+01, 2.817874580455423370e+02,
         0.000000000000000000e+00, 4.149240710037339852e+01, 0.000000000000000000e+00],
        [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00,
         8.670817856178307181e+01, 0.000000000000000000e+00, 1.193825928636376688e+01],
        [4.431718910590255689e+00, 1.404182374948849699e+01, 4.149240710037339852e+01,
         0.000000000000000000e+00, 7.138439837622649975e+01, 0.000000000000000000e+00],
        [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00,
         1.193825928636376688e+01, 0.000000000000000000e+00, 9.655568727596958922e+01]])
    np.testing.assert_allclose(result, reference)

